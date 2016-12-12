#!/usr/bin/env python

import os
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib
# Plotting without X
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
from sklearn.cluster.bicluster import SpectralCoclustering
from sklearn.cluster.bicluster import SpectralBiclustering
import scipy.sparse as sparse
import ClusteringUtils as cu

pipelinePath = 'FortunatoPipelinePath.txt'
# pipelinePath = 'co-cluster-test-path.txt'

def invert_permutation(p):
    '''Returns an array s, where s[i] gives the index of i in p.
    The argument p is assumed to be some permutation of 0, 1, ..., len(p)-1. 
    '''
    n = len(p)
    s = np.zeros(n, dtype = np.int32)
    i = np.arange(n, dtype = np.int32)
    np.put(s, p, i) # s[p[i]] = i 
    return s

def build_paths():
    '''Build list of paths for tp files'''

    # Path to parent folder holding all tp files
    with open(pipelinePath, 'r') as fh:
        fortunatoPath = fh.readline().rstrip()

    # Regex to get file and capture cluster number/slice number
    f_re = re.compile('.*slice_out_oslom/slice_(\d+)/results_(\d+)/tp$')

    # Get tp paths, their slice number and their cluster number
    paths = {}
    for root, dirs, files in os.walk(fortunatoPath):
        for f in files:
            match = f_re.match(os.path.join(root,f))
            if match:
                p = match.group(0) #Path
                s = match.group(1) #Slice
                c = match.group(2) #Cluster

                if s in paths.keys():
                    paths[s].append({'path':p, 'cluster':c})
                else:
                    paths[s] = [{'path':p, 'cluster':c}]
    return paths

def build_matrix(args):
    '''Build matrix from tp column vectors'''

    M = None

    for s in args.slices:
        nResults = 0 # Track number of results
        try:
            for p in args.paths[s]:
                # Create matrix
                if M is None:
                    M, CL = tp_to_col(p['path'], args, [])
                # Append columns to matrix
                else:
                    C, CL = tp_to_col(p['path'], args, CL)
                    M = append_cols(C, M)

                # Break loop if a maximum number of results has been set
                if args.results != 'All':
                    nResults += 1
                    if nResults == int(args.results):
                        break
        except KeyError, e:
            print '%s not a valid slice'%e

    return M, int(np.median(map(int, CL)))

def build_matrix2(args):
    all_clusterings = [cu.readComms(p['path'])[1] for s in args.slices for p in args.paths[s]]
    all_clusters = [c for clustering in all_clusterings for c in clustering]
    all_nodes = set([n for c in all_clusters for n in c])
    ind2node = {i : node for (i, node) in enumerate(all_nodes)}
    node2ind = {node : i for (i, node) in ind2node.items()}

    M = sparse.dok_matrix((len(all_nodes), len(all_clusters)))
    for i in range(len(all_clusters)):
        for n in all_clusters[i]:
            M[node2ind[n], i] = 1

    median_size = np.median([len(c) for c in all_clusterings])
    args.ind2node = ind2node
    args.node2ind = node2ind
    print('Matrix built')
    return M, median_size

def append_cols(N,M):
    '''Append column to matrix. If their sizes differ, buffer by zeros'''

    # Matrix/column dimensions
    shapeM = M.shape
    shapeN = N.shape

    # Test if there is a difference in number of rows
    diff = shapeM[0] - shapeN[0]

    # Append column to matrix as new column
    if diff == 0:
        M = np.hstack([M, N])
    # Add extra zero-rows if new column has less rows than matrix.
    #  Append to matrix as new column
    elif diff > 0:
        extraRows = np.zeros([diff, shapeN[1]])
        N = np.vstack([N, extraRows])
        M = np.hstack([M, N])
    # Add extra zero-rows if matrix has less rows than new column.
    #  Append column to matrix as new column
    elif diff < 0:
        extraRows = np.zeros([abs(diff), shapeM[1]])
        M = np.vstack([M, extraRows])
        M = np.hstack([M, N])

    return M

def tp_to_col(p, args, clusterLengths):
    '''
    Read tp file column vector. Module information
    is stored in case it becomes useful later.
    '''

    # Dict of modules for tp file
    modules = {}

    # Regex to capture module information
    s = re.compile('module\s+(\d+)\s+size:\s+(\d+)\s+bs:\s+(.*)')

    # Read tp file into dict holding cluster information
    with open(p, 'r') as fh:
        line = fh.readline()
        max_val = 0

        while line:
            # Capture module information
            groups = s.search(line)
            members = map(int,fh.readline().rstrip().split(' '))

            # Record max value for determining number of rows
            max_val = max(max_val, max(members))

            modules[int(groups.group(1))] = {
                'size':groups.group(2),
                'bs':groups.group(3),
                'nodes':members}

            line = fh.readline()

    if args.set_module is not None:
        return modules
    else:
        # Create column vector
        nClusters = len(modules.keys())
        M = np.zeros([max_val, nClusters], dtype=float)
        for k in modules.keys():
            for i in modules[k]['nodes']:
                M[i - 1,k] = 1
        clusterLengths.append(nClusters)
        return M, clusterLengths

def Spectral_CoClustering(args):
    '''Function to perform bipartite clustering'''
    # Create model
    try:
        if args.arpack:
            model = SpectralCoclustering(
                n_clusters=args.nClusters, svd_method='arpack')
        else:
            model = SpectralCoclustering(
                n_clusters=args.nClusters)
    except:
        print '-r 1 may cause problems when svd_method has been set to arpack'
    print('Running coclustering')
    model.fit(args.M.tocsc())
    print('Coclustering done')

    # Fit to data
    # fit_data = args.M[np.argsort(model.row_labels_)]
    # fit_data = fit_data[:, np.argsort(model.column_labels_)]
    fit_data = args.M.tocoo()
    fit_data.row = invert_permutation(np.argsort(model.row_labels_))[fit_data.row]
    fit_data.col = invert_permutation(np.argsort(model.column_labels_))[fit_data.col]

    save_clusters(model, fit_data, args, '_CoClustering')

    return model, fit_data

def Spectral_BiClustering(M, args):
    '''Function to perform bipartite clustering'''
    # Create model
    try:
        if args.arpack:
            model = SpectralBiclustering(
                n_clusters=args.nClusters, svd_method='arpack')
        else:
            model = SpectralBiclustering(
                n_clusters=args.nClusters)
    except:
        print '-r 1 may cause problems when svd_method has been set to arpack'

    print('Running biclustering')
    model.fit(M.tocsc())
    print('Biclustering done')

    # Fit to data
    # fit_data = M[np.argsort(model.row_labels_)]
    # fit_data = fit_data[:, np.argsort(model.column_labels_)]
    fit_data = M.tocoo()
    fit_data.row = invert_permutation(np.argsort(model.row_labels_))[fit_data.row]
    fit_data.col = invert_permutation(np.argsort(model.column_labels_))[fit_data.col]

    save_clusters(model, fit_data, args, '_BiClustering', 1)

    return model, fit_data

def save_clusters(model, fit_data, args, suffix, idx=0):
    # Save fitted data if specified
    print('Saving')
    if args.fout is not None:
        # print(np.argsort(model.row_labels_))
        # coo = args.M.tocoo()
        # pdx = pd.DataFrame(
        #         {'data' : coo,
        #     'index' : np.argsort(model.row_labels_)})
            # 'columns' : np.argsort(model.column_labels_)})
        # Save entire dataframe to file. Commented out for speed
        # pdx.to_csv(args.df + '/' + args.fout + suffix + '.csv', index=True)

        # Save tp file
        with open(args.df + '/' + args.fout + suffix + '_tp.csv', 'w') as fh:
            n = 0
            for i in xrange(args.nClusters):
                subM = model.get_indices(i)[idx]
                if list(subM):
                    comment = '#module {n} size: {l}'.format(n=n, l=len(subM))
                    # cluster = map(str, map(lambda x:x+1, subM))
                    cluster = map(lambda x: args.ind2node[x], subM)
                    fh.write(comment + '\n')
                    fh.write(' '.join(cluster) + '\n')
                    n += 1
    print('Saving done')

def plot_spectral(data, fout, args, title):
    '''Function to plot bipartite cluster'''

    print('Plotting')
    # Set figure size
    # plt.figure(figsize=args.psize)

    # Heatmap
    # sns.heatmap(data, xticklabels=[], yticklabels=[])
    # Spy plot
    plt.spy(data, markersize = 3)

    # Annotations
    plt.title(title)
    plt.savefig(args.pf + fout + '.png')
    plt.clf(); # Clear all plots
    print('Plotting done')

if __name__=='__main__':

    # Get file paths and number of clusters
    paths = build_paths()

    # PARSE ARGUMENTS
    parser = argparse.ArgumentParser(description = 'Parse input network files')
    parser.add_argument(
        '-c',
        dest = 'nClusters',
        help = 'Number of clusters. (Default = median number of cluster)',
        default = None,
        type = int)
    parser.add_argument(
        '-arpack',
        dest = 'arpack',
        help = ' '.join([
            'Change svd_method for clustering to "arpack"',
            '(slower for large matrices). Default = "randomized"']),
        nargs = '?',
        default = False,
        const = True,
        type = bool)
    parser.add_argument(
        '-s',
        dest = 'slices',
        nargs = '+',
        help = ' '.join([
            'Slices to work on.',
            'Takes space-separated inputs, i.e. "1 3 14".',
            '(Default = all)']),
        default = ['All'],
        type = str)
    parser.add_argument(
        '-r',
        dest = 'results',
        help = ' '.join([
            'Number of results per slice to include',
            '(Default = all; minimum value = 1)']),
        default = 'All',
        type = str)
    parser.add_argument(
        '-fout',
        dest = 'fout',
        help = 'Output file if saving fitted data',
        type = str)
    parser.add_argument(
        '-df',
        dest = 'df',
        help = 'Folder to save data. (Default = "data")',
        default = 'data',
        type = str)
    parser.add_argument(
        '-p',
        dest = 'plot',
        help = 'Set to plot heatmaps',
        nargs = '?',
        default = False,
        const = True,
        type = bool)
    parser.add_argument(
        '-pf',
        dest = 'pf',
        help = 'Plot folder to save plots. (Default = "plots/")',
        default = 'plots/',
        type = str)
    parser.add_argument(
        '-po',
        dest = 'pOrigName',
        help = ' '.join([
            'Filename for plot of original data before clustering',
            '(Default = "tp_original")']),
        default = 'tp_original',
        type = str)
    parser.add_argument(
        '-pc',
        dest = 'pCluName',
        help = ' '.join([
            'Filename for plot of clustered data',
            '(Default = "tp_clustered")']),
        default = 'tp_clustered',
        type = str)
    parser.add_argument(
        '-ps',
        dest = 'psize',
        nargs = 2,
        help = ' '.join([
            'Width, height of plot (space separated values).',
            '(Default = "10 20")']),
        default = [10, 20],
        type = int)


    args = parser.parse_args()

    # Update args
    args.set_module = None # Necessary ofr other scripts
    args.paths = paths

    if 'All' in args.slices:
        args.slices = args.paths.keys()

    # BUILD MATRIX
    print 'Build Matrix'
    # args.M, nClusters = build_matrix(args)
    args.M, nClusters = build_matrix2(args)
    print('Median cluster size: {}'.format(nClusters))

    if args.nClusters is None:
        args.nClusters = nClusters
    print 'Clusters used = %i'%args.nClusters


    # CLUSTERING
    print 'Clustering'
    print 'Dimensions =', args.M.shape
    cc_model, cc_fit = Spectral_CoClustering(args)
    bc_model, bc_fit = Spectral_BiClustering(args.M.T, args)

    # Plot co-clusters
    if args.plot:
        print 'Plotting'
        plot_spectral(args.M, args.pOrigName + 'no_clustering', args,
            'M: Before Clustering')
        plot_spectral(cc_fit, args.pCluName + 'CoClustering', args,
            '\n'.join(['M: After CoClustering; rearranged to show CoClusters',
                      '{n} clusters used'.format(n=args.nClusters)]))
        plot_spectral(bc_fit.T, args.pCluName + 'BiClustering', args,
            '\n'.join(['M: After BiClustering; rearranged to show BiClusters',
                      '{n} clusters used'.format(n=args.nClusters)]))







