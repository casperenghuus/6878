#!/usr/bin/env python

import os
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
from sklearn.cluster.bicluster import SpectralCoclustering

pipelinePath = 'FortunatoPipelinePath.txt'

def build_paths():
    '''Build list of paths for tp files'''

    # Path to parent folder holding all tp files
    with open(pipelinePath, 'r') as fh:
        fortunatoPath = fh.readline().rstrip()

    # Regex to get file and capture cluster number/slice number
    f_re = re.compile('.*/slice_(\d+)/(\d+)/tp$')

    # Number of clusters
    nClusters = 0

    # Get tp paths, their slice number and their cluster number
    paths = {}
    for root, dirs, files in os.walk(fortunatoPath):
        for f in files:
            match = f_re.match(os.path.join(root,f))
            if match:
                nClusters += 1
                p = match.group(0) #Path
                s = match.group(1) #Slice
                c = match.group(2) #Cluster

                if s in paths.keys():
                    paths[s].append({'path':p, 'cluster':c})
                else:
                    paths[s] = [{'path':p, 'cluster':c}]
    return paths, nClusters

def build_matrix(args):
    '''Build matrix from tp column vectors'''

    M = None
    col_idx = []

    for s in args.slices:
        nResults = 0 # Track number of results
        try:
            for p in args.paths[s]:
                # Name for column based on tp file
                colName = 's' + s + 'c' + p['cluster']
                col_idx.append(colName)
                # Create matrix
                if M is None:
                    M = tp_to_col(p['path'])
                # Append columns to matrix
                else:
                    C = tp_to_col(p['path'])
                    M = append_col(C, M)

                # Break loop if a maximum number of results has been set
                if args.results != 'All':
                    nResults += 1
                    if nResults == int(args.results):
                        break

        except KeyError, e:
            print '%s not a valid slice'%e

    return M, col_idx

def append_col(C,M):
    '''Append column to matrix. If their sizes differ, buffer by zeros'''

    # Matrix/column dimensions
    shapeM = M.shape
    sizeCol = C.size

    # Test if there is a difference in number of rows
    diff = shapeM[0] - sizeCol

    # Append column to matrix as new column
    if diff == 0:
        M = np.hstack([M, C])
    # Add extra zero-rows if new column has less rows than matrix.
    #  Append to matrix as new column
    elif diff > 0:
        extraRows = np.ones([diff, 1])
        C = np.vstack([C, extraRows])
        M = np.hstack([M, C])
    # Add extra zero-rows if matrix has less rows than new column.
    #  Append column to matrix as new column
    elif diff < 0:
        extraRows = np.ones([abs(diff), shapeM[1]])
        M = np.vstack([M, extraRows])
        M = np.hstack([M, C])

    return M

def tp_to_col(p):
    '''
    Read tp file column vector. Module information
    is stored in case it becomes useful later.
    '''

    # Dict of modules for tp file
    modules = {}

    # Regex to capture module information
    s = re.compile('module\s+(\d+)\s+size:\s+(\d+)\s+bs:\s+(\d+\.{0,1}\d*)')

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

    # Create column vector
    col = np.ones([max_val, 1], dtype=float)
    for k in modules.keys():
        for i in modules[k]['nodes']:
            col[i - 1,0] = 10

    return col

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
    model.fit(args.M)

    # Fit to data
    fit_data = args.M[np.argsort(model.row_labels_)]
    fit_data = args.M[:, np.argsort(model.column_labels_)]

    # Save fitted data if specified
    if args.fout is not None:
        pdx = pd.DataFrame(
            fit_data,
            index=np.argsort(model.row_labels_),
            columns=np.array(args.col_idx)[np.argsort(model.column_labels_)]
        )
        pdx.to_csv(args.fout + '.csv', index=True)

    return model, fit_data

def plot_spectral(data, fout, col_idx, args, title):
    '''Function to plot bipartite cluster'''

    # Set figure size
    plt.figure(figsize=args.psize)

    # Heatmap
    sns.heatmap(data, xticklabels=col_idx)

    # Annotations
    plt.xticks(rotation='vertical')
    plt.title(title)
    plt.savefig(args.pf + fout + '.png')
    plt.clf(); # Clear all plots

if __name__=='__main__':

    # Get file paths and number of clusters
    paths, nClusters = build_paths()

    # PARSE ARGUMENTS
    parser = argparse.ArgumentParser(description = 'Parse input network files')
    parser.add_argument(
        '-c',
        dest = 'nClusters',
        help = 'Number of clusters. (Default = number of tp files)',
        default = nClusters,
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
    args.paths = paths
    if 'All' in args.slices:
        args.slices = args.paths.keys()

    # BUILD MATRIX
    args.M, args.col_idx = build_matrix(args)

    # CO-CLUSTERING
    cc_model, cc_fit = Spectral_CoClustering(args)

    # Plot co-clusters
    if args.plot:
        plot_spectral(args.M, args.pOrigName, args.col_idx, args,
            'M: Before coclustering')
        plot_spectral(cc_fit, args.pCluName, [], args,
            'M: After coclustering; rearranged to show coclusters')







