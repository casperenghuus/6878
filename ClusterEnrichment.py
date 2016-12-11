#! /usr/bin/env python2

import ClusteringUtils as cu
import argparse
import scipy.stats as stats
import pandas as pd
import numpy as np
import os
import igraph as ig
import scipy.interpolate as interp
from joblib import Parallel, delayed

parser = argparse.ArgumentParser(description = 'Parse input files')
parser.add_argument('clusterfiles', help = 'clusters to analyze',
        type = str, nargs = '+')
parser.add_argument('--outfile', help = 'where to store the p values',
        type = str)
parser.add_argument('--msigprefix', help = 'prefix for msig database files',
        type = str)
parser.add_argument('--nodes', help = 'file to identify the nodes',
        type = str)
parser.add_argument('--graphs', help = 'files to the graphs to compute statistics like conductance',
        type = str, nargs = '+')
parser.add_argument('--thresholds', help = 'file p value thresholds',
        type = str)
ns = parser.parse_args()

with open(ns.nodes, 'r') as f:
    cu.loadNodes(f)

# Read in graph file
gs = []
for fname in ns.graphs:
    with open(fname, 'r'):
        gs.append(ig.Graph().Read_Ncol(fname, directed = False))

# Read db files to analyze
db = cu.readMSig('Gene4x/msigdb/', ['chr', 'cp', 'biocarta', 'reactome', 'go', 'kegg', 'mir'])
for s in db.keys():
    db[s] = {d: set(c) for (d, c) in db[s].items()}
total_length = sum(len(s) for s in db)

gene_count = 45956
# or 54127 :)

summary = pd.DataFrame()
threshold = 0.005
nodes = [str(ind) for ind in cu.ind2node.keys()]

# Threshold interpolation
thresholds = ns.thresholds
thresh_df = pd.read_csv(thresholds)
t_mat = thresh_df.as_matrix()
t_lengths = t_mat[:,0]
t_cutoffs_log = np.log(t_mat[:,1])
# t_interp = interp.interp1d(t_lengths, t_cutoffs_log)
t_interp = interp.InterpolatedUnivariateSpline(t_lengths, t_cutoffs_log, k = 1)

# Read clusters
# for (k, clusterfile) in enumerate(ns.clusterfiles):
def analyze(k, clusterfile, cu, db, thresholds, t_interp, gs, total_length, gene_count, nodes):
    print('File: {}'.format(clusterfile))
    (memb, clusters) = cu.readComms(clusterfile)
    (memb, clusters) = cu.addSingletonClusters(memb, clusters, nodes)
    # Throw away slice nodes
    clusters = [c for c in clusters if len(c) > 0]
    
    clusters_nodes = cu.comms2nodes(clusters)
    
    ps = {dbname: np.zeros((len(clusters_nodes), len(db[dbname]))) for dbname in db.keys()}
    conductances = np.zeros((len(clusters_nodes), len(gs)))
    # ps = np.zeros((len(clusters), len(db)))
    
    # Compute statistics
    summary_row = {}

    # Cluster size
    summary_row['clusters'] = len(clusters)
    cluster_lengths = np.array([len(c) for c in clusters])
    summary_row['cluster length min'] = np.min(cluster_lengths)
    summary_row['cluster length max'] = np.max(cluster_lengths)
    summary_row['cluster length mean'] = np.mean(cluster_lengths)
    summary_row['cluster length median'] = np.median(cluster_lengths)

    # Compute thresholds
    thresholds = np.exp(t_interp(cluster_lengths))

    # Modularity
    for (j, g) in enumerate(gs):
        for v in g.vs:
            v['cluster'] = memb[v['name']]
        # vc = ig.VertexClustering.FromAttribute(g, 'cluster',)
        summary_row['modularity ' + ns.graphs[j]] = g.modularity(g.vs['cluster'], weights = g.es['weight'] if g.is_weighted() else None)

    most_sigs = []
    print('Clusters: {}'.format(len(clusters)))
    for i in range(len(clusters)):
        print(i)
        c_set = clusters_nodes[i]
        c = clusters[i]
        most_sig = {}
        most_sig['threshold'] = thresholds[i]

        # Hypergeometric test
        for (dbname, dbcontent) in db.items():
            for (j, (dname, dgenes)) in enumerate(dbcontent.items()):
                enriched = len(dgenes.intersection(c_set))
                verified = len(dgenes)
                size = len(c_set)
                ps[dbname][i, j] = min(stats.hypergeom.sf(enriched - 1, gene_count, verified, size) * total_length, 1)

            # Filter for most significant hit
            min_ind = np.argmin(ps[dbname][i, :])
            if ps[dbname][i, min_ind] < 1:
                most_sig[dbname + ' entry'] = dbcontent.keys()[min_ind]
                most_sig[dbname + ' p'] = ps[dbname][i, min_ind]
                most_sig[dbname + ' sig'] = ps[dbname][i, min_ind] < most_sig['threshold']

        most_sigs.append(most_sig)

        # Conductances
        for (j, g) in enumerate(gs):
            conductances[i, j] = cu.conductance(g, c)

    summary_row.update({dbname: np.mean(ps[dbname].min(1) < thresholds) for dbname in db.keys()})
    for (j, g) in enumerate(gs):
        summary_row['conductance min ' + ns.graphs[j]] = conductances[:,j].min()
        summary_row['conductance max ' + ns.graphs[j]] = conductances[:,j].max()
        summary_row['conductance mean ' + ns.graphs[j]] = conductances[:,j].mean()
        summary_row['conductance median ' + ns.graphs[j]] = np.median(conductances[:,j])
    
    all_ps = np.concatenate(ps.values(), axis = 1)
    all_vals = np.concatenate([all_ps, conductances], axis = 1)
    all_labels = [key for dbcontent in db.values() for key in dbcontent.keys()]
    all_labels.extend(['conductance ' + g for g in ns.graphs])
    df = pd.DataFrame(all_vals, columns = all_labels)
    basename = os.path.splitext(clusterfile)[0]
    df.to_csv(basename + '_enr.csv')
    most_sigs = pd.DataFrame(most_sigs)
    most_sigs.to_csv(basename + '_enr_high.csv')
    cu.writeComms(clusters_nodes, basename + '_sym.txt')

    return summary_row

summary_rows = Parallel(n_jobs = 1)(delayed(analyze)(k, clusterfile, cu, db, thresholds, t_interp, gs, total_length, gene_count, nodes) for (k, clusterfile) in enumerate(ns.clusterfiles))
summary = pd.DataFrame(summary_rows, index = ns.clusterfiles)
summary.to_csv(ns.outfile)
