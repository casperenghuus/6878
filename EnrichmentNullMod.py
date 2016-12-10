#! /usr/bin/env python2

import ClusteringUtils as cu
import argparse
import scipy.stats as stats
import pandas as pd
import numpy as np
import os
import igraph as ig
from itertools import chain

parser = argparse.ArgumentParser(description = 'Parse input files')
parser.add_argument('community_sizes', help = 'different sizes to compute the p-value for',
        type = int, nargs = '*')
# parser.add_argument('--comm_size_max', help = 'max sizes to compute the p-value for',
#         type = int)
# parser.add_argument('--comm_size_step', help = 'step to compute the p-value for',
#         type = int)
parser.add_argument('--outfile', help = 'where to store the p values',
        type = str)
parser.add_argument('--iters', help = 'how many iterations to run',
        type = int)
parser.add_argument('--msigprefix', help = 'prefix for msig database files',
        type = str)
parser.add_argument('--nodes', help = 'file to identify the nodes',
        type = str)
ns = parser.parse_args()

with open(ns.nodes, 'r') as f:
    nodes = cu.loadNodes(f)

# Read db files to analyze
db = cu.readMSig('Gene4x/msigdb/', ['chr', 'cp', 'biocarta', 'reactome', 'go', 'kegg', 'mir'])
for s in db.keys():
    db[s] = {d: set(c) for (d, c) in db[s].items()}
total_length = sum(len(s) for s in db)

gene_count = 45956
# or 54127 :)

threshold = 0.05
cut_offs = []

# Generate p values
for length in ns.community_sizes:
# size_range = list(chain(range(1, 10), range(ns.comm_size_step, ns.comm_size_max + 1, ns.comm_size_step)))
# for length in size_range:
    print('Cluster length: {}'.format(length))
    # Generate cluster
    ps = {dbname: np.zeros((ns.iters, len(db[dbname]))) for dbname in db.keys()}
    # ps = np.zeros((len(clusters), len(db)))
    
    # Compute statistics
    for i in range(ns.iters):
        print('Iteration: {}'.format(i))
        c_set = np.random.choice(np.array(nodes), size = length, replace = False)

        # Hypergeometric test
        for (dbname, dbcontent) in db.items():
            for (j, (dname, dgenes)) in enumerate(dbcontent.items()):
                enriched = len(dgenes.intersection(c_set))
                verified = len(dgenes)
                size = len(c_set)
                ps[dbname][i, j] = min(stats.hypergeom.sf(enriched - 1, gene_count, verified, size) * total_length, 1)

    all_ps = np.concatenate(ps.values(), axis = 1)
    cut_off = np.percentile(all_ps, threshold)
    cut_offs.append(cut_off)

# df = pd.DataFrame(cut_offs, columns = ['cut-off'], index = size_range)
df = pd.DataFrame(cut_offs, columns = ['cut-off'], index = ns.community_sizes)
df.index.name = 'length'
df.to_csv(ns.outfile)
