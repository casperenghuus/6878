#! /usr/bin/env python2

import argparse
import ClusteringUtils as cu
import sklearn.metrics
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description = 'Parse input files')
parser.add_argument('clusterfiles', help = 'clusters to compare',
        type = str, nargs = '+')
parser.add_argument('--outfile', help = 'clusters to compare',
        type = str)
ns = parser.parse_args()

clusterings = []
membs = []
for cfname in ns.clusterfiles:
    (memb, clusters) = cu.readComms(cfname)
    clusterings.append(clusters)
    membs.append(memb)

# if set(memb1.keys()) != set(memb2.keys()):
#     print('Warning! Node sets are not the same.')

n = len(clusterings)
nmis = np.zeros((n, n))

for i in range(n):
    for j in range(i+1, n):
        print(i, j)
        memb1 = membs[i]
        memb2 = membs[j]
        memb1_list = []
        memb2_list = []
        for k in memb1.keys():
            if k in memb2.keys():
                memb1_list.append(memb1[k])
                memb2_list.append(memb2[k])
        nmis[i, j] = sklearn.metrics.normalized_mutual_info_score(memb1_list, memb2_list)

nmis += nmis.T
np.fill_diagonal(nmis, 1)

df = pd.DataFrame(nmis, columns = ns.clusterfiles, index = ns.clusterfiles)
df.to_csv(ns.outfile)
