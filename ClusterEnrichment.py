#! /usr/bin/env python2

import ClusteringUtils as cu
import argparse
import scipy.stats as stats
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser(description = 'Parse input files')
parser.add_argument('clusterfiles', help = 'clusters to analyze',
        type = str, nargs = '+')
parser.add_argument('--outfile', help = 'where to store the p values',
        type = str)
parser.add_argument('--msigprefix', help = 'prefix for msig database files',
        type = str)
parser.add_argument('--nodes', help = 'file to identify the nodes',
        type = str)
ns = parser.parse_args()

with open(ns.nodes, 'r') as f:
    cu.loadNodes(f)

# Read db files to analyze
db = cu.readMSig('Gene4x/msigdb/', ['chr', 'biocarta'])
for s in db.keys():
    db[s] = {d: set(c) for (d, c) in db[s].items()}
total_length = sum(len(s) for s in db)

gene_count = 45956
# or 54127 :)

fracs = pd.DataFrame()
threshold = 0.05

# Read clusters
for (k, clusterfile) in enumerate(ns.clusterfiles):
    (memb, clusters) = cu.readComms(clusterfile)
    
    clusters = cu.comms2nodes(clusters)
    
    ps = {dbname: np.zeros((len(clusters), len(db[dbname]))) for dbname in db.keys()}
    # ps = np.zeros((len(clusters), len(db)))
    
    # Hypergeometric test
    for (i, c) in enumerate(clusters):
        c_set = c
        for (dbname, dbcontent) in db.items():
            for (j, (dname, dgenes)) in enumerate(dbcontent.items()):
                enriched = len(dgenes.intersection(c_set))
                verified = len(dgenes)
                size = len(c_set)
                ps[dbname][i, j] = min(stats.hypergeom.sf(enriched - 1, gene_count, verified, size) * total_length, 1)

    fracs = fracs.append(pd.DataFrame([{dbname: np.mean(ps[dbname].min(1) < threshold) for dbname in db.keys()}], index = ['clusterfile']))
    
    all_ps = np.concatenate(ps.values(), axis = 1)
    all_labels = [key for dbcontent in db.values() for key in dbcontent.keys()]
    df = pd.DataFrame(all_ps, columns = all_labels)
    df.to_csv(os.path.splitext(clusterfile)[0] + '_enr.csv')

fracs.to_csv(ns.outfile)
