import ClusteringUtils as cu
import pandas as pd
import numpy as np
import time
import argparse
import os

parser = argparse.ArgumentParser(description = 'Parse input network files')
parser.add_argument('--prefix', help = 'prefix of network files to analyze',
        type = str)
parser.add_argument('--outfile', help = 'where to store the results',
        type = str)
parser.add_argument('--dbfolder', help = 'folder containing the database files',
        type = str)
parser.add_argument('--nodes', help = 'file to identify the nodes',
        type = str)
ns = parser.parse_args()

t0 = time.time()
# Read in nodes file
with open(ns.nodes) as nodeFile:
    cu.loadNodes(nodeFile)

# Read in database files
dbFolder = ns.dbfolder
dbNames = ['BioGrid', 'cell', 'hprd',
        'intact', 'mint', 'nat', 'react',]
# dbNames = ['imid']
dbFnames = [os.path.join(dbFolder, dbName + '.csv') for dbName in dbNames]
dbdf = pd.DataFrame()
for dbFname in dbFnames:
    df = pd.read_csv(dbFname)
    dbdf = dbdf.append(df[['V1', 'V2']])
confEdges = cu.convertToTuples(dbdf)
confNodes = set(dbdf.iloc[:,0].append(dbdf.iloc[:,1]))

linksizes = []
nodesizes = []
densities = []
pvals = []

# Read expression files
for i in range(10):
    print(i)
    expFname = '{}{}.txt'.format(ns.prefix, i)
    try:
        df = pd.read_csv(expFname, header = None, delimiter = ' ')
        df = df.iloc[:,0:2]
        df = cu.convertToNodes(df)
        expEdges = set(cu.convertToTuples(df))
        expNodes = set(df.iloc[:,0].append(df.iloc[:,1]))
    
        p = cu.enrichmentTest(confNodes, confEdges, expNodes, expEdges)
        
        # Fisher test
        nPredicted =len(expEdges)
        nodes = len(expNodes)
        nUniverse = nodes * (nodes - 1) / 2
        linksizes.append(nPredicted)
        nodesizes.append(len(expNodes))
        densities.append(float(nPredicted)/nUniverse)
        pvals.append(p)
    except pd.io.common.EmptyDataError:
        nodesizes.append(0)
        linksizes.append(0)
        densities.append(0)
        pvals.append(1)

alphas = [0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
# alphas = alphas[0:2]
df = pd.DataFrame({'alpha' : alphas, 'Links' : linksizes, 'p-value': pvals, 'nodes' : nodesizes, 'density' : densities})
df = df.reindex_axis(['alpha', 'Links', 'p-value', 'nodes', 'density'], axis=1)
df.to_csv(ns.outfile, sep = '\t', index = False)

t1 = time.time()
print('Time for filtering {}: {} minutes'.format(ns.outfile, (t1-t0)/60))
