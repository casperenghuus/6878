import ClusteringUtils as cu
import pandas as pd
import numpy as np
import scipy.stats as stats
import time

# Read in nodes file
ut = cu.Utils()
with open('data/nodes.txt') as nodeFile:
    ut.loadNodes(nodeFile)

# Read in database files
dbFolder = 'data/Reference/'
# dbNames = ['BCI', 'BioGrid', 'cell', 'cyc', 'hprd',
#         'imid', 'intact', 'mint', 'nat', 'react', 'UniHI']
dbNames = ['BioGrid', 'cell', 'hprd',
        'intact', 'mint', 'nat', 'react',]
# dbNames = ['imid']
dbFnames = [dbFolder + dbName + '.csv' for dbName in dbNames]
dbdf = pd.DataFrame()
for dbFname in dbFnames:
    df = pd.read_csv(dbFname)
    dbdf = dbdf.append(df[['V1', 'V2']])
confEdges = ut.convertToTuples(dbdf)
confNodes = set(dbdf.iloc[:,0].append(dbdf.iloc[:,1]))

linksizes = []
nodesizes = []
densities = []
pvals = []

# Read expression files
for i in range(10):
    print(i)
    expFname = 'Filtering/test.edgelist_exp{}.txt'.format(i)
    df = pd.read_csv(expFname, header = None, delimiter = ' ')
    df = df.iloc[:,0:2]
    df = ut.convertToNodes(df)
    expEdges = set(ut.convertToTuples(df))
    expNodes = set(df.iloc[:,0].append(df.iloc[:,1]))
    
    # Filter database
    focusedSet = set([(a, b) for (a, b) in confEdges if a in expNodes and b in expNodes])
    
    # Fisher test
    nPredicted =len(expEdges)
    nOverlap = len(expEdges.intersection(focusedSet))
    nFocusedSet = len(focusedSet)
    nodes = len(expNodes)
    nUniverse = nodes * (nodes - 1) / 2
    linksizes.append(nPredicted)
    nodesizes.append(len(expNodes))
    densities.append(float(nPredicted)/nUniverse)
    fe1 = nOverlap
    fe2 = nPredicted - nOverlap
    fe3 = nFocusedSet - nOverlap
    fe4 = nUniverse - nPredicted - nFocusedSet + nOverlap 
    table = np.array([[fe1, fe2], [fe3, fe4]])
    (_, p) = stats.fisher_exact(table, alternative = 'greater')
    pvals.append(p)
    print(table)
    print(p)

alphas = [0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
# alphas = alphas[0:2]
df = pd.DataFrame({'alpha' : alphas, 'Links' : linksizes, 'p-value': pvals, 'nodes' : nodesizes, 'density' : densities})
df = df.reindex_axis(['alpha', 'Links', 'p-value', 'nodes', 'density'], axis=1)
df.to_csv('Enrichment/selection_links_EXP.txt', sep = '\t', index = False)
