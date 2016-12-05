import re
import scipy.stats as stats
import numpy as np
import time

ind2node = {}
node2ind = {}

def loadNodes(nodeFile):
    global ind2node, node2ind
    """ Load node indices in nodes.txt file """
    i = 1
    for line in nodeFile:
        if line != "":
            line_s = line.strip()
            line_s = line_s.strip('"')
            ind2node[i] = line_s
            node2ind[line_s] = i
        i += 1

def convertToNodes(df):
    """ Convert string entries in data frame to indices """
    return df.applymap(lambda x: ind2node[x])

def convertToInds(df):
    """ Convert indices in data frame to strings """
    return df.applymap(lambda x: node2ind[x])

def sortedTuples(s):
    """ Sort the tuples in the given set """
    return map(lambda x: tuple(sorted(x)), s)

def convertToTuples(df):
    """ Convert edge dataframe to edge tuples for intersection testing """
    return df.apply(lambda x: tuple(sorted(x)), axis = 1).values

def nx2dot(infname, outfname):
    prog = re.compile(r'(\d*) (\d*) {\'weight\': (\d*\.\d*)}')

    with open(infname, 'r') as infile, open(outfname, 'w') as outfile:
        outfile.write('strict graph {\n')
        for line in infile:
            if line != "":
                m = prog.match(line)
                if m:
                    outfile.write('  ' + m.group(1) + ' -- ' + m.group(2) + ' [weight = ' + m.group(3) + ']\n')

        outfile.write('}\n')

def edge2dot(infname, outfname):
    with open(infname, 'r') as infile, open(outfname, 'w') as outfile:
        outfile.write('strict graph {\n')
        for line in infile:
            if line != "":
                m = line.strip('\n')
                m = m.split('\t')
                if m:
                    outfile.write('  ' + m[0] + ' -- ' + m[1] + ' [weight = ' + m[2] + ']\n')

        outfile.write('}\n')

def writeEdgeFile(g, outfname):
    """ Write edge file from graph-tools graph g """
    weight = g.edge_properties['weight']
    names = g.vertex_properties['vertex_name']

    with open(outfname, 'w') as outfile:
        for e in g.edges():
            outfile.write('{} {} {{\'weight\': {}}}\n'.format(names[e.source()], names[e.target()], weight[e]))

def writeEdgeFileIG(g, outfname):
    """ Write edge file from igraph graph g """

    with open(outfname, 'w') as outfile:
        for e in g.es:
            outfile.write('{} {} {{\'weight\': {}}}\n'.format(g.vs[e.source]['name'], g.vs[e.target]['name'], e['weight']))

def enrichmentTest(confNodes, confEdges, testNodes, testEdges):
    # Filter database
    focusedSet = set([(a, b) for (a, b) in confEdges if a in testNodes and b in testNodes])
    
    # Fisher test
    nPredicted =len(testEdges)
    nOverlap = len(testEdges.intersection(focusedSet))
    nFocusedSet = len(focusedSet)
    nodes = len(testNodes)
    nUniverse = nodes * (nodes - 1) / 2
    fe1 = nOverlap
    fe2 = nPredicted - nOverlap
    fe3 = nFocusedSet - nOverlap
    fe4 = nUniverse - nPredicted - nFocusedSet + nOverlap 
    table = np.array([[fe1, fe2], [fe3, fe4]])
    (_, p) = stats.fisher_exact(table, alternative = 'greater')
    return p

