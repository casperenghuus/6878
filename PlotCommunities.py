#! /usr/bin/env python2

import igraph as ig
import networkx as nx
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import numpy as np
import ClusteringUtils as cu

def spAdjMat(graph):
    xs, ys = map(array, zip(*graph.get_edgelist()))
    if not graph.is_directed():
        xs, ys = np.hstack((xs, ys)).T, np.hstack((ys, xs)).T
    else:
        xs, ys = xs.T, ys.T
    return coo_matrix((np.ones(xs.shape), (xs, ys)))

def readComms(fname):
    membership = {}
    clusters = []
    counter = 0
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip('\n')
            if not line.startswith('#'):
                cluster = []
                clusters.append(cluster)
                for node in line.split(' '):
                    if not node.startswith('slice'):
                        membership[node] = counter
                        cluster.append(node)
                counter += 1
    return (membership, clusters)

def clusterToInd(g, c):
    ret = []
    for name in c:
        try:
            v = g.vs.find(name)
            ret.append(v.index)
        except ValueError:
            pass

    return ret

def layoutByAttr(g, memb, cluster_strength=100, layout = lambda x: x.layout('fr')):
    g_new = g.copy() # create a lightweight copy of graph w/o the attributes.
    g_new.es['weight'] = 1
    cluster_edges = [(key, 'cluster{}'.format(val)) for (key, val) in memb.items()]
    cluster_vertices = ['cluster{}'.format(val) for val in set(memb.values())]
    g_new.add_vertices(cluster_vertices)
    for e in cluster_edges:
        try:
            g_new.add_edge(*e, weight = cluster_strength)
        except ValueError: 
            pass
    l = layout(g_new)[0:len(g.vs)]
    return l

# def conductance(g, c):
#     """ Compute conductance of cluster c in graph g """

#     # New
#     inds = clusterToInd(g, c)
#     A = spAdjMat(g).tocsc()
#     c_mask = np.zeros(A.shape[0], dtype = bool)
#     c_mask[inds] = True
#     comp_mask = np.logical_not(c_mask)

#     intra_w = A[c_mask, :][:, comp_mask].sum()
#     inter_c_w = A[c_mask, :][:, c_mask].sum()
#     inter_comp_w = A[comp_mask, :][:, comp_mask].sum()
#     # print(intra)
#     # print(inter_c)
#     # print(inter_comp)
#     if intra_w == 0:
#         cond = 0
#     elif inter_c_w == 0 or inter_comp_w == 0:
#         cond = float('inf')
#     else:
#         cond = intra_w/min(inter_c_w, inter_comp_w) 
#     return cond

(memb, clusters) = readComms('data/MultiClusterTensor/clusters_renamed.txt')
g = ig.Graph().Read_Ncol('data/net_mirna_final.txt', directed = False)
visual_style = {}
# visual_style['layout'] = g.layout('fr')
visual_style['vertex_size'] = 5
for v in g.vs:
    v['cluster'] = memb[v['name']]
vc = ig.VertexClustering.FromAttribute(g, 'cluster')
# ig.plot(vc)
# ig.plot(g, **visual_style)

# g = nx.read_weighted_edgelist('data/net_exp_final.txt')

# Have a look at http://stackoverflow.com/questions/28715736/how-to-spread-out-community-graph-made-by-using-igraph-package-in-r
