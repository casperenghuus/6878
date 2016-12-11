#! /usr/bin/env python2

import igraph as ig
import networkx as nx
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import numpy as np
import ClusteringUtils as cu

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

cluster_file = 'data/full-sd/slice_out_oslom/consensus/tp'
cluster_file = 'data/full-sd/biCl_full_sd_c13_BiClustering_tp.csv'
graph_file = 'data/full-sd/net_exp_final.txt'
graph_file = 'data/full-sd/net_mirna_final.txt'
graph_file = 'data/full-sd/net_tf_final.txt'
graph_file = 'data/full-sd/net_ppi_unique_nodes_int.txt'

(memb, clusters) = cu.readComms(cluster_file)
g = ig.Graph().Read_Ncol(graph_file, directed = False)
visual_style = {}
# visual_style['layout'] = g.layout('fr')
visual_style['layout'] = layoutByAttr(g, memb, cluster_strength = 100, layout = lambda x: x.layout('fr', weights = 'weight'))
visual_style['vertex_size'] = 5
for v in g.vs:
    v['cluster'] = memb[v['name']]
vc = ig.VertexClustering.FromAttribute(g, 'cluster')
ig.plot(vc, **visual_style)
# ig.plot(g, **visual_style)
