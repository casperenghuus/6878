#! /usr/bin/env python2

# import graph_tool.all as grat
import igraph as ig
import time
import ClusteringUtils as cu
import argparse
from scipy import integrate

def compute_sig(g):
    # sig = g.new_edge_property('double')
    # g.edge_properties['sig'] = sig
    # weight = g.edge_properties['weight']

    g.es['sig'] = float('inf')

    for v in g.vs:
       nbhd = v.neighbors() 
       k_n = v.degree()
       if k_n > 1:
           conn_edges = g.es[g.incident(v, mode = ig.ALL)]
           sum_w = sum(float(e['weight']) for e in conn_edges)
           for e in conn_edges:
               edgeWeight = float(e['weight'])
               p_ij = edgeWeight/sum_w
               f = lambda x: (1-x)**(k_n - 2)
               alpha_ij = 1 - (k_n - 1) * integrate.quad(f, 0, p_ij)[0]
               e['sig'] = min(e['sig'], alpha_ij)

def extract_backbone(g, alpha):
    backbone = g.es.select(sig_lt = alpha).subgraph()
    return backbone

parser = argparse.ArgumentParser(description = 'Parse input network files')
parser.add_argument('edgelist', help = 'filename of the network file',
        type = str)
parser.add_argument('dotfile', help = 'filename of temporary storage',
        type = str)
parser.add_argument('prefix', help = 'filename prefix of output',
        type = str)

t0 = time.time()
ns = parser.parse_args()
# # nxfname = "./Filtering/net_exp_unique_nodes_int.txt"
# nxfname = ns.edgelist
# dotfname = ns.dotfile
# cu.edge2dot(nxfname, dotfname)
g = ig.Graph().Read_Ncol(ns.edgelist, directed = False)
g.simplify(combine_edges = 'first')
compute_sig(g)

vet=[0.005,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.4,0.5]

for i,alpha in enumerate(vet):
    cu.writeEdgeFileIG(extract_backbone(g, alpha), ns.prefix + str(i) + ".txt")

t1 = time.time()
print('Time for filtering {}: {} minutes'.format(ns.edgelist, (t1-t0)/60))
