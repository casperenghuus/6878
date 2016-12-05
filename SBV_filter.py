#! /usr/bin/env python2

import graph_tool.all as grat
import ClusteringUtils as cu
import argparse
from scipy import integrate

def compute_sig(g):
    sig = g.new_edge_property('double')
    g.edge_properties['sig'] = sig
    weight = g.edge_properties['weight']

    for e in g.edges():
        sig[e] = float('inf')

    for v in g.vertices():
       nbhd = v.out_neighbours() 
       k_n = v.out_degree()
       if k_n > 1:
           sum_w = sum(float(weight[e]) for e in v.out_edges())
           for e in v.out_edges():
               edgeWeight = float(weight[e])
               p_ij = edgeWeight/sum_w
               f = lambda x: (1-x)**(k_n - 2)
               alpha_ij = 1 - (k_n - 1) * integrate.quad(f, 0, p_ij)[0]
               sig[e] = min(sig[e], alpha_ij)

def extract_backbone(g, alpha):
    backbone = grat.Graph()
    sig = g.edge_properties['sig']
    new_weight = backbone.new_edge_property('string')
    new_names = backbone.new_vertex_property('string')
    names = g.vertex_properties['vertex_name']
    weight = g.edge_properties['weight']
    backbone.edge_properties['weight'] = new_weight
    backbone.vertex_properties['vertex_name'] = names

    for e in g.edges():
        if sig[e] < alpha:
            e_new = backbone.add_edge(e.source(), e.target())
            new_weight[e_new] = weight[e]
            new_names[e_new.source()] = names[e.source()]
            new_names[e_new.target()] = names[e.target()]

    return backbone

parser = argparse.ArgumentParser(description = 'Parse input network files')
parser.add_argument('edgelist', help = 'filename of the network file',
        type = str)
parser.add_argument('dotfile', help = 'filename of temporary storage',
        type = str)
parser.add_argument('prefix', help = 'filename prefix of output',
        type = str)

ns = parser.parse_args()
# nxfname = "./Filtering/net_exp_unique_nodes_int.txt"
nxfname = ns.edgelist
dotfname = ns.dotfile
cu.edge2dot(nxfname, dotfname)
g = grat.load_graph(dotfname)
compute_sig(g)

vet=[0.005,0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.4,0.5]

for i,alpha in enumerate(vet):
    cu.writeEdgeFile(extract_backbone(g, alpha), ns.prefix + str(i) + ".txt")
