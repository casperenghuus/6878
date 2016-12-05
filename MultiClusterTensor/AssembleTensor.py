#!/usr/bin/env python2

from __future__ import print_function

import argparse
import os

outfolder = '../data/MultiClusterTensor'
# print(outfolder)
if not os.path.exists(outfolder):
    os.makedirs(outfolder)

parser = argparse.ArgumentParser(description = 'Parse input network files')
parser.add_argument('files', nargs = '+', help = 'filenames of files corresponding to each layer',
        type = argparse.FileType('r'))

ns = parser.parse_args()
# ns = parser.parse_args(["data/net_exp_final.txt"])

# Just assemble 3-tensor
# Not enough, need to symmetrize

# for (i, f) in enumerate(ns.files):
#     for line in f:
#         print('{} {}'.format(i + 1, line), end = '')

# Read different slices into dictionaries
scount = len(ns.files)
slices = [{} for i in range(scount)]
for (i, f) in enumerate(ns.files):
    for line in f:
        if line != "":
            vals = line.split(" ")
            if len(vals) == 3:
                slices[i][tuple(map(int, vals[0:2]))] = float(vals[2])
            elif len(vals) == 2:
                slices[i][tuple(map(int, vals[0:2]))] = float(1)
    f.close()

# Find node set, make new indices
nodes = set()
for i in range(scount):
    for k in slices[i]:
        nodes.add(k[0])
        nodes.add(k[1])
# print(nodes)

ind2node = { i + 1 + scount: node for (i, node) in enumerate(nodes) }
node2ind = {}
for (ind, node) in ind2node.items():
    node2ind[node] = ind

with open(os.path.join(outfolder, 'nodes.txt'), 'w') as f:
    for (ind, node) in ind2node.items():
        f.write("{} {}\n".format(ind, node))

with open(os.path.join(outfolder, 'slices.txt'), 'w') as f:
    for i in range(scount):
        f.write("{} slice_{}\n".format(i+1, i+1))

# Symmetrize slices
slices_new = [{} for i in range(scount)]
for i in range(scount):
    for ((j, k), val) in slices[i].items():
        perm_val = slices[i].get((k,j))
        if perm_val:
            new_val = (val + perm_val)/2
        else:
            new_val = val
        new_k = node2ind[k]
        new_j = node2ind[j]
        slices_new[i][(new_k,new_j)] = new_val
        slices_new[i][(new_j,new_k)] = new_val

slices = slices_new

# Symmetrize tensor, use 1:(scount+1) as indices for the layers
for i in range(scount):
    for ((j, k), val) in slices[i].items():
        print("  ".join(map(str, [i + 1, j, k, val])))
        print("  ".join(map(str, [j, i + 1, k, val])))
        print("  ".join(map(str, [j, k, i + 1, val])))
