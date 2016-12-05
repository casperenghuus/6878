#!/usr/bin/env python2

from __future__ import print_function

import argparse
import os

outfolder = 'data/MultiClusterTensor'
if not os.path.exists(outfolder):
    os.makedirs(outfolder)

parser = argparse.ArgumentParser(description = 'Parse input network files')
parser.add_argument('-c', '--clusters', help = 'file with cluster results',
        type = argparse.FileType('r'))
parser.add_argument('-n', '--nodes', help = 'file with node indices',
        type = argparse.FileType('r'))
parser.add_argument('-s', '--slices', help = 'file with slice indices',
        type = argparse.FileType('r'))
parser.add_argument('-o', '--out', help = 'file to store the results',
        type = argparse.FileType('w'))

ns = parser.parse_args()
# ns = parser.parse_args(['-c', os.path.join(outfolder, 'nodes.txt')])
# ns = parser.parse_args(["data/net_exp_final.txt"])

conv = {}
for line in ns.nodes:
    if line != "":
        vals = line.split()
        conv[vals[0]] = vals[1]
ns.nodes.close()

for line in ns.slices:
    if line != "":
        vals = line.split()
        conv[vals[0]] = vals[1]
ns.slices.close()

for line in ns.clusters:
    line = line.strip('\n')
    if line.startswith('#'):
        ns.out.write(line+'\n')
    else:
        ns.out.write(" ".join([conv[s] for s in line.split()]) + "\n")
ns.clusters.close()
ns.out.close()
