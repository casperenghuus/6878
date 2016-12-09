#! /usr/bin/env python2

import os

msigprefix = '../Gene4x/msigdb'
outprefix = '../Gene4x/msigdb_filtered'
cut_off = 500

for fname in os.listdir(msigprefix):
    with open(os.path.join(msigprefix, fname), 'r') as f, open(os.path.join(outprefix, fname), 'w') as fout:
        for line in f:
            if len(line.split('\t')) <= cut_off + 2:
                fout.write(line)
