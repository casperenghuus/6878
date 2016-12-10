#! /bin/bash

opdir=data/full-sd
outname=p_cut_offs_0.csv
# pvals=$(seq 1 9; seq 10 10 100; seq 100 50 500;)
pvals=$(seq 1 9)
echo $pvals
python2 EnrichmentNullMod.py $pvals --outfile $opdir/p_cut_offs_0.csv --iters 500 --msigprefix Gene4x/msigdb_filtered --nodes $opdir/nodes.txt
