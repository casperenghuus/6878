#! /bin/bash

opdir=data/full-sd
outname=p_cut_offs.csv
# pvals=$(seq 1 9; seq 10 10 100; seq 100 50 500;)
# pvals=$(seq 1 9)
pvals="1 2"
echo $pvals
python2 EnrichmentNullMod.py $pvals --outfile $opdir/$outname --iters 500 --msigprefix Gene4x/msigdb_filtered --nodes $opdir/nodes.txt
