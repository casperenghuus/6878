#! /bin/bash

opdir=data/red-sd
python2 EnrichmentNullMod.py --comm_size_max 500 --comm_size_step 10 --outfile $opdir/p_cut_offs.csv --iters 500 --msigprefix Gene4x/msigdb_filtered --nodes $opdir/nodes.txt
