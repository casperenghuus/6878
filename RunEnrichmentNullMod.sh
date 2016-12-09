#! /bin/bash

opdir=data/red-sd
python2 EnrichmentNullMod.py --comm_size_max 20 --comm_size_step 10 --outfile $opdir/p_cut_offs.csv --iters 5 --msigprefix Gene4x/msigdb_filtered --nodes $opdir/nodes.txt
