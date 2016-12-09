#! /bin/bash

opdir=red-sd
./CompareClusters.py data/$opdir/slice_out_oslom/slice_0/results_{1,2,3}/tp --outfile data/$opdir/nmis.csv
