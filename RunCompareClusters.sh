#! /bin/bash

opdir=data/red-new
# clusterfiles=$(find $opdir \( -wholename '*multi_cluster/clusters_*.txt' -o -iname 'biCl*.csv' -o -wholename '*slice_out_oslom/consensus/tp' \) ! \( -iname '*sym.txt' -o -iname '*enr*.csv' \) -exec bash -c 'printf "%q " "$@"' dummy '{}' +)
clusterfiles=$(find $opdir \( -wholename '*multi_cluster/clusters_*phi_0800.txt' -o -iname 'biCl*900*.csv' -o -iname 'biCl*300*.csv' -o -wholename '*slice_out_oslom/consensus/tp' -o -wholename '*slice_out_oslom/slice_*/results_1/tp' \) ! \( -iname '*sym.txt' -o -iname '*enr*.csv' \) -exec bash -c 'printf "%q " "$@"' dummy '{}' +)
./CompareClustersNMI.py $clusterfiles --outfile $opdir/nmi.csv
