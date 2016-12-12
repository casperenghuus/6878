#! /bin/bash
opdir=data/red-new
# clusterfiles=$(find $opdir \( -wholename '*multi_cluster/clusters_*.txt' -o -iname 'biCl*.csv' -o -wholename '*slice_out_oslom/consensus/tp' \) ! \( -iname '*sym.txt' -o -iname '*enr*.csv' \) -exec bash -c 'printf "%q " "$@"' dummy '{}' +)
clusterfiles=$(find $opdir \( -iname 'biCl*.csv' -o -wholename '*slice_out_oslom/consensus/tp' \) ! \( -iname '*sym.txt' -o -iname '*enr*.csv' \) -exec bash -c 'printf "%q " "$@"' dummy '{}' +)
# clusterfiles=$(find $opdir \( -iname 'biCl*.csv' -o -wholename '*slice_out_oslom/consensus/tp' \) ! \( -iname '*sym.txt' -o -iname '*enr*.csv' \) -exec bash -c 'printf "%q " "$@"' dummy '{}' +)
# clusterfiles=$opdir/biCl_red_sd_c25_BiClustering_tp.csv
# echo $clusterfiles
python2 ClusterEnrichment.py $clusterfiles --nodes $opdir/nodes.txt --msigprefix Gene4x/msigdb/ --outfile $opdir/overview_2.csv --graphs $opdir/net_* --thresholds $opdir/p_cut_offs.csv
