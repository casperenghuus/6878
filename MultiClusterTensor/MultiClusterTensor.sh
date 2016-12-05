#!/bin/bash

opdir="../data/MultiClusterTensor"
[ -d $opdir ] || mkdir -p $opdir
echo "Assembling tensor..."
python2 AssembleTensor.py ../data/net_exp_final.txt ../data/net_tf_final.txt ../data/net_mirna_final.txt > ${opdir}/multi_network.txt

julia4 MultiClusterTensor.jl ${opdir}/multi_network.txt ${opdir}/clusters.txt
python2 ConvertIDs.py -c ${opdir}/clusters.txt -n ${opdir}/nodes.txt -s ${opdir}/slices.txt -o ${opdir}/clusters_renamed.txt
