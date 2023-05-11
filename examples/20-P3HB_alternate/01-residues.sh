#!/bin/bash

echo "P3H 0 0" >residues_tmp.dat

nres=28
for (( ires=1; ires<$nres; ires+=2 )); do
    iresp=`echo $ires+1|bc -l`
    echo "P3A $ires $ires" >>residues_tmp.dat
    echo "P3B $iresp $iresp" >>residues_tmp.dat
done
last=`echo $nres+1|bc -l`
echo "P3T $last $last" >>residues_tmp.dat
