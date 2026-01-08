#!bin/bash
nchains=10
natch=950
# First chain
idxs=175
idxe=758
PAT=PE_04br_xsd_headtail.dat

echo "nmols: $nchains" >$PAT
echo "#ich idx-head idx-tail (indexes start at 0)" >>$PAT
for (( i=0; i<$nchains; i++ )) do
    echo "$i $idxs $idxe" >>$PAT
    idxs=`echo $idxs+$natch|bc -l`
    idxe=`echo $idxe+$natch|bc -l`
done
