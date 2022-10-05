#!bin/bash
nchains=40
natch=121
# First chain
idxs=0
idxe=120

echo "nmols: $nchains" >PE24iPP24_40Ch_headtail.dat
echo "#ich idx-head idx-tail (indexes start at 0)" >>PE24iPP24_40Ch_headtail.dat
for (( i=0; i<$nchains; i++ )) do
    echo "$i $idxs $idxe" >>PE24iPP24_40Ch_headtail.dat
    idxs=`echo $idxs+$natch|bc -l`
    idxe=`echo $idxe+$natch|bc -l`
done
