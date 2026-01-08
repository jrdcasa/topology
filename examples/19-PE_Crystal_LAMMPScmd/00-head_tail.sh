#!bin/bash
nchains=8
nchA=8
nchB=0
natchA=1000
# First kind chain
idxs=0
idxe=999
HEADFILENAME=PE_C1000_008Ch_headtail.dat

echo "nmols: $nchains" >${HEADFILENAME}
echo "#ich idx-head idx-tail (indexes start at 0)" >>${HEADFILENAME}
for (( i=0; i<$nchA; i++ )) do
    echo "$i $idxs $idxe" >>${HEADFILENAME}
    idxs=`echo $idxs+$natchA|bc -l`
    idxe=`echo $idxe+$natchA|bc -l`
done

# Second kind chain
natchB=0
idxsB=60561
idxeB=62829
for (( i=$nchA; i<$nchains; i++ )) do
    echo "$i $idxsB $idxeB" >>${HEADFILENAME}
    idxsB=`echo $idxsB+$natchB|bc -l`
    idxeB=`echo $idxeB+$natchB|bc -l`
done

