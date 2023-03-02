#!bin/bash
nchains=67
nchA=40
nchB=27
natchA=1514
# First kind chain
idxs=1
idxe=1508
HEADFILENAME=PEiPP_film_252mon_40-27Ch_headtail.dat

echo "nmols: $nchains" >${HEADFILENAME}
echo "#ich idx-head idx-tail (indexes start at 0)" >>${HEADFILENAME}
for (( i=0; i<$nchA; i++ )) do
    echo "$i $idxs $idxe" >>${HEADFILENAME}
    idxs=`echo $idxs+$natchA|bc -l`
    idxe=`echo $idxe+$natchA|bc -l`
done

# Second kind chain
natchB=2273
idxsB=60561
idxeB=62829
for (( i=$nchA; i<$nchains; i++ )) do
    echo "$i $idxsB $idxeB" >>${HEADFILENAME}
    idxsB=`echo $idxsB+$natchB|bc -l`
    idxeB=`echo $idxeB+$natchB|bc -l`
done

