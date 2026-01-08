#!bin/bash
nchains=30
nchA=10
nchB=10
nchC=10

natchA=866
# First kind chain
idxs=1
idxe=861
HEADFILENAME=headtail.dat

echo "nmols: $nchains" >${HEADFILENAME}
echo "#ich idx-head idx-tail (indexes start at 0)" >>${HEADFILENAME}
for (( i=0; i<$nchA; i++ )) do
    echo "$i $idxs $idxe" >>${HEADFILENAME}
    idxs=`echo $idxs+$natchA|bc -l`
    idxe=`echo $idxe+$natchA|bc -l`
    echo "A, $i"
done

# Second kind chain
natchB=866
idxsB=8660
idxeB=9522
nchAB=`echo $nchA +$nchB|bc -l`
for (( i=$nchA; i<$nchAB; i++ )) do
    echo "$i $idxsB $idxeB" >>${HEADFILENAME}
    idxsB=`echo $idxsB+$natchB|bc -l`
    idxeB=`echo $idxeB+$natchB|bc -l`
    echo "B, $i"
done

# Third kind chain
natchC=866
idxsC=17326
idxeC=18182
nchAB=`echo $nchA +$nchB|bc -l`
for (( i=$nchAB; i<$nchains; i++ )) do
    echo "$i $idxsC $idxeC" >>${HEADFILENAME}
    idxsC=`echo $idxsC+$natchC|bc -l`
    idxeC=`echo $idxeC+$natchC|bc -l`
    echo "C, $i"
done
