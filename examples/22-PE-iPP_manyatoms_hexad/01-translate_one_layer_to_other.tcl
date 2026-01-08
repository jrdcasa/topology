#=============== MAIN PROGRAM ====================
source 00-Aux.tcl
display projection orthographic
#axes location off
#color change rgb red2 0.31 0.0 0.36
#color Display Background "red2"
display depthcue off

# Load molecule 1 ===================
mol new PE-iPP_more_99999_atoms.pdb
set imol1 [molinfo top]
mol delrep 0 $imol1
set rep1 0
#newRep "index <=48959" "VDW 0.200000 12.000000" "colorID 0" $rep1 $imol1
newRep "index <=48959" "lines" "colorID 0" $rep1 $imol1
set rep1 [expr $rep1+1]
#newRep "index >48950" "VDW 0.200000 12.000000" "colorID 1" $rep1 $imol1
newRep "index >48959" "lines" "colorID 1" $rep1 $imol1

rotate x by -90
pbc box

set layertop [atomselect $imol1 "index >48959"]
set layerbot [atomselect $imol1 "index <=48959"]

$layertop moveby {0 0 +10}
$layerbot moveby {0 0 -5}

animate write pdb PE-iPP_more_99999_atoms_nooverlap.pdb
