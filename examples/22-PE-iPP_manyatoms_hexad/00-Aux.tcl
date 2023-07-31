proc newRep { sel type color rep imol} {
  
    mol selection $sel
    mol representation $type
    mol addrep $imol
    mol showrep $imol $rep on
    mol modcolor $rep $imol $color

}

proc center_of_mass {selection} {
    
    # Returns the center of mass of a selection    

    # some error checking
    if {[$selection num] <= 0} {
        error "center_of_mass: needs a selection with atoms"
    }
    # set the center of mass to 0
    set com [veczero]
    # set the total mass to 0
    set mass 0
    # [$selection get {x y z}] returns the coordinates {x y z}
    # [$selection get {mass}] returns the masses
    # so the following says "for each pair of {coordinates} and masses,
	#  do the computation ..."
    foreach coord [$selection get {x y z}] m [$selection get mass] {
       # sum of the masses
       set mass [expr $mass + $m]
       # sum up the product of mass and coordinate
       set com [vecadd $com [vecscale $m $coord]]
    }
    # and scale by the inverse of the number of atoms
    if {$mass == 0} {
        error "center_of_mass: total mass is zero"
    }
    # The "1.0" can't be "1", since otherwise integer division is done
    return [vecscale [expr 1.0/$mass] $com]
}

proc translate_com_to_origin {imol} {
   
    # Translate the center of mass of a molecule
    # to 0,0,0
    #   imol: id of the molecule, i.e: [molinfo top] or 0 or top

    set all [atomselect $imol "all"]
    set com [center_of_mass $all]
    $all moveby [vecscale -1 $com]
    display resetview

}

proc translate_atom_to_origin {imol idx1} {
   
    # Translate a molecule
    # to 0,0,0 using the atom idx1
    #   imol: id of the molecule, i.e: [molinfo top] or 0 or top

    set all [atomselect $imol "all"]
    set xyz [atomselect $imol "index $idx1" ]
    set com [center_of_mass $xyz]
    puts "cJ: $com"
    $all moveby [vecscale -1 $com]
    display resetview

}

proc rotate_bond_parallel_to_axis {idx1 idx2 axis imol} {

    # Rotate a molecule to be alligened to x, y or z axis
    #   idx1: Index atom 1 in the bond
    #   idx2: Index atom 2 in the bond
    #   axis: "x", "y" or "z"
    #   imol: id of the molecule, i.e: [molinfo top] or 0 or top
    
#    translate_com_to_origin $imol
    set atm1 [atomselect $imol "index $idx1"]
    set atm2 [atomselect $imol "index $idx2"]
    set vecx [expr [$atm1 get x] - [$atm2 get x]]
    set vecy [expr [$atm1 get y] - [$atm2 get y]]
    set vecz [expr [$atm1 get z] - [$atm2 get z]]
    set sel [atomselect top all]
    $sel move [transvecinv "$vecx $vecy $vecz"]

    if { $axis == "y" } {
        puts "Bond $idx1-$idx2 is alligned to Y axis"
        set sel [atomselect $imol "all"]
        $sel move [transaxis z 90]

    } elseif { $axis == "z"} {
        puts "Bond $idx1-$idx2 is alligned to Z axis"
        set sel [atomselect $imol "all"]
        $sel move [transaxis y -90]
    } else {
        puts "Bond $idx1-$idx2 is alligned to X axis"
    }

    display resetview
}

proc create_td {latoms imol {colorid 1} {material Transparent}} {

    # Create the polyhedra for a tetraedal arrangement

    # latoms: List of atoms in the vertex Td
    # imol: id of the molecule, i.e: [molinfo top] or 0 or top

    # Check the input list
    set n [llength $latoms]
    set r 3
    if {$n != 4} {
        puts "List of atoms for Td is not 4"
        return 0
    }
    set idx1 [lindex $latoms 0]
    set idx2 [lindex $latoms 1]
    set idx3 [lindex $latoms 2]
    set idx4 [lindex $latoms 3]


    graphics $imol color $colorid
    graphics $imol material $material

    # Create Td planes
    set atm1 [atomselect $imol "index $idx1"]
    set atm2 [atomselect $imol "index $idx2"]
    set atm3 [atomselect $imol "index $idx3"]
    foreach x1 [$atm1 get {x y z}] x2 [$atm2 get {x y z}] x3 [$atm3 get {x y z}] {
       graphics $imol triangle $x1 $x2 $x3 
    }
    set atm1 [atomselect $imol "index $idx1"]
    set atm2 [atomselect $imol "index $idx2"]
    set atm3 [atomselect $imol "index $idx4"]
    foreach x1 [$atm1 get {x y z}] x2 [$atm2 get {x y z}] x3 [$atm3 get {x y z}] {
       graphics $imol triangle $x1 $x2 $x3 
    }
    set atm1 [atomselect $imol "index $idx1"]
    set atm2 [atomselect $imol "index $idx3"]
    set atm3 [atomselect $imol "index $idx4"]
    foreach x1 [$atm1 get {x y z}] x2 [$atm2 get {x y z}] x3 [$atm3 get {x y z}] {
       graphics $imol triangle $x1 $x2 $x3 
    }
    set atm1 [atomselect $imol "index $idx2"]
    set atm2 [atomselect $imol "index $idx3"]
    set atm3 [atomselect $imol "index $idx4"]
    foreach x1 [$atm1 get {x y z}] x2 [$atm2 get {x y z}] x3 [$atm3 get {x y z}] {
       graphics $imol triangle $x1 $x2 $x3 
    }

}    

proc center_triangle {idx1 idx2 idx3 imol} {

    # Return the center of a triangle
    set at1 [atomselect $imol "index $idx1"] 
    set at2 [atomselect $imol "index $idx2"] 
    set at3 [atomselect $imol "index $idx3"] 
    set cx [expr ([$at1 get x] + [$at2 get x] + [$at3 get x])/3.0]
    set cy [expr ([$at1 get y] + [$at2 get y] + [$at3 get y])/3.0]
    set cz [expr ([$at1 get z] + [$at2 get z] + [$at3 get z])/3.0]
    #graphics $imol sphere [list $cx $cy $cz] radius 0.1

    return [list $cx $cy $cz] 

}

proc mid_point {idx1 idx2 imol} {

    # Return the midpoint between two atoms 
    set at1 [atomselect $imol "index $idx1"] 
    set at2 [atomselect $imol "index $idx2"] 
    set cx [expr ([$at1 get x] + [$at2 get x] )/2.0]
    set cy [expr ([$at1 get y] + [$at2 get y] )/2.0]
    set cz [expr ([$at1 get z] + [$at2 get z] )/2.0]
    #graphics $imol sphere [list $cx $cy $cz] radius 0.1

    return [list $cx $cy $cz] 
}

proc create_symmetry_axis {a b imol {s1 -0.4} {s2 0.4} {colorid 1}} {

    # Create an axis of symmetry with two points for visualization
    # p1, p2: Points of the symmetry axis
    # imol: id of the molecule, i.e: [molinfo top] or 0 or top
    # s1 and s2 scale values for the axis
    # color id of the axis
    
    set ab [vecsub $b $a]

    set p1 [vecadd $a [vecscale $s1 $ab]]
    set p2 [vecadd $b [vecscale $s2 $ab]]

    graphics $imol color $colorid
    graphics $imol cylinder $p1 $p2 radius .02 filled yes
    graphics $imol color 0

    return [vecsub $p1 $p2]

} 

proc rotate_axis {vec deg {molid top}} {

    #DESCRIPTION:
    #   This procedure will rotate the display rather than
    #   the physical coordinates of atoms.
    #PROCEDURES:
    #   rotate_axis - takes as arguments a vector to rotate about and
    #   the number of degees to rotate, as well as the molecule to use
    #   as a reference (defaults to top)

    #EXAMPLE USAGE:
    #   # around the X axis
    #   rotate_axis {1 0 0} 10
    #   # around the 45 degree angle in the x-y plane
    #   rotate_axis {1 1 0} 10
    #   or specify a molecule number (the default is "top")
    #   rotate_axis {1 2 3} 10 5

    # get the current matrix
    lassign [molinfo $molid get rotate_matrix] curr
    # the transformation matrix
    set r [trans axis $vec $deg]
    # get the new matrix
    set m [transmult $r $curr]
    # and apply it to the molecule
    molinfo $molid set rotate_matrix "{ $m }"
}

proc label_atom {selection_string label_string imol {colorid 16} {shift {0.0 0.3 0.0}}} {
    set sel [atomselect top $selection_string]
    if {[$sel num] != 1} {
        error "label_atom: '$selection_string' must select 1 atom"
    }
    # get the coordinates of the atom
    lassign [$sel get {x y z}] coord
    puts "H, $coord, $colorid, $label_string"
    # and draw the text
    graphics $imol color $colorid
    graphics $imol text [vecadd $shift $coord] $label_string size 1 thickness 4
}

proc create_plane_triangle {p1 p2 p3 imol} {

    # Vector paralell to the straight line p1-p2, position vector is p1
    set v1 [vecsub $p1 $p2]
    # Midpoint between p1 and p2
    set m [vecscale 0.5 [vecadd $p1 $p2]]
    # graphics $imol sphere $m radius 0.1
    
    # Straigh line r = $m + t * $v1
    set p1_0 [ vecadd $m [vecscale +1.5 $v1]]
    set p2_0 [ vecadd $m [vecscale -1.5 $v1]]

    # Vector paralell to the straight line m-p3, position vector is m
    set v2 [vecsub $m $p3]
    set p3_0 [ vecadd $m [vecscale -1.5 $v2]]

    graphics $imol triangle $p1_0 $p2_0 $p3_0

}
