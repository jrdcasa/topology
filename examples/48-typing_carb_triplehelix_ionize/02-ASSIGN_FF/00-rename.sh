#!/bin/bash

sed -e 's/BGA/BGL/g' \
    -e 's/BGB/BGL/g' \
    -e 's/BGC/BGL/g' ../01-TOPOLOGY/03-TripleHelix_12M_O6ion_connect_topo_residues.pdb >03-TripleHelix_12M_O6ion_connect_topo_residues_edit.pdb

