#!/bin/bash

sed -e 's/BGA/BGL/g' \
    -e 's/BGB/BGL/g' ../01-TOPOLOGY/03-Triple_Helix_3M_residues.pdb >03-Triple_Helix_3M_edit.pdb
