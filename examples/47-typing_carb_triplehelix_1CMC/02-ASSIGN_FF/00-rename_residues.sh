#!/bin/bash

sed -e 's/BGA/BGL/g' \
    -e 's/BGB/BGL/g' \
    -e 's/BGC/BGL/g' ../01-TOPOLOGY/03-Triple_Helix_3M_1CMC_residues.pdb >03-Triple_Helix_3M_1CMC_edit.pdb
