#!/bin/bash

# cp ../../02-CHARMM_GUI/charmm-gui-0808503697/gromacs/step3_input.gro .
# Remove water with an editor
# ans save as 36CURDLAN_extended.gro
# cp -R ../../02-CHARMM_GUI/charmm-gui-0808503697/gromacs/top* .
# cp ../../../01-24CURDLAN13/03-GROMACS/01a-RUN01_45d/minimization.mdp .

GRO_INP=../02-ASSIGN_FF/03-Triple_Helix_3M_1CMC_edit.pdb
TOP=../02-ASSIGN_FF/topopoly.top
MDP=minimization.mdp

GRO_BOX=03-Triple_Helix_3M_1CMC_dodec.gro
GRO_SOLV=03-Triple_Helix_3M_1CMC_dodec_solv.gro

# GROMACS 2024
source /opt/gromacs/gromacs-2024_x/bin/GMXRC.bash

# EDITCONF
gmx editconf -f ${GRO_INP} -o ${GRO_BOX} -d 0.5 -bt dodecahedron -rotate 0 45. 0 -c

# SOLVATE
gmx solvate -cp ${GRO_BOX} -o ${GRO_SOLV} -p ${TOP}

# Change SOLV by TIP3
sed -e "/^TIP3/d" ${TOP} >tmp
#sed -e "s/SOL/TIP3/g" tmp >${TOP}
rm tmp

# Change SOLV by TIP3 and atom types of water
#sed -e "s/SOL/TIP3/g" \
#    -e "s/  OW/OH2/g" \
#    -e "s/HW1/H1/g" \
#    -e "s/HW2/H2/g" ${GRO_SOLV} >tmp
#mv tmp ${GRO_SOLV}

# Generate topol.tpr for minimization
gmx grompp -c ${GRO_SOLV} -f ${MDP} -p ${TOP}  -r ${GRO_SOLV} -maxwarn 2

# Visualization
GRO_VMD=`echo ${GRO_SOLV%.*}`_vmd.gro
echo 0 |gmx trjconv -f ${GRO_SOLV} -pbc mol -ur compact -o ${GRO_VMD}

# Delete backup files
rm \#*

