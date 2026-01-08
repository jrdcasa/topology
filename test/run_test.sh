#!/bin/bash
WK=`pwd`
source /home/jramos/Programacion/sandboxes/sandbox_common/bin/activate

cd Test-indigo/
python -m unittest test_indigo.py
cd $WK

cd Test-internalcoordinates
python -m unittest test_internalcoordinates.py
cd $WK

cd Test-MolecularGraphs
python -m unittest test_molecular_graph.py
cd $WK 

cd Test-ReadMol2
python -m unittest Test_readmol2.py
cd $WK 

cd Test-ReadPDB
python -m unittest Test_readpdb.py
cd $WK 

cd Test-ReadXSD
python -m unittest Test_readxsd.py
cd $WK

cd Test-Segment
python -m unittest test_Segment.py
cd $WK

cd Test-Topology
python -m unittest test_topology.py
cd $WK
