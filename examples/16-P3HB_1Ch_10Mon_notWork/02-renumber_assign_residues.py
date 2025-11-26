import datetime
from topology.readmol.readxsdformat import ReadXsdFormat
from topology.readmol.readpdbformat import ReadPdbFormat


filenamepdb = "./P3HB_10mon_01c_xsd.pdb"
pattern="./P3HB_10mon_01c_xsd"
headinfo_file = "./P3HB_10mon_01c_headtail.dat"
residueinfo_file="./P3HB_10mon_01c_residues.dat"

# Renumber PDB
pdb = ReadPdbFormat(filenamepdb)
head_atoms, tail_atoms = pdb.read_head_tail_info(headinfo_file)
test = pdb.write_renumber_pdb(head_idx_atom=head_atoms, tail_idx_atom=tail_atoms)
newfile=pattern+"_renumber.pdb"

# Assign residues
pdb_new = ReadPdbFormat(newfile)
pdb_new.assign_residues_chains(residueinfo_file)


print ("Job Done!!!!")



