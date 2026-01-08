import datetime
from topology.readmol.readxsdformat import ReadXsdFormat
from topology.readmol.readpdbformat import ReadPdbFormat
from topology.readmol.readdataLAMMPSformat import ReadDataLammpsFormat
import utils

filename = "./restart.data.500000.data"
filemap = "./filemap.dat"
pattern="./lammps_read"
#headinfo_file = "./butane_headtail.dat"
#residueinfo_file = "./butane_residues.dat"

# Logger
filelog = pattern+".log"
log = utils.init_logger("Output", fileoutput=filelog, append=False, inscreen=True)
m = "\n\t***************** PE C102 40Chains LOPLS *****************"
print(m) if log is None else log.info(m)
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\t\tStarting: \t {}\n".format(now))

# Create the xsd object
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Reading {}...({})".format(filename, now)
print(m) if log is None else log.info(m)
xsd = ReadDataLammpsFormat(filename, filemap)
m = "\t\t Writing PSF {}...({})".format(filename, now)
print(m) if log is None else log.info(m)
xsd.write_psf()
m = "\t\t Writing PDB {}...({})".format(filename, now)
print(m) if log is None else log.info(m)
xsd.write_pdb()

# # Write pdb file
# now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
# m = "\t\t Writing pdb file...({})".format(now)
# print(m) if log is None else log.info(m)
# filenamepdb = "{}.pdb".format(pattern)
# xsd.write_pdb(filename_pdb=filenamepdb, separate_chains=False)
#
# # Renumber PDB
# now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
# m = "\t\t Read PDB and renumber pdb file...({})".format(now)
# print(m) if log is None else log.info(m)
# pdb = ReadPdbFormat(filenamepdb)
# head_atoms, tail_atoms = pdb.read_head_tail_info(headinfo_file)
# test = pdb.write_renumber_pdb(head_idx_atom=head_atoms, tail_idx_atom=tail_atoms)
#
# # Assign residues PDB
# filenamepdb=pattern+"_renumber.pdb"
# now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
# m = "\t\t Read PDB and assign residues...({})".format(now)
# print(m) if log is None else log.info(m)
# pdb_new = ReadPdbFormat(filenamepdb)
# pdb_new.assign_residues_chains(residueinfo_file)


# Logger
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\n\t\tFinishing: \t {}\n".format(now))
m = "\t============== END   ==============================="
print(m) if log is None else log.info(m)

