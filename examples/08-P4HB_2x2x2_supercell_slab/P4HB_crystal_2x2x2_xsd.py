import datetime
from topology.readmol.readxsdformat import ReadXsdFormat
import utils

filename = "./P4HB_crystal_2x2x2.xsd"
pattern = "P4HB_crystal_2x2x2_xsd"

# Logger
filelog = "P4HB_crystal_2x2x2_xsd.log"
log = utils.init_logger("Output", fileoutput=filelog, append=False, inscreen=True)
m = "\n\t***************** P4HB 8 chains X monomers *****************"
print(m) if log is None else log.info(m)
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\t\tStarting: \t {}\n".format(now))

# Create the xsd object
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Reading {}...({})".format(filename, now)
xsd = ReadXsdFormat(filename)
print(m) if log is None else log.info(m)

# Write topology
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing xyz file...({})".format(now)
print(m) if log is None else log.info(m)
xsd.write_xyz(filename_xyz="{}.xyz".format(pattern))

now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing gro file...({})".format(now)
print(m) if log is None else log.info(m)
xsd.write_gro(filename_gro="{}.gro".format(pattern))

# now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
# m = "\t\t Writing pdb file...({})".format(now)
# print(m) if log is None else log.info(m)
# xsd.write_pdb(filename_pdb="{}.pdb".format(pattern), separate_chains=True)

# Logger
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\t\tFinishing: \t {}\n".format(now))
m = "\t============== END   ==============================="
print(m) if log is None else log.info(m)

