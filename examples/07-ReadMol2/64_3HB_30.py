import datetime
from topology.readmol.readmol2format import ReadMol2Format
import utils

filename = "./64_3HB_alfa_30.mol2"
pattern = "64_3HB_alfa_30"

# Logger
filelog = "64_3HB_alfa_30.log"
log = utils.init_logger("Output", fileoutput=filelog, append=False, inscreen=True)
m = "\n\t***************** P3HB i64 chain 30 monomers *****************"
print(m) if log is None else log.info(m)
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\t\tStarting: \t {}\n".format(now))

# Create the xsd object
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Reading {}...({})".format(filename, now)
mol2 = ReadMol2Format(filename)
print(m) if log is None else log.info(m)

# Write topology
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing xyz file...({})".format(now)
print(m) if log is None else log.info(m)
mol2.write_xyz(filename_xyz="{}.xyz".format(pattern))

now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing gro file...({})".format(now)
print(m) if log is None else log.info(m)
mol2.write_gro(filename_gro="{}.gro".format(pattern))

now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing pdb file...({})".format(now)
print(m) if log is None else log.info(m)
xsd.write_pdb(filename_pdb="{}.pdb".format(pattern), separate_chains=True)

# Logger
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\t\tFinishing: \t {}\n".format(now))
m = "\t============== END   ==============================="
print(m) if log is None else log.info(m)

