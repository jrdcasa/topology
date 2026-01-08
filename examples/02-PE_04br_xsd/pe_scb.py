import datetime
from topology.readmol.readxsdformat import ReadXsdFormat
import utils

# Logger
filelog = "example02.log"
log = utils.init_logger("Output", fileoutput=filelog, append=False, inscreen=True)
m = "\n\t***************** EXAMPLE 02 PE-SCB with 4 branches *****************"
print(m) if log is None else log.info(m)
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\t\tStarting: \t {}\n".format(now))

# Create the xsd object
filename = "./01-Min_Frame_001.xsd"
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Reading {}...({})".format(filename, now)
xsd = ReadXsdFormat(filename)
print(m) if log is None else log.info(m)

# Write topology
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing xyz file...({})".format(now)
print(m) if log is None else log.info(m)
xsd.write_xyz(filename_xyz="PE300_04br04.xyz")

now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing gro file...({})".format(now)
print(m) if log is None else log.info(m)
xsd.write_gro(filename_gro="PE300_04br04.gro")

now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing pdb file...({})".format(now)
print(m) if log is None else log.info(m)
xsd.write_pdb(filename_pdb="PE300_04br04.pdb", separate_chains=True)

# Logger
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\t\tFinishing: \t {}\n".format(now))
m = "\t============== END   EXAMPLE02 ================================"
print(m) if log is None else log.info(m)
