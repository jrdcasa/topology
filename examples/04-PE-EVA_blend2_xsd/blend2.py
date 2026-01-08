import datetime
from topology.readmol.readxsdformat import ReadXsdFormat
import utils

# Logger
filelog = "example04    .log"
log = utils.init_logger("Output", fileoutput=filelog, append=False, inscreen=True)
m = "\n\t***************** EXAMPLE 04 PE+EVA no consecutives *****************"
print(m) if log is None else log.info(m)
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\t\tStarting: \t {}\n".format(now))

# Create the xsd object
filename = "./blend2.xsd"
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Reading {}...({})".format(filename, now)
xsd = ReadXsdFormat(filename)
print(m) if log is None else log.info(m)

# Write topology
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing xyz file...({})".format(now)
print(m) if log is None else log.info(m)
xsd.write_xyz(filename_xyz="blend2.xyz")

now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing gro file...({})".format(now)
print(m) if log is None else log.info(m)
xsd.write_gro(filename_gro="blend2.gro")

now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
m = "\t\t Writing pdb file...({})".format(now)
print(m) if log is None else log.info(m)
xsd.write_pdb(filename_pdb="blend2.pdb", separate_chains=True)

# Logger
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
log.info("\t\tFinishing: \t {}\n".format(now))
m = "\t============== END   EXAMPLE04 ================================"
print(m) if log is None else log.info(m)
