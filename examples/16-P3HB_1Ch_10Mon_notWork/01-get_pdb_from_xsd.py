import datetime
from topology.readmol.readxsdformat import ReadXsdFormat

filename = "./P3HB_10mon_01c.xsd"
pattern = "P3HB_10mon_01c_xsd"

# Create the xsd object
xsd = ReadXsdFormat(filename)

# Write pdb file
xsd.write_pdb(filename_pdb="{}.pdb".format(pattern), separate_chains=False)

# Logger
print("Job Done!!!")

