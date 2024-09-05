import datetime
import sys
import os
from utils.Logger import init_logger
from topology.readmol.readpdbformat import ReadPdbFormat

#pdbfile = "Helix6_capped.pdb"
pdbfile="3Helix12_link.pdb"
isconect = True

# Init logger =============================================================
starttime = datetime.datetime.now()
logger = init_logger("Output", fileoutput="assign_name_helix.log",
                     append=False, inscreen=True)
m1 = ""
for item in sys.argv[1:]:
    m1 += " {}".format(item)
m = "\n\t\tCommand line: \n"
m += "\t\t\tpython {}".format(os.path.split(sys.argv[0])[1])
m += m1 + "\n"
m += "\t\t\t         or\n"
m += "\t\t\t{} ".format(os.path.split(sys.argv[0])[1])
m += m1 + "\n"
print(m) if logger is None else logger.info(m)

pdbobj = ReadPdbFormat(pdbfile, isconect=isconect, logger=logger)

# Loop over all atoms in topology
for idx_at in range(pdbobj._topology.natoms):
    neigh_list_idx = pdbobj._topology._graphdict[idx_at]
    if pdbobj._topology.elements[idx_at] == "H":
        if pdbobj._topology._names[idx_at] == "H61" or\
           pdbobj._topology._names[idx_at] == "H62":
            continue
        for ineigh_idx in neigh_list_idx:
            if pdbobj._topology.elements[ineigh_idx] == "O":
                try:
                    n = pdbobj._topology._names[ineigh_idx][1]
                    pdbobj._topology._names[idx_at] = "H{}O".format(n)
                except IndexError:
                    print("ERROR atom {}".format(idx_at))
                    print("ERROR ineigh {}".format(ineigh_idx))
                    exit()
            else:
                try:
                    n = pdbobj._topology._names[ineigh_idx][1]
                    pdbobj._topology._names[idx_at] = "H{}".format(n)
                except IndexError:
                    print("ERROR atom {}".format(idx_at))
                    print("ERROR ineigh {}".format(ineigh_idx))
                    exit()

pdbobj.write_pdb("Helix6_capped_opt.pdb", separate_chains=False)

m = "\n\t\tPDB Created: {}\n".format("Helix6_capped_opt.pdb")
print(m) if logger is None else logger.info(m)


# Final message ===========================================================
now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
endtime = datetime.datetime.now()
m = "\n\t\tJob  Done at {} ============\n".format(now)
print(m) if logger is None else logger.info(m)
m = "\t\tTotal time: {0:.2f} seconds".format((endtime - starttime).total_seconds())
print(m) if logger is None else logger.info(m)