import datetime
from collections import defaultdict
from utils.Logger import init_logger
from utils.parse_arguments import print_header, parse_arguments_check
from topology.Topology import Topology


# =============================================================================
def bonds_from_conect_pdb(pdbfname):

    dict_bonds = defaultdict(list)

    # Open file and find number of atoms and conect section
    natoms = 0
    listtotalbonds = []
    with open(pdbfname, 'r') as fpdb:
        lines = fpdb.readlines()
        for iline in lines:
            if iline.find("CONECT") != -1:
                ini = 6
                delta = 5
                iat = int(iline[ini:ini + delta]) - 1
                lbonds = []
                sentinel = ini + delta
                while True:
                    try:
                        jat = int(iline[sentinel:sentinel+delta]) - 1
                        lbonds.append(jat)
                        listtotalbonds.append(sorted([iat, jat]))
                    except ValueError:
                        break
                    sentinel += delta
                dict_bonds[iat] = sorted(lbonds)
            elif iline.find("ATOM") != -1 or iline.find("HETATM") != -1:
                natoms += 1

    # Remove repeat bonds in listtotalbonds
    # using set() + sorted()
    # removing duplicate sublist
    res = list(set(tuple(sorted(sub)) for sub in listtotalbonds))
    listtotalbonds = []
    for i in res:
        listtotalbonds.append(list(i))

    topo_pdb = Topology(natoms, listbonds=listtotalbonds)

    return topo_pdb


# =============================================================================
def bonds_from_conect_psf(psffname):

    listtotalbonds = []
    with open(psffname, 'r') as fpsf:
        lines = fpsf.readlines()
        for idx, iline in enumerate(lines):
            if iline.find("NATOM") != -1:
                natoms = int(iline.split()[0])
            elif iline.find("NBOND") != -1:
                nbonds = int(iline.split()[0])
                idxbond = idx
        idxbond += 1
        listtotalbonds = []
        sentinel = False
        while True:
            ibondline = lines[idxbond]
            ibond_list_len = int(len(ibondline.split())/2)
            ibond_list = [int(i) for i in ibondline.split()]
            for j in range(0, ibond_list_len+4, 2):
                try:
                    iat = ibond_list[j] - 1
                    jat = ibond_list[j+1] - 1
                    listtotalbonds.append(sorted([iat, jat]))
                except IndexError:
                    sentinel = True
            if sentinel:
                break
            idxbond += 1

    topo_psf = Topology(natoms, listbonds=listtotalbonds)

    return topo_psf


# =============================================================================
def main_app(version):

    # Init logger =============================================================
    starttime = datetime.datetime.now()
    logger = init_logger("Output", fileoutput="InfoCheckPDBPSF.log",
                         append=False, inscreen=True)
    print_header(version, logger)

    # Parse arguments in the command line =====================================
    opts = parse_arguments_check()

    # Comparing =====================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** COMPARING PDB and PSF files ({}) ****** \n".format(now)
    m += "\t\t PDB file:  {}\n".format(opts.pdbfile)
    m += "\t\t PSF file:  {}".format(opts.psffile)
    print(m) if logger is None else logger.info(m)

    topo_pdb = bonds_from_conect_pdb(opts.pdbfile)
    topo_psf = bonds_from_conect_psf(opts.psffile)

    # Check number of atoms
    equal_bonds = True
    if topo_pdb.natoms != topo_psf.natoms:
        m = "\t\t\t Warning!!! Number of atoms in PDB ({}) and PSF ({}) are different!!!!".\
            format(topo_pdb.natoms, topo_psf.natoms)
        print(m) if logger is None else logger.info(m)
        equal_bonds = False
    else:
        m = "\t\t\t Number of atoms in PDB ({}) and PSF ({}) are equal!!!!".\
            format(topo_pdb.natoms, topo_psf.natoms)
        print(m) if logger is None else logger.info(m)

    if len(topo_pdb._bonds) != len(topo_psf._bonds):
        m = "\t\t\t Warning!!! Number of bonds in PDB ({}) and PSF ({}) are different!!!!".\
            format(len(topo_pdb._bonds), len(topo_psf._bonds))
        print(m) if logger is None else logger.info(m)
        equal_bonds = False
    else:
        m = "\t\t\t Number of bonds in PDB ({}) and PSF ({}) are equal!!!!".\
            format(len(topo_pdb._bonds), len(topo_psf._bonds))
        print(m) if logger is None else logger.info(m)

    # Check bondlist ==============================================================
    m = "\n"
    ndiffpdb = 0
    for item in topo_pdb._bonds:
        if item not in topo_psf._bonds:
            m += "\t\t\t Warning!!! Bond {} in PDB is not in the PSF file\n".format(item)
            equal_bonds = False
            ndiffpdb += 1
    print(m) if logger is None else logger.info(m)

    m = ""
    ndiffpsf = 0
    for item in topo_psf._bonds:
        if item not in topo_pdb._bonds:
            m += "\t\t\t Warning!!! Bond {} in PSF is not in the PDB file\n".format(item)
            equal_bonds = False
            ndiffpsf += 1
    print(m) if logger is None else logger.info(m)

    if equal_bonds:
        m += "\t\t\t Bonds in PDB and PSF seem to be equal!!!!\n".format(item)
        print(m) if logger is None else logger.info(m)
    else:
        m += "\n\t\t\t Warning!!!! Number of different bonds in PDB = {}\n".format(ndiffpdb)
        m += "\t\t\t Warning!!!! Number of different bonds in PSF = {}\n".format(ndiffpsf)
        print(m) if logger is None else logger.info(m)

    # Final message ===========================================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    endtime = datetime.datetime.now()
    m = "\n\t\tJob  Done at {} ============\n".format(now)
    print(m) if logger is None else logger.info(m)
    m = "\t\tTotal time: {0:.2f} seconds".format((endtime-starttime).total_seconds())
    print(m) if logger is None else logger.info(m)


# =============================================================================
if __name__ == "__main__":

    __version__ = "1.1"
    main_app(__version__)
