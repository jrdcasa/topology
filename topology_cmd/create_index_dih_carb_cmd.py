import datetime
import argparse
import os
import re
from collections import defaultdict
from utils.Logger import init_logger
from utils.parse_arguments import command_info
import MDAnalysis as mda
import numpy as np


# =============================================================================
def print_header(version, logger_log=None):

    msg = """
    ***********************************************************************
                     Create the index file for phi-psi 
                           angles in carbohydrates
              ----------------------------------------------

                                Version {}

                              Dr. Javier Ramos
                      Macromolecular Physics Department
                Instituto de Estructura de la Materia (IEM-CSIC)
                               Madrid (Spain)

        This program generates an index file listing the atomic indices that 
        define dihedral angles in carbohydrates, using the same format 
        as GROMACS

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)


# =============================================================================
def parse_arguments():

    import time

    desc = """ Create an index file for link dihedrals in a polysaccharide.\n"""

    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--inp", dest="inputfile",
                        help="A file containing the coordinates information."
                             "It should be a GRO, PDB or MOL2 file",
                        action="store", metavar="INPUT_FILE", required=True)

    parser.add_argument("-t", "--topo", dest="topologyfile",
                        help="A file containing the topology information."
                             "It should be a PSF or TPR file",
                        action="store", metavar="TOPOLOGY_FILE", required=True)

    parser.add_argument("--links", dest="linklist", nargs="+",
                        help="A list of the type of links in the polysaccharide. \n"
                             "For example: 1-2 1-3 for polysaccharide with 1-2 and 1-3 links.",
                        action="store", metavar="LINK_LIST", required=True)

    args = parser.parse_args()

    # Check for existing files:
    if not os.path.isfile(args.inputfile):
        print(desc)
        time.sleep(.25)
        parser.error("Input file must exist!!!! ({})".format(os.path.abspath(args.inputfile)))

    if not os.path.isfile(args.topologyfile):
        print(desc)
        time.sleep(.25)
        parser.error("Topology file must exist!!!! ({})".format(os.path.abspath(args.topologyfile)))

    for item in args.linklist:
        if not bool(re.fullmatch(r'\d+-\d+', item)):
            print(desc)
            time.sleep(.25)
            parser.error("Link list must be [1-9]-[1-9] format not {}".format(item))

    return args


# =============================================================================
def main_app(version=None):

    # Init logger =============================================================
    starttime = datetime.datetime.now()
    dict_psi = defaultdict()
    dict_phi = defaultdict()
    logger = init_logger("Output", fileoutput="InfoCreateIndex.log",
                         append=False, inscreen=True)
    print_header(version, logger)

    # Parse arguments in the command line =====================================
    opts = parse_arguments()
    # Print command info
    command_info(opts, logger=logger)

    # Create an universe
    msg = "\t\tMDAnalysis version: {}".format(mda.__version__)
    print(msg) if logger is None else logger.info(msg)
    u_mda = mda.Universe(opts.topologyfile, opts.inputfile)

    # Links
    for item in opts.linklist:
        dict_psi[item] = []
        dict_phi[item] = []

    number_links = 0
    for item in opts.linklist:

        jlink = int(item.split("-")[0])
        klink = int(item.split("-")[1])
        jname = "C{0:1d}".format(jlink)
        kname = "O{0:1d}".format(klink)
        lname = "C{0:1d}".format(klink)
        if jname == "C1":
            # ilink = 5
            iname = "O5"
        else:
            iname = ""
            msg = "\n\t\tERROR: {} is not implemented".format(jname)
            msg += "\n\t\tERROR: iname-{}-{}-{} cannot be found".format(jname, kname, lname)
            print(msg) if logger is None else logger.error(msg)
            exit()
        if lname == "C3":
            # mlink = 2
            mname = "C2"
        elif lname == "C4":
            # mlink = 3
            mname = "C3"
        else:
            mname = ""
            msg = "\n\t\tERROR: {} is not implemented".format(lname)
            msg += "\n\t\tERROR: {}-{}-{}-{} cannot be found".format(iname, jname, kname, lname)
            print(msg) if logger is None else logger.error(msg)
            exit()

        sel_ini = u_mda.select_atoms("name {}".format(jname))
        #
        # For each iname atom
        for jatom in sel_ini:
            # idx-jdx-kdx-ldx
            # jdx-kdx-ldx-mdx
            jdx = jatom.id + 1

            # iatom
            try:
                ii = np.where(jatom.bonded_atoms.names == iname)[0][0]
            except IndexError:
                continue
            idx = jatom.bonded_atoms.indices[ii] + 1
            # iatom = u_mda.select_atoms("index {}".format(idx-1))[0]

            # katom
            try:
                ii = np.where(jatom.bonded_atoms.names == kname)[0][0]
            except IndexError:
                continue
            kdx = jatom.bonded_atoms.indices[ii] + 1
            katom = u_mda.select_atoms("index {}".format(kdx-1))[0]

            # latom
            try:
                ii = np.where(katom.bonded_atoms.names == lname)[0][0]
            except IndexError:
                msg = "\n\t\tERROR: {} cannot be found".format(lname)
                msg += "\n\t\tERROR: {}-{} cannot be found".format(jname, kname)
                print(msg) if logger is None else logger.error(msg)
                exit()
            ldx = katom.bonded_atoms.indices[ii] + 1
            latom = u_mda.select_atoms("index {}".format(ldx-1))[0]

            # matom
            try:
                ii = np.where(latom.bonded_atoms.names == mname)[0][0]
            except IndexError:
                msg = "\n\t\tERROR: {} cannot be found".format(mname)
                msg += "\n\t\tERROR: {}-{} cannot be found".format(jname, kname)
                print(msg) if logger is None else logger.error(msg)
                exit()
            mdx = latom.bonded_atoms.indices[ii] + 1

            label_phi = "{}-{}-{}-{}".format(iname, jname, kname, lname)
            label_psi = "{}-{}-{}-{}".format(jname, kname, lname, mname)
            dict_phi[item].append([idx, jdx, kdx, ldx, label_phi])
            dict_psi[item].append([jdx, kdx, ldx, mdx, label_psi])

            # DEBUG
            # print("psi:", idx, jdx, kdx, ldx, "|||", "{}-{}-{}-{}".format(iname, jname, kname, lname))
            # print("phi:", jdx, kdx, ldx, mdx, "|||", "{}-{}-{}-{}".format(jname, kname, lname, mname))
            # print("===========")
            number_links += 1

    # Write index files
    line_phi = "[ phi ]\n"
    line_phi_test = "[ phi ]\n"
    line_psi = "[ psi ]\n"
    line_psi_test = "[ psi ]\n"
    line_phipsi = " [ phipsi ]\n"
    line_phipsi_test = " [ phipsi ]\n"
    for key, values in dict_phi.items():
        for item in values:
            line_phi += "{} {} {} {}\n".format(item[0], item[1], item[2], item[3])
            line_phi_test += "{} {} {} {} #{}\n".format(item[0], item[1], item[2], item[3], item[4])
            line_phipsi += "{} {} {} {}\n".format(item[0], item[1], item[2], item[3])
            line_phipsi_test += "{} {} {} {} #{}\n".format(item[0], item[1], item[2], item[3], item[4])

    for key, values in dict_psi.items():
        for item in values:
            line_psi += "{} {} {} {}\n".format(item[0], item[1], item[2], item[3])
            line_psi_test += "{} {} {} {} #{}\n".format(item[0], item[1], item[2], item[3], item[4])
            line_phipsi += "{} {} {} {}\n".format(item[0], item[1], item[2], item[3])
            line_phipsi_test += "{} {} {} {} #{}\n".format(item[0], item[1], item[2], item[3], item[4])

    with open('index.ndx', 'w') as f1:
        f1.writelines(line_phi)
        f1.writelines(line_psi)
        f1.writelines(line_phipsi)
        for key, values in dict_phi.items():
            f1.writelines("[ phi{} ]\n".format(key))
            for item in values:
                f1.writelines("{0} {1} {2} {3}\n".format(item[0], item[1], item[2], item[3]))
        for key, values in dict_psi.items():
            f1.writelines("[ psi{} ]\n".format(key))
            for item in values:
                f1.writelines("{0} {1} {2} {3}\n".format(item[0], item[1], item[2], item[3]))

    with open('index_check.ndx', 'w') as f2:
        f2.writelines(line_phi_test)
        f2.writelines(line_psi_test)
        f2.writelines(line_phipsi_test)
        for key, values in dict_phi.items():
            f2.writelines("[ phi{} ]\n".format(key))
            for item in values:
                f2.writelines("{0} {1} {2} {3} {4}\n".format(item[0], item[1], item[2], item[3], item[4]))
        for key, values in dict_psi.items():
            f2.writelines("[ psi{} ]\n".format(key))
            for item in values:
                f2.writelines("{0} {1} {2} {3} {4}\n".format(item[0], item[1], item[2], item[3], item[4]))

    msg = "\n\t\t============= SUMMARY =============\n"
    msg += "\t\t Number of links             = {}\n".format(number_links)
    for key, values in dict_phi.items():
        msg += "\t\t Number of {} phi dihedrals = {}\n".format(key, len(values))
    for key, values in dict_psi.items():
        msg += "\t\t Number of {} psi dihedrals = {}\n".format(key, len(values))
    print(msg) if logger is None else logger.info(msg)

    msg = "\n\t\t Indices have been written in index.ndx file\n"
    msg += "\t\t Indices to check have been written in index_check.ndx file\n"
    print(msg) if logger is None else logger.info(msg)

    # Final message ===========================================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    endtime = datetime.datetime.now()
    m = "\n\t\tJob  Done at {} ============".format(now)
    print(m) if logger is None else logger.info(m)
    m = "\t\tTotal time: {0:.2f} seconds".format((endtime-starttime).total_seconds())
    print(m) if logger is None else logger.info(m)


# =============================================================================
if __name__ == "__main__":

    version = "1.1"
    main_app(version)
