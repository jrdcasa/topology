import argparse
import textwrap
import os
import datetime
import sys


# =============================================================================
def print_header(version, logger_log=None):
    msg = """
    ***********************************************************************
                     Manipulate topology for polymers
              ----------------------------------------------

                                Version {}

                              Dr. Javier Ramos
                      Macromolecular Physics Department
                Instituto de Estructura de la Materia (IEM-CSIC)
                               Madrid (Spain)

        Topology is an open-source python library to quickly modify topology 
        of polymer molecules

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)


# =============================================================================
def command_info(opts, logger=None):

    m1 = ""
    for item in sys.argv[1:]:
        m1 += " {}".format(item)
    m = "\n\t\tCommand line: \n"
    m += "\t\t\tpython {}".format(os.path.split(sys.argv[0])[1])
    m += m1+"\n"
    m += "\t\t\t         or\n"
    m += "\t\t\t{} ".format(os.path.split(sys.argv[0])[1])
    m += m1+"\n"
    print(m) if logger is None else logger.info(m)


# =============================================================================
def parse_arguments():

    import time

    desc = """ Modify topology of a polymer or molecule.\n"""

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    # # group1 = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-i", "--inp", dest="inputfile",
                        help="A file containing the input information. "
                             "It should be a XSD, PDB or MOL2 file",
                        action="store", metavar="INPUT_FILE", required=True)
    parser.add_argument("-p", "--pattern", dest="pattern",
                        help="A string pattern to name the new files",
                        action="store", metavar="STR_PATTERN", required=False, default="topology")

    parser.add_argument("-r", "--renumber_pdb", dest="renumberpdb",
                        help=textwrap.dedent("Renumber pdb in accordance with the head and tail info"
                                             " in the provided file. The format is:\n"
                                             "  nmols: 40\n"
                                             "  #ich idx-head idx-tail (indexes start at 0)\n"
                                             "  0 1 456\n"
                                             "  1 465 920\n"
                                             "  ...\n"
                                             "  39 18097 18552\n"),
                        action="store", metavar="HEAD_TAIL_file", required=False)

    parser.add_argument("-a", "--assign_residues", dest="assignresidues",
                        help=textwrap.dedent("Assign residues in accordance with residue file."),
                        action="store", metavar="SETUP_RESIDUE_file", required=False)

    parser.add_argument("--filemap", dest="filemap",
                        help=textwrap.dedent("A file to match LAMMPS type with name and element. The format is:\n"
                                             "4\n"
                                             "# # TypeLammps Name Element\n" 
                                             "#   1 CA C\n"
                                             "#   2 HA H\n"
                                             "#   3 CB C\n"
                                             "#   4 HB H\n"),
                        action="store", metavar="FILEMAP", required=False)

    parser.add_argument("--separate_chains", dest="separatechains",
                        help="It creates a pdb file for each chain.",
                        action="store_true", required=False)

    parser.add_argument("-w", "--isunwrap", dest="isunwrap",
                        help="If true, the coordinates are unwrapped in the final structure.",
                        action="store_true", required=False, default=False)

    parser.add_argument("--guess_improper", dest="guessimproper",
                        help="Try to guess the improper angles in the system.",
                        action="store_true", required=False, default=False)

    args = parser.parse_args()

    #
    # Check for existing files:
    if not os.path.isfile(args.inputfile):
        print(desc)
        time.sleep(.25)
        parser.error("Input file must exist!!!! ({})".format(os.path.abspath(args.inputfile)))
    if args.renumberpdb is not None and not os.path.isfile(args.renumberpdb):
        print(desc)
        time.sleep(.25)
        parser.error("Head-Tail info file must exist!!!! ({})".format(os.path.abspath(args.renumberpdb)))
    if args.assignresidues is not None and not os.path.isfile(args.assignresidues):
        print(desc)
        time.sleep(.25)
        parser.error("Residues info file must exist!!!! ({})".format(os.path.abspath(args.assignresidues)))

    return args


# =============================================================================
def parse_arguments_check():

    import time

    desc = """ Check if the CONECT of a PDB and NBONDS in a PSF are the same.\n"""

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    # # group1 = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-p", "--pdb", dest="pdbfile",
                        help="A PDB file containing CONECT section.",
                        action="store", metavar="PDB_FILE", required=True)
    parser.add_argument("-t", "--psf", dest="psffile",
                        help="A PSF file containing a NBONDS section.",
                        action="store", metavar="PSF_FILE", required=True)

    args = parser.parse_args()

    # Check for existing files:
    if not os.path.isfile(args.pdbfile):
        print(desc)
        time.sleep(.25)
        parser.error("PDB file must exist!!!! ({})".format(os.path.abspath(args.pdbfile)))
    if not os.path.isfile(args.psffile):
        print(desc)
        time.sleep(.25)
        parser.error("PSF file must exist!!!! ({})".format(os.path.abspath(args.psffile)))

    return args


# =============================================================================
def parse_arguments_label():

    import time

    desc = """ Label a PDB file, designating head and tail atoms as well as backbone atoms.
     Additionally, generate a list of atoms  \n"""

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-p", "--pdb", dest="pdbfile",
                        help="A PDB file to be labelled",
                        action="store", metavar="PDB_FILE", required=True)

    parser.add_argument("-f", "--headfile", dest="headfile",
                        help=textwrap.dedent("Renumber pdb in accordance with the head and tail info"
                                             " in the provided file. The format is:\n"
                                             "  nmols: 40\n"
                                             "  #ich idx-head idx-tail (indexes start at 0)\n"
                                             "  0 1 456\n"
                                             "  1 465 920\n"
                                             "  ...\n"
                                             "  39 18097 18552\n"),
                        action="store", metavar="HEAD_TAIL_file", required=True)

    parser.add_argument("-a", "--assign_residues", dest="assignresidues",
                        help=textwrap.dedent("Assign residues in accordance with residue file."),
                        action="store", metavar="SETUP_RESIDUE_file", required=False)

    args = parser.parse_args()

    # Check for existing files:
    if not os.path.isfile(args.pdbfile):
        print(desc)
        time.sleep(.25)
        parser.error("PDB file must exist!!!! ({})".format(os.path.abspath(args.pdbfile)))
    if not os.path.isfile(args.headfile):
        print(desc)
        time.sleep(.25)
        parser.error("Head-Tail file must exist!!!! ({})".format(os.path.abspath(args.headfile)))

    return args
