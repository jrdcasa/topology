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


# =============================================================================
def parse_arguments_typing_carb():

    import time

    desc = """ Type a PDB containing a polysaccharide.\n"""

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-p", "--pdb", dest="pdbfile",
                        help="A PDB file contaning a polysaccharide for typing",
                        action="store", metavar="PDB_FILE", required=True)

    parser.add_argument("-b", "--box", dest="boxdimensions", nargs="+",
                        help="Box dimensions in the GRO format, "
                             "lx, ly, lz",
                        action="store", metavar="PDB_FILE", required=False, default=None)

    args = parser.parse_args()

    # Check for existing files:
    if not os.path.isfile(args.pdbfile):
        print(desc)
        time.sleep(.25)
        parser.error("PDB file must exist!!!! ({})".format(os.path.abspath(args.pdbfile)))

    # Check dimensions
    if args.boxdimensions is not None:
        if len(args.boxdimensions) not in [1, 3]:
            print("ERROR in box dimensions")
            exit()

    return args


# =============================================================================
def parse_arguments_type_gas():

    import time

    desc = """ Type a PDB containing a gas to perfomr MD in GROMACS and MC in RASPA.\n"""

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-p", "--pdb", dest="pdbfile",
                        help="A PDB file contaning a molecule for typing. "
                             "The PDB file must be edited manually: atom names should correspond to "
                             "the atom types defined in the force field, and the residue name field "
                             "must be completed as well.",
                        action="store", metavar="PDB_FILE", required=True)

    parser.add_argument("-f", "--ffname", dest="ffname", choices=["OPLSAA"],
                        help="Name of the force field",
                        action="store", metavar="FFNAME", required=True)

    parser.add_argument("-c", "--charges", dest="chargefile",
                        help="List of atomic charges, provided in the same order as the atoms appear in the PDB file.",
                        action="store", metavar="CHARGE_FILE", required=True)

    parser.add_argument("-b", "--box", dest="boxdimensions", nargs="+",
                        help="Box dimensions in the GRO format, "
                             "lx, ly, lz",
                        action="store", metavar="BOX_FILE", required=False, default=None)

    parser.add_argument("--impropers", dest="impropers",
                        help="A file with improper angles and parameters in GROMACS TOP format.",
                        action="store", metavar="BOX_FILE", required=False, default=None)

    parser.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="InfoTypeGas.log")

    args = parser.parse_args()

    # Check for existing files:
    if not os.path.isfile(args.pdbfile):
        print(desc)
        time.sleep(.25)
        parser.error("PDB file must exist!!!! ({})".format(os.path.abspath(args.pdbfile)))

    if not os.path.isfile(args.chargefile):
        print(desc)
        time.sleep(.25)
        parser.error("CHARGE file must exist!!!! ({})".format(os.path.abspath(args.chargefile)))

    # Check dimensions
    if args.boxdimensions is not None:
        if len(args.boxdimensions) not in [1, 3]:
            print("ERROR in box dimensions")
            exit()

    return args


# =============================================================================
def parse_arguments_multiplepdb():

    import time

    desc = """ Build a XTC trajectory from PDB.\n"""

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-d", "--dir", dest="pdbdir",
                        help="A directory containing the PDB files.",
                        action="store", metavar="PDB_FILE", required=True)

    parser.add_argument("-r", "--renumber_pdb", dest="renumberpdb",
                        help=textwrap.dedent("Renumber pdb in accordance with the head and tail info"
                                             " in the provided file. The format is:\n"
                                             "  nmols: 40\n"
                                             "  #ich idx-head idx-tail (indexes start at 0)\n"
                                             "  0 1 456\n"
                                             "  1 465 920\n"
                                             "  ...\n"
                                             "  39 18097 18552\n"),
                        action="store", metavar="HEAD_TAIL_file", required=True
                        )

    parser.add_argument("-a", "--assign_residues", dest="assignresidues",
                        help=textwrap.dedent("Assign residues in accordance with residue file."),
                        action="store", metavar="SETUP_RESIDUE_file", required=False)

    parser.add_argument("-n", "--ncpus", dest="ncpus", type=int,
                        help="Number of cpus to be used\n",
                        action="store", metavar="NCPUS", required=False
                        )

    parser.add_argument("-p", "--pattern", dest="pattern",
                        help="A string pattern to name the new files",
                        action="store", metavar="STR_PATTERN", required=False, default="trajectory")

    args = parser.parse_args()

    # Check for existing files:
    if not os.path.isdir(args.pdbdir):
        print(desc)
        time.sleep(.25)
        parser.error("Directory with PDB files must exist!!!! ({})".format(os.path.abspath(args.pdbdir)))

    return args


# =============================================================================
def print_header_type_cb(version, logger_log=None):
    msg = """
    ***********************************************************************
                     Typing carbohydrates to run MD
              ----------------------------------------------

                                Version {}

                              Dr. Javier Ramos
                      Macromolecular Physics Department
                Instituto de Estructura de la Materia (IEM-CSIC)
                               Madrid (Spain)

        This program processes a PDB file containing a carbohydrate generated 
        by either CarbBuilder (https://people.cs.uct.ac.za/~mkuttel/Downloads.html) 
        or DoGlycans (https://bitbucket.org/biophys-uh/doglycans/src/main/). 
        It attempts to type the molecule using the CHARMM36 force field for sugars.

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)


# =============================================================================
def print_header_multiplepdb(version, logger_log=None):
    msg = """
    ***********************************************************************
                     Typing carbohydrates to run MD
              ----------------------------------------------

                                Version {}

                              Dr. Javier Ramos
                      Macromolecular Physics Department
                Instituto de Estructura de la Materia (IEM-CSIC)
                               Madrid (Spain)

        This program processes a set of PDB files, typically generated 
        by Materials Studio, representing a molecular trajectory and 
        converts them into GRO and XTC formats compatible with GROMACS.

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)


# =============================================================================
def print_header_type_gas(version, logger_log=None):
    msg = """
    ***********************************************************************
                  Typing small molecules to run MD and MC
              ----------------------------------------------

                                Version {}

                              Dr. Javier Ramos
                      Macromolecular Physics Department
                Instituto de Estructura de la Materia (IEM-CSIC)
                               Madrid (Spain)

        This program processes a set of PDB files, typically generated 
        by Materials Studio, representing a molecular trajectory and 
        converts them into GRO and XTC formats compatible with GROMACS.

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
    
    The PDB file must contain the atom type in the atom name field and 
    the residue name in the residue name field. Once the PDB has been generated, 
    it must be properly formatted to ensure compatibility. 
    In particular, each atom type must match a valid atom type 
    defined in the force field being used.
    It is essential to verify that the file is correctly formatted prior to use.
    
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)
