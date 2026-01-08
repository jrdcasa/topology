import datetime
import os.path

from utils.Logger import init_logger
from topology.readmol.readpdbformat import ReadPdbFormat
from utils.parse_arguments import print_header, parse_arguments_label, command_info


# =============================================================================
def main_app(version=None):

    # Init logger =============================================================
    starttime = datetime.datetime.now()
    logger = init_logger("Output", fileoutput="InfoLabelPDB.log",
                         append=False, inscreen=True)
    print_header(version, logger)

    # Parse arguments in the command line =====================================
    opts = parse_arguments_label()
    # Print command info
    command_info(opts, logger=logger)

    # Comparing =====================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** LABELLING PDB file ({}) ****** \n".format(now)
    m += "\t\t PDB file:  {}\n".format(opts.pdbfile)
    print(m) if logger is None else logger.info(m)

    # Check the type of input file and read it ================================
    with open(opts.pdbfile) as f:
        lines = str(f.readlines())
        if lines.count("CONECT") != 0:
            isconect = True
        else:
            isconect = False

    pdbobj = ReadPdbFormat(opts.pdbfile, isconect=isconect, logger=logger)
    head_atoms, tail_atoms = pdbobj.read_head_tail_info(opts.headfile)
    filenamepdb_labeled = os.path.splitext(os.path.split(opts.pdbfile)[-1])[0]+"_labeled.pdb"
    if opts.assignresidues is not None:
        residueinfo_file = opts.assignresidues
    else:
        residueinfo_file = None

    pdbobj.write_label_pdb(head_idx_atom=head_atoms, tail_idx_atom=tail_atoms,
                           assign_residues_info=residueinfo_file, fnameout=filenamepdb_labeled,
                           isunwrap=False)

    # Get files for analysis ================================
    with open(filenamepdb_labeled) as f:
        lines = str(f.readlines())
        if lines.count("CONECT") != 0:
            isconect = True
        else:
            isconect = False
    pdbobj = ReadPdbFormat(filenamepdb_labeled, isconect=isconect, logger=logger)
    pdbobj.read_head_tail_info(opts.headfile)
    pdbobj.get_backbone_headtail_files_for_analysis()
    filename_psf = os.path.splitext(filenamepdb_labeled)[0]+".psf"
    pdbobj.write_psf(filename_psf=filename_psf, improper_angles=None)

    # Final message ===========================================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    endtime = datetime.datetime.now()
    m = "\n\t\tJob  Done at {} ============\n".format(now)
    print(m) if logger is None else logger.info(m)
    m = "\t\tTotal time: {0:.2f} seconds".format((endtime-starttime).total_seconds())
    print(m) if logger is None else logger.info(m)


# =============================================================================
if __name__ == "__main__":

    __version__ = "0.1"
    main_app(__version__)
