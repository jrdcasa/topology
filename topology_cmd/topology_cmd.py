import datetime
import os.path
import numpy as np
from utils.Logger import init_logger
from utils.parse_arguments import print_header, parse_arguments, command_info
from topology.readmol.readxsdformat import ReadXsdFormat
from topology.readmol.readpdbformat import ReadPdbFormat
from topology.readmol.readmol2format import ReadMol2Format
from topology.readmol.readdataLAMMPSformat import ReadDataLammpsFormat


# =============================================================================
def summary_initial_log(topobj, log=None):

    # Print info to log
    nows = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    mm1 = "\n\t\tSummarize of the input system ({})\n".format(nows)
    mm = "\t\t" + len(mm1) * "=" + "\n"
    mm += "\t\tInput file                : {}\n".format(topobj._fnamepath)
    mm += "\t\tNumber of molecules       : {}\n".format(topobj._nmols)
    mm += "\t\tNumber of residues        : {}\n".format(topobj._nres)
    mm += "\t\tNumber of atoms           : {}\n".format(topobj._natoms)
    mm += "\t\tNumber of bonds           : {}\n".format(topobj._nbonds)
    if len(topobj._topology._bonds) > 0:
        mm += "\t\tBond info available       : True\n"
    else:
        mm += "\t\tBond info available       : False\n"
    if topobj._isthere_boxdimension:
        mm += "\t\tSimulation box dimension  : {0:.4f} {1:.4f} {2:.4f} nm\n".format(float(topobj._boxlength[0]/10.),
                                                                                    float(topobj._boxlength[1]/10.),
                                                                                    float(topobj._boxlength[2]/10.))
        mm += "\t\tSimulation box angles     : {0:.2f} {1:.2f} {2:.2f} degrees\n".\
            format(float(topobj._boxangle[0]*180/np.pi), float(topobj._boxangle[1]*180/np.pi),
                   float(topobj._boxangle[2]*180/np.pi))
    else:
        mm += "\t\tSimulation box dimension  : {0:.4f} {1:.4f} {2:.4f} nm\n".format(0.0, 0.0, 0.0)
        mm += "\t\tSimulation box angles     : {0:.2f} {1:.2f} {2:.2f} degrees\n".format(0.0, 0.0, 0.0)
    mm += "\t\t" + len(mm1) * "=" + "\n"
    print(mm1 + mm) if log is None else log.info(mm1 + mm)


# =============================================================================
def main_app(version):

    # Init logger =============================================================
    starttime = datetime.datetime.now()
    logger = init_logger("Output", fileoutput="InfoTopology.log",
                         append=False, inscreen=True)
    print_header(version, logger)

    # Parse arguments in the command line =====================================
    opts = parse_arguments()
    # Print command info
    command_info(opts, logger=logger)

    # Check the type of input file and read it ================================
    ext_inputfile = str.upper(opts.inputfile.split(".")[-1])
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    logger.info("\t\tStarting: \t {}\n".format(now))
    if ext_inputfile == "XSD":
        m = "\t\tInput is a XSD format.\n"
        m += "\t\tReading... ({})".format(now)
        print(m) if logger is None else logger.info(m)
        obj = ReadXsdFormat(opts.inputfile)
    elif ext_inputfile == "PDB":
        m = "\t\tInput is a PDB format.\n"
        m += "\t\tReading... ({})".format(now)
        print(m) if logger is None else logger.info(m)
        obj = ReadPdbFormat(opts.inputfile)
    elif ext_inputfile == "MOL2":
        m = "\t\tInput is a MOL2 format.\n"
        m += "\t\tReading... ({})".format(now)
        print(m) if logger is None else logger.info(m)
        obj = ReadMol2Format(opts.inputfile)
    elif ext_inputfile == "DATA":
        m = "\t\tInput is a LAMMPS DATA format.\n"
        m += "\t\tReading... ({})".format(now)
        print(m) if logger is None else logger.info(m)
        if opts.filemap is None:
            m = "\t\tA filemap is needed for DATA-LAMMPS format.\n"
            print(m) if logger is None else logger.error(m)
            exit()
        obj = ReadDataLammpsFormat(opts.inputfile, opts.filemap)
    else:
        obj = None
        m = "Extension {} cannot be used. " \
            "Allowed formats are XSD, PDB, MOL2 and LAMMPS_data.".format(ext_inputfile)
        print(m) if logger is None else logger.error(m)
        exit()
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tEnd read   ({})".format(now)
    print(m) if logger is None else logger.info(m)
    summary_initial_log(obj, logger)

    # Write PDB file from object ==============================================
    filenamepdb = "{}.pdb".format(opts.pattern)
    filenamepsf = "{}.psf".format(opts.pattern)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** WRITE PDB and PSF files ({}) ****** \n".format(now)
    m += "\t\t Writing PDB file...{} ({})".format(filenamepdb, now)
    print(m) if logger is None else logger.info(m)
    obj.write_pdb(filename_pdb=filenamepdb, separate_chains=opts.separatechains)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    if obj._natoms < 99999:
        m = "\t\t CONECT table has been written.\n".format(now)
        isconect = True
    else:
        m = "\t\t CONECT table has NOT been written (natoms > 99999).\n".format(now)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m += "\t\t Writing PSF file...{} ({})\n".format(filenamepsf, now)
        print(m) if logger is None else logger.info(m)
        obj.write_psf(filename_psf=filenamepsf, improper_angles=None)
        isconect = False
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** END WRITE PDB and PSF files ({}) ****** ".format(now)
    print(m) if logger is None else logger.info(m)

    # Renumber topology =======================================================
    headinfo_file = opts.renumberpdb
    if headinfo_file:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        if opts.assignresidues is not None:
            residueinfo_file = opts.assignresidues
        else:
            residueinfo_file = None
        if residueinfo_file is None:
            m = "\n\t\t****** RENUMBERING PDB and PSF files ({}) ****** ".format(now)
            filenamepdb_renumber = "{}_renumber.pdb".format(opts.pattern)
            filenamepsf_renumber = "{}_renumber.psf".format(opts.pattern)
        else:
            m = "\n\t\t****** RENUMBERING AND ASSIGN RESIDUES PDB and PSF files ({}) ****** ".format(now)
            filenamepdb_renumber = "{}_residues.pdb".format(opts.pattern)
            filenamepsf_renumber = "{}_residues.psf".format(opts.pattern)
        m += "\n\t\t Renumbering PDB file... {} ({})".format(filenamepdb, now)
        print(m) if logger is None else logger.info(m)

        if isconect:
            pdbobj = ReadPdbFormat(filenamepdb, isconect=isconect)
        else:
            pdbobj = ReadPdbFormat(filenamepdb, filenamepsf, isconect=isconect)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t Head-Tail info from ... {} ({})".format(opts.renumberpdb, now)
        print(m) if logger is None else logger.info(m)
        head_atoms, tail_atoms = pdbobj.read_head_tail_info(headinfo_file)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t Write renumbered file ... {} ({})".format(filenamepdb_renumber, now)
        print(m) if logger is None else logger.info(m)
        pdbobj.write_renumber_pdb(head_idx_atom=head_atoms, tail_idx_atom=tail_atoms,
                                  assign_residues_info=residueinfo_file, fnameout=filenamepdb_renumber,
                                  isunwrap=opts.isunwrap)
        if opts.isunwrap:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\t\t Coordinates have been unwrapped... ({})".format(now)
            print(m) if logger is None else logger.info(m)

        if not isconect:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\t\t Writing PSF file...{} ({})\n".format(filenamepsf_renumber, now)
            print(m) if logger is None else logger.info(m)
            pdbobj.write_psf(filename_psf=filenamepsf_renumber, improper_angles=None)
        else:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\t\t Writing PSF file...{} ({})\n".format(filenamepsf_renumber, now)
            print(m) if logger is None else logger.info(m)
            pdbobj.write_psf(filename_psf=filenamepsf_renumber, improper_angles=None)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t****** END RENUMBERING PDB and PSF files ({}) ****** ".format(now)
        print(m) if logger is None else logger.info(m)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\t****** SUMMARY ({}) ****** \n".format(now)
        m += "\t\t PDB file: {}\n".format(filenamepdb_renumber)
        m += "\t\t GRO file: {}\n".format(os.path.splitext(filenamepdb_renumber)[0]+".gro")
        if not isconect:
            m += "\t\t PSF file: {}\n".format(filenamepsf_renumber)
        else:
            m += "\t\t PSF file: {}\n".format(filenamepsf_renumber)
        m += "\t\t****** END SUMMARY ({}) ****** ".format(now)
        print(m) if logger is None else logger.info(m)
    else:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\t****** SUMMARY ({}) ****** \n".format(now)
        m += "\t\t PDB file: {}\n".format(filenamepdb)
        m += "\t\t****** END SUMMARY ({}) ****** ".format(now)
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