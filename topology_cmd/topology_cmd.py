import datetime
import os.path
import numpy as np
from collections import defaultdict
from utils.Logger import init_logger
from utils.parse_arguments import print_header, parse_arguments, command_info
from topology.readmol.readxsdformat import ReadXsdFormat
from topology.readmol.readpdbformat import ReadPdbFormat
from topology.readmol.readmol2format import ReadMol2Format
from topology.readmol.readdataLAMMPSformat import ReadDataLammpsFormat


# =============================================================================
def summary_initial_log(topobj, dict_messages, log=None):

    # Print info to log
    nows = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    mm1 = "\n\t\tSummarize of the input system ({})\n".format(nows)
    mm = "\t\t" + len(mm1) * "=" + "\n"
    mm += "\t\tInput file                  : {}\n".format(topobj._fnamepath)
    mm += "\t\tCONECT table input file     : {}\n".format(dict_messages['CONECT_in_Table'])
    mm += "\t\tHead info file              : {}\n".format(dict_messages['headinfo_file'])
    mm += "\t\tResidue info file           : {}\n".format(dict_messages['residueinfo_file'])
    mm += "\t\t\n"
    mm += "\t\tNumber of molecules         : {}\n".format(topobj._nmols)
    mm += "\t\tNumber of residues          : {}\n".format(topobj._nres)
    mm += "\t\tNumber of atoms             : {}\n".format(topobj._natoms)
    mm += "\t\tNumber of bonds             : {}\n".format(topobj._nbonds)
    if len(topobj._topology._bonds) > 0:
        mm += "\t\tBond info available         : True\n"
    else:
        mm += "\t\tBond info available         : False\n"
    if topobj._isthere_boxdimension:
        mm += "\t\tSimulation box dimension    : {0:.4f} {1:.4f} {2:.4f} nm\n".format(float(topobj._boxlength[0]/10.),
                                                                                    float(topobj._boxlength[1]/10.),
                                                                                    float(topobj._boxlength[2]/10.))
        mm += "\t\tSimulation box angles       : {0:.2f} {1:.2f} {2:.2f} degrees\n".\
            format(float(topobj._boxangle[0]*180/np.pi), float(topobj._boxangle[1]*180/np.pi),
                   float(topobj._boxangle[2]*180/np.pi))
    else:
        mm += "\t\tSimulation box dimension    : {0:.4f} {1:.4f} {2:.4f} nm\n".format(0.0, 0.0, 0.0)
        mm += "\t\tSimulation box angles       : {0:.2f} {1:.2f} {2:.2f} degrees\n".format(0.0, 0.0, 0.0)
    mm += "\t\t\n"
    mm += "\t\tOutput raw pdb file         : {}\n".format(dict_messages['Output_PDB_file_initial'])
    mm += "\t\tCONECT table output pdb file: {}\n".format(dict_messages['CONECT_out_Table'])
    mm += "\t\t" + len(mm1) * "=" + "\n"
    print(mm1 + mm) if log is None else log.info(mm1 + mm)

    mm = ''
    nw = 1
    if dict_messages['headinfo_file'] is None:
        mm += "\t\t  WARNING {}: Important information about 'Head and Tails' is missing.\n".format(nw)
        mm += "\t\t  WARNING {}: Please review the output files carefully.\n".format(nw)
        nw += 1
    if not dict_messages['CONECT_in_Table']:
        mm += "\n"
        mm += "\t\t  WARNING {}: Topology is absent in the input file.\n".format(nw)
        mm += "\t\t  WARNING {}: The topology has been estimated from atom distances,\n".format(nw)
        mm += "\t\t  WARNING {}: which may introduce errors if overlapping exists in the system.\n".format(nw)
        nw += 1
    mm = mm[:-1]
    print(mm) if log is None else log.info(mm)


# =============================================================================
def main_app(version):

    # Init logger =============================================================
    starttime = datetime.datetime.now()
    dict_messages = defaultdict()
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
        dict_messages['CONECT_in_Table'] = True
    elif ext_inputfile == "PDB":
        m = "\t\tInput is a PDB format.\n"
        m += "\t\tReading... ({})".format(now)
        print(m) if logger is None else logger.info(m)
        # Check if Conect section is present in the PDB
        with open(opts.inputfile) as f:
            lines = str(f.readlines())
            if lines.count("CONECT") != 0:
                isconect = True
                dict_messages['CONECT_in_Table'] = True
            else:
                isconect = False
                dict_messages['CONECT_in_Table'] = False
        obj = ReadPdbFormat(opts.inputfile, isconect=isconect, logger=logger)
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
        dict_messages['CONECT_in_Table'] = True
        obj = ReadDataLammpsFormat(opts.inputfile, opts.filemap)
    else:
        obj = None
        m = "Extension {} cannot be used. " \
            "Allowed formats are XSD, PDB, MOL2 and LAMMPS_data.".format(ext_inputfile)
        print(m) if logger is None else logger.error(m)
        exit()

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tEnd read   ({})\n".format(now)
    print(m) if logger is None else logger.info(m)

    # Write PDB file from object ==============================================
    filenamepdb = "{}.pdb".format(opts.pattern)
    filenamepsf = "{}.psf".format(opts.pattern)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** WRITE PDB and PSF files ({}) ****** \n".format(now)
    m += "\t\t Writing PDB file...{} ({})".format(filenamepdb, now)
    dict_messages['Output_PDB_file_initial'] = filenamepdb
    print(m) if logger is None else logger.info(m)

    obj.write_pdb(filename_pdb=filenamepdb, separate_chains=opts.separatechains)
    obj.check_read_structure()

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    if obj._natoms < 99999:
        m = "\t\t CONECT table has been written.".format(now)
        isconect = True
        dict_messages['CONECT_out_Table'] = True
        print(m) if logger is None else logger.info(m)
    else:
        m = "\t\t CONECT table has NOT been written (natoms > 99999).\n".format(now)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m += "\t\t Writing PSF file...{} ({})\n".format(filenamepsf, now)
        print(m) if logger is None else logger.info(m)

        obj.write_psf(filename_psf=filenamepsf, improper_angles=None)
        isconect = False
        dict_messages['CONECT_Table'] = False
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** END WRITE PDB and PSF files ({}) ****** ".format(now)
    print(m) if logger is None else logger.info(m)

    # Renumber topology =======================================================
    headinfo_file = opts.renumberpdb
    dict_messages['headinfo_file'] = headinfo_file
    dict_messages['residueinfo_file'] = opts.assignresidues
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
        dict_messages['Residue_Files'] = [filenamepdb_renumber, filenamepsf_renumber]
        m += "\n\t\t Renumbering PDB file... {} ({})".format(filenamepdb, now)
        print(m) if logger is None else logger.info(m)

        if isconect:
            pdbobj = ReadPdbFormat(filenamepdb, isconect=isconect, logger=logger)
        else:
            pdbobj = ReadPdbFormat(filenamepdb, filenamepsf, isconect=isconect, logger=logger)

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

        # ToDO: IMPROPER GUESSING
        # pdbobj.guessing_impropers()

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

    summary_initial_log(obj, dict_messages, logger)

    # Final message ===========================================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    endtime = datetime.datetime.now()
    m = "\n\t\tJob  Done at {} ============".format(now)
    print(m) if logger is None else logger.info(m)
    m = "\t\tTotal time: {0:.2f} seconds".format((endtime-starttime).total_seconds())
    print(m) if logger is None else logger.info(m)


# =============================================================================
if __name__ == "__main__":

    __version__ = "1.1"
    main_app(__version__)
