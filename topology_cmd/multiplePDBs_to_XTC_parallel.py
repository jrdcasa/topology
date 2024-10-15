import datetime
import glob
import os.path
import MDAnalysis as mda
from collections import defaultdict
import multiprocessing as mp


from utils.Logger import init_logger
from utils.parse_arguments import print_header_multiplepdb, parse_arguments_multiplepdb, command_info
from topology.readmol.readpdbformat import ReadPdbFormat


# ==================================================
def write_etoe_bb_info(trjobj):

    fee = open("listendtoend_replicate.dat", 'w')
    fbb = open("backbone_idx_replicate.dat", 'w')
    faa = open("allatom_idx_replicate.dat", 'w')
    fee.writelines("# ich head tail\n")

    ich = 0
    itail = -1
    ihead = -1
    for item in trjobj._topology._nmols:
        fbb.writelines("[mol{}]\n".format(ich))
        faa.writelines("[mol{}]\n".format(ich))
        aa_ch = []
        for idx in item:
            if trjobj._atom3d_occupancy[idx] == 1:
                ihead = idx
            elif trjobj._atom3d_occupancy[idx] == 2:
                itail = idx
            aa_ch.append(idx)

        # Get the backbone atoms
        queue = [ihead]
        bb_ch = [ihead]
        visit_bb = [ihead]
        while queue:
            ibb = queue.pop()
            neighs = trjobj._topology._graphdict[ibb]
            for jbb in neighs:
                if trjobj._atom3d_bfactor[jbb] == 0 and jbb not in visit_bb:
                    bb_ch.append(jbb)
                    visit_bb.append(jbb)
                    queue.append(jbb)

        for ibb in bb_ch:
            fbb.writelines("{}\n".format(ibb))
        for iaa in aa_ch:
            faa.writelines("{}\n".format(iaa))
        try:
            fee.writelines('{} {} {}\n'.format(ich, ihead, itail))
        except Exception as e:
            print(e)
            continue
        ich += 1

    fee.close()
    fbb.close()


# ============================================================================
# Function to check each PDB file for TER and CONECT labels
def check_pdb_file(idx_ipdb):

    idx, ipdb = idx_ipdb
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = ""

    with open(ipdb) as f:
        lines = str(f.readlines())
        if lines.count("TER") != 0:
            ister = True
            m += "\t\t\t {0:5d} TER label is found\n".format(idx)
        else:
            ister = False
            m += "\t\t\t {0:5d} TER label is not found\n".format(idx)
        if lines.count("CONECT") != 0:
            isconect = True
            m += "\t\t\t {0:5d} CONECT label is found\n ".format(idx)
        else:
            isconect = False
            m += "\t\t\t {0:5d} CONECT label is not found\n".format(idx)

    # Create ReadPdbFormat object (assuming this is already defined)
    obj = ReadPdbFormat(ipdb, isconect=isconect, ister=ister, logger=None)

    # Return the updated PDB path and the log message for this file
    return obj._fnamepath, m


# Function to process a single PDB file (renumber and check)
# =========================================================================
def renumber_pdb_file(args):
    idx, ipdb, opts = args
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\t****** RENUMBERING PDB FILE ({}) ****** \n".format(now)

    with open(ipdb) as f:
        lines = str(f.readlines())
        isconect = lines.count("CONECT") != 0

    obj = ReadPdbFormat(ipdb, isconect=isconect, logger=None)

    # Check for large atom counts
    if obj._natoms > 99999:
        isconect = False
        filenamepsf = "{0:s}_{1:04d}.psf".format(opts.pattern, 0)
        obj.write_psf(filename_psf=filenamepsf, improper_angles=None)

    # Renumber topology =======================================================
    headinfo_file = opts.renumberpdb

    if opts.assignresidues is not None:
        residueinfo_file = opts.assignresidues
    else:
        residueinfo_file = None

    if isconect:
        pdbobj = ReadPdbFormat(ipdb, isconect=isconect, logger=None)
    else:
        pdbobj = ReadPdbFormat(ipdb, filenamepsf, isconect=isconect, logger=None)

    head_atoms, tail_atoms = pdbobj.read_head_tail_info(headinfo_file)
    ffname_tmp = os.path.splitext(os.path.split(ipdb)[-1])[0]
    filenamepdb_renumber = "{}_renumber.pdb".format(ffname_tmp)

    # Write the renumbered PDB file
    pdbobj.write_renumber_pdb(head_idx_atom=head_atoms, tail_idx_atom=tail_atoms,
                              assign_residues_info=residueinfo_file, fnameout=filenamepdb_renumber,
                              isunwrap=False)

        # Write end-to-end and backbone information
    if idx == 0:
        write_etoe_bb_info(pdbobj)

    # Return the new PDB path and the log message
    return filenamepdb_renumber, m


# =============================================================================
def main_app(version=None):

    # Init logger =============================================================
    logger = init_logger("Output", fileoutput="InfoPDBtoXTC.log",
                         append=False, inscreen=True)
    print_header_multiplepdb(version, logger)
    opts = parse_arguments_multiplepdb()
    starttime = datetime.datetime.now()
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    logger.info("\t\tStarting: \t {}\n".format(now))
    dict_messages = defaultdict()
    if opts.ncpus is None:
        ncpus = mp.cpu_count()
    else:
        ncpus = opts.ncpus

    m = "\t\tUsing {} processors \n".format(ncpus)
    print(m) if logger is None else logger.info(m)

    # Print command info
    command_info(opts, logger=logger)

    # PDB files Check for TER and CONECT LABELS =========================================================
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\t****** CHECK FOR TER AND CONECT LABELS ({}) ****** \n".format(now)
    print(m) if logger is None else logger.info(m)
    # Get the list of PDB files
    pdblist = sorted(glob.glob(opts.pdbdir + "/*.pdb"))

    # Prepare arguments as list of tuples (index, pdbfile)
    idx_ipdb_list = list(enumerate(pdblist))

    # Create a pool of workers to process the PDB files in parallel
    with mp.Pool(processes=ncpus) as pool:
        results = pool.map(check_pdb_file, idx_ipdb_list)

    print(m) if logger is None else logger.info(m)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** END CHECK FOR TER AND CONECT LABELS ({}) ****** \n".format(now)
    print(m) if logger is None else logger.info(m)
    #
    # # Renumber all pdb files =========================================================
    # Get the current time for the log header
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** RENUMBER ALL PDB FILES ({}) ****** \n".format(now)
    # Log the message (either print or use logger)
    print(m) if logger is None else logger.info(m)

    # Prepare arguments as a list of tuples (index, pdbfile, opts)
    pdblist = []
    for i in results:
        pdblist.append(i[0])
    idx_ipdb_list = [(idx, ipdb, opts) for idx, ipdb in enumerate(pdblist)]

    # Use multiprocessing to process the files in parallel
    pdblist_new = []
    with mp.Pool(processes=ncpus) as pool:
        results = pool.map(renumber_pdb_file, idx_ipdb_list)

    # Collect the results
    for idx, (renumbered_pdb_path, log_message) in enumerate(results):
        pdblist_new.append(renumbered_pdb_path)  # Store new renumbered PDB path
        m += log_message  # Collect log messages

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** END RENUMBER ALL PDB FILES ({}) ****** \n".format(now)
    print(m) if logger is None else logger.info(m)

    # Append all pdb files
    # Output file where concatenated content will be stored
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** APPEND ALL PDB FILES ({}) ****** \n".format(now)
    output_file_trj = '{}_trj.pdb'.format(opts.pattern)
    output_file_top = '{}_topo.pdb'.format(opts.pattern)

    # Open the output file in write mode
    with open(output_file_trj, 'w') as ftrj:
        # Loop through each PDB file
        for idx, ipdb in enumerate(pdblist_new):
            # Open each file in read mode
            with open(ipdb, 'r') as infile:
                # Write the content of each file to the output file
                ftrj.write("MODEL{0:>9d}\n".format(idx+1))
                lines = infile.readlines()
                for iline in lines:
                    if iline.count("ATOM") != 0 or iline.count("HETATM") != 0 :
                        ftrj.write(iline)
                    elif iline.count("CRYST1") != 0 or \
                            iline.count("ORIGX") != 0 or\
                            iline.count("SCALE"):
                        ftrj.write(iline)
                # Optionally, add a newline between concatenated PDB files
                ftrj.write('ENDMDL\n')
                if idx == 0:
                    with open(output_file_top, 'w') as ftop:
                        for iline in lines:
                            ftop.write(iline)
    print(m) if logger is None else logger.info(m)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** END APPEND ALL PDB FILES ({}) ****** \n".format(now)
    print(m) if logger is None else logger.info(m)

    # Transform PDB trajectory to XTC one
    updb = mda.Universe(output_file_trj)
    output_file_xtc = os.path.splitext(output_file_trj)[0]+".xtc"

    # Create a writer for the XTC file
    with mda.Writer(output_file_xtc, n_atoms=updb.atoms.n_atoms) as writer:
        # Loop through the trajectory and write each frame to the XTC file
        for ts in updb.trajectory:
            writer.write(updb.atoms)

    # Clear files:
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** CLEAR FILES ({}) ****** \n".format(now)
    clear_files = glob.glob("*_mod_*")
    for iclear in clear_files:
        os.remove(iclear)
    clear_files = glob.glob("*_mod.*pdb")
    for iclear in clear_files:
        os.remove(iclear)
    print(m) if logger is None else logger.info(m)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t****** END CLEAR FILES ({}) ****** \n".format(now)
    print(m) if logger is None else logger.info(m)

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