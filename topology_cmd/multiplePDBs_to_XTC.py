import datetime
import glob
import os.path
import MDAnalysis as mda
from collections import defaultdict
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

    # Print command info
    command_info(opts, logger=logger)

    # PDB files
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\t****** CHECK FOR TER AND CONECT LABELS ({}) ****** \n".format(now)
    pdblist = sorted(glob.glob(opts.pdbdir + "/*.pdb"))
    # Check for TER and CONECT labels in the PDBs and correct the numbers
    for idx, ipdb in enumerate(pdblist):
        with open(ipdb) as f:
            lines = str(f.readlines())
            if lines.count("TER") != 0:
                ister = True
                m += "\t\t\t {0:5d} TER label is found\n".format(idx, now)
            else:
                ister = False
                m += "\t\t\t {0:5d} TER label is not found\n".format(idx, now)
            if lines.count("CONECT") != 0:
                isconect = True
                m += "\t\t\t {0:5d} CONECT label is found\n ".format(idx, now)
            else:
                isconect = False
                m += "\t\t\t {0:5d} CONECT label is not found\n".format(idx, now)
            obj = ReadPdbFormat(ipdb, isconect=isconect, ister=ister, logger=logger)
        # New pdb list, in the case of TER labels new files are saved
        pdblist[idx] = obj._fnamepath

    print(m) if logger is None else logger.info(m)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\t****** END CHECK FOR TER AND CONECT LABELS ({}) ****** \n".format(now)
    print(m) if logger is None else logger.info(m)

    # Renumber all pdb files
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\t****** RENUMBER ALL PDB FILES ({}) ****** \n".format(now)
    print(m) if logger is None else logger.info(m)
    pdblist_new = []
    for idx, ipdb in enumerate(pdblist):
        with open(ipdb) as f:
            lines = str(f.readlines())
            if lines.count("CONECT") != 0:
                isconect = True
            else:
                isconect = False
            obj = ReadPdbFormat(ipdb, isconect=isconect, logger=logger)
            if obj._natoms > 99999:
                isconect = False
                filenamepsf = "{0:s}_{1:04d}.psf".format(opts.pattern, 0)
                obj.write_psf(filename_psf=filenamepsf, improper_angles=None)

        # Renumber topology =======================================================
        headinfo_file = opts.renumberpdb
        if isconect:
            pdbobj = ReadPdbFormat(ipdb, isconect=isconect, logger=logger)
        else:
            pdbobj = ReadPdbFormat(ipdb, filenamepsf, isconect=isconect, logger=logger)
        head_atoms, tail_atoms = pdbobj.read_head_tail_info(headinfo_file)
        ffname_tmp = os.path.splitext(os.path.split(ipdb)[-1])[0]
        filenamepdb_renumber = "{}_renumber.pdb".format(ffname_tmp)
        pdbobj.write_renumber_pdb(head_idx_atom=head_atoms, tail_idx_atom=tail_atoms,
                                  assign_residues_info=None, fnameout=filenamepdb_renumber,
                                  isunwrap=False)
        pdblist_new.append(filenamepdb_renumber)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\t****** END RENUMBER ALL PDB FILES ({}) ****** \n".format(now)
    print(m) if logger is None else logger.info(m)

    # Append all pdb files
    # Output file where concatenated content will be stored
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\t****** APPEND ALL PDB FILES ({}) ****** \n".format(now)
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
                    write_etoe_bb_info(pdbobj)

    print(m) if logger is None else logger.info(m)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\t****** END APPEND ALL PDB FILES ({}) ****** \n".format(now)
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
    m = "\n\t\t****** CLEAR FILES ({}) ****** \n".format(now)
    clear_files = glob.glob("*_mod_*")
    for iclear in clear_files:
        os.remove(iclear)
    clear_files = glob.glob("*_mod.*pdb")
    for iclear in clear_files:
        os.remove(iclear)
    print(m) if logger is None else logger.info(m)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\t****** END CLEAR FILES ({}) ****** \n".format(now)
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