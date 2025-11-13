import datetime
import MDAnalysis as Mda
from topology.GromacsFormat import GromacsFormat
from topology.RaspaFormat import RaspaFormat
from utils.Logger import init_logger
from utils.parse_arguments import parse_arguments_type_gas, \
    print_header_type_gas, command_info


# =============================================================================
def main_app(version=None):

    """

    Args:
        version:

    Returns:

    """

    # Init logger =============================================================
    opts = parse_arguments_type_gas()
    logger = init_logger("Output", fileoutput=opts.log,
                         append=False, inscreen=True)
    print_header_type_gas(version, logger)
    starttime = datetime.datetime.now()
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    logger.info("\t\tStarting: \t {}\n".format(now))

    # Print command info
    command_info(opts, logger=logger)

    # Create an universe
    u = Mda.Universe(opts.pdbfile, opts.pdbfile)

    # Gromacs files
    m = "\t\t =========== WRITE TOPOLOGY FILES FOR GROMACS ===========\n"
    print(m) if logger is None else logger.info(m)
    g = GromacsFormat(u, logger=logger)
    g.write_gro()
    g.write_mdp()
    g.write_top(opts.ffname, opts.chargefile, opts.impropers)

    # RASPA files
    m = "\t\t =========== WRITE TOPOLOGY FILES FOR RASPA ===========\n"
    print(m) if logger is None else logger.info(m)
    r = RaspaFormat(u, logger=logger)
    r.write_def_molecules(opts.ffname, opts.chargefile, opts.impropers)
    r.write_json_files(opts.ffname, opts.chargefile, opts.impropers)
    r.write_pseudo_atoms(opts.ffname)
    r.write_ff_vdw(opts.ffname)
    m  = "\t\t                   ATTENTION                          \n"
    m += "\t\t Please note that the files generated for the RASPA program\n"
    m += "\t\t are intended to be used as templates.\n"
    m += "\t\t User intervention is always required to edit them."
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
