"""
Class to handle trajectories. This make a heavy use of MDAnalysis.
In fact, it works as a wrapper to the MDAnalysis library
"""
import MDAnalysis
import datetime
import os
import numpy as np
import topology as top
from copy import copy


class ExtTrajectory(object):

    """
    Extended Trajectory class.

    A list of trajectory paths are given in order to create a Extended Trajectory object

    """

    __slots__ = ['trjtype', 'trjpath', 'dt', 'topology',
                 'logger', 'universe', 'nframes', 'natoms', 'trjlist']

    # ***********************************************************************************
    def __init__(self, trjlist=None, topfile=None, templatefile=None, logger=None, is_unitedatom=False):

        """
        Trajectory Constructor.

        Parameters:
            * ``trjlist`` (list of string): Path list of the trajectories
            * ``topfile`` string : Path to the topology file. If None the topology is empty
            * ``logger`` (Logger instance): Log the results. If logger=None, the output is on the screen

        Attributes:
            * ``self.typetrj`` (list) --> Type of trajectories, example: ['xtc','xtc','xtc']
            * ``self.trjpath`` (list) --> List of the trajectory paths
            * ``self.trjlist`` (list) --> List if the Trajectory objects (MDAnalysis)
            * ``self.dt`` (float) --> Time between frames in ps
            * ``self.topology`` (Topology) --> Topology of the trajectories
            * ``self.logger`` (Logger instance) --> Log the results
            * ``self.universe`` () --> MDAnalysis Universe
            * ``self.nframes (int) --> Number of frames of the trajectory

        """

        start_time = datetime.datetime.now()

        if not isinstance(trjlist, list):
            trjlist = [trjlist]

        if logger is not None:
            self.logger = logger
            msg = "\t*** Creating a trajectory from..."
            print(msg) if self.logger is None else self.logger.info(msg)
        else:
            self.logger = None

        self.universe = None
        self.trjtype = []
        self.trjpath = trjlist
        self.natoms = 0
        self.nframes = 0
        self.topology = None
        self.trjlist = None

        msg = ""
        for item in trjlist:
            msg += "\t"+item+"\n"
        print(msg) if self.logger is None else self.logger.info(msg)

        try:
            for item in self.trjpath:
                ext = item.split(".")[-1]
                self.trjtype.append(ext)
        except TypeError:
            pass

        if not topfile:
            t = top.Topology(logger=self.logger)
            self.set_top_universe(t)
        else:
            t = top.Topology(logger=self.logger)
            t.get_bonds_topologyMDAnalysis(topfile, filecoord=trjlist[0],
                                           filetemplate=templatefile, is_unitedatom=is_unitedatom)
            self.set_top_universe(t)

        try:
            self.dt = self.universe.trajectory.dt
        except AttributeError:
            self.dt = 0.0

        if self.universe is not None:
            msg = "\tNumber of frames: {}".format(self.nframes)
            print(msg) if self.logger is None else self.logger.info(msg)

        self.natoms = len(self.universe.atoms)

        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        msg = "\tTIME(Read_Trajectories): {0:s} seconds\n".format(str(elapsed_time.total_seconds()))
        msg += "\t*** End Creating a trajectory \n"
        # print(msg) if self.logger is None else self.logger.info(msg)

    # ***********************************************************************************
    def __copy__(self):

        trj = ExtTrajectory(logger=self.logger)
        trj.trjtype = self.trjtype[:]
        trj.trjpath = self.trjpath[:]
        trj.topology = copy(self.topology)
        trj.dt = self.dt
        trj.logger = self.logger
        trj.universe = copy(self.universe)

        return trj

    # ***********************************************************************************
    def add_trjs(self, fullnametrj, append=True):

        """
        Add trajectories from a list (or set) and check that its type is equal to the others added trajectories.
        The trajectory names have to be an order. It is supposed that the first trajectory is the first name
        in an alphabetically order and so on.

        Modify self.typetrj, self.trjpath

        Parameters:
            * ``fullnametrjset`` (string): Path of the new trajectory file
            * ``append`` (boolean, default=True) : Append (True) or First (False)

        Returns:
            * ``Ǹone``

        """

        if append:
            self.trjpath.append(fullnametrj)
            self.trjtype.append(fullnametrj.split(".")[-1])
        else:
            self.trjpath.insert(0, fullnametrj)
            self.trjtype.insert(0, fullnametrj.split(".")[-1])

        try:
            self.universe.load_new(self.trjpath, continuous=True)
            self.nframes = self.get_numframes()
            msg = "\tTotal number of frames: {0}".format(self.nframes)
            print(msg) if self.logger is None else self.logger.info(msg)
            msg = "\tAdding the trajectory {0:s}".format(fullnametrj)
            print(msg) if self.logger is None else self.logger.info(msg)
        except (FileNotFoundError, TypeError) as e:
            msg = "\n\t{0:s} is not found or is not an allowed file\n".format(fullnametrj)
            msg += "\tTrajectory object is empty.\n"
            msg += "{}".format(e)
            print(msg) if self.logger is None else self.logger.error(msg)
            pass

        return None

    # ***********************************************************************************
    def get_numframes(self):

        """
        Get the total number of frames

        Parameters:
            * ``None``

        Returns:
            * ``Ǹone``

        """
        self.nframes = len(self.universe.trajectory)
        return self.nframes

    # ***********************************************************************************
    def num_atoms(self):

        """
        Get the total number of frames

        Parameters:
            * ``None``

        Returns:
            * ``None``

        """
        self.natoms = len(self.universe.atoms)
        return self.natoms

    # ***********************************************************************************
    def trajectory_reader_mdanalysis(self, trjname):

        """
        Reader for XTC, DCD and TRR trajectories.

        Set up the attribute self.trjlist

        Parameters:
            * ``trjpath``: Name of the trajectory to load
            * ``removefirst``: Default: True --> Remove the first frame of the trajectory

        Return:
            * ``None``
        """

        try:
            # Try to open the file to see if it exists or not.
            open(trjname, 'r').close()
            if self.trjtype[-1] == "xtc":
                self.trjlist.append(MDAnalysis.coordinates.XTC.XTCReader(trjname))
            elif self.trjtype[-1] == "trr":
                self.trjlist.append(MDAnalysis.coordinates.TRR.TRRReader(trjname))
            elif self.trjtype[-1] == "dcd":
                self.trjlist.append(MDAnalysis.coordinates.DCD.DCDReader(trjname))
        except (FileNotFoundError, OSError):
            msg = "\n\t{0:s} format file is wrong.\n".format(self.trjpath)
            msg += "\tTrajectory object is empty.\n"
            print(msg) if self.logger is None else self.logger.error(msg)

    # ***********************************************************************************
    def set_top_universe(self, topology):

        """
        Add topology to the Trajectory object.

        This method changes self.topology attributr

        Parameters:
            * ``top`` (Topology) : A topology instance

        Return:
            * ``None``

        """
        start_time = datetime.datetime.now()

        msg = "\t*** Adding topology \n\t{0:s}".format(topology.get_topologyfile())
        print(msg) if self.logger is None else self.logger.info(msg)

        self.topology = topology

        if self.universe is None:
            if len(self.topology.get_topologyfile()) != 0:
                ext = os.path.splitext(self.trjpath[0])[-1]
                # Only XTC and TRR formats can be treated as virtual trajectories in MDAnalysis (continuous=True)
                if ext.upper() == ".TRR" or ext.upper() == ".XTC":
                    self.universe = MDAnalysis.Universe(self.topology.get_topologyfile(),
                                                        self.trjpath, continuous=True)
                else:
                    # Be careful if the last frame of the trajectory file is the same that the
                    # first frame of the next trajectory, the frame will be repeat in the universe.
                    self.universe = MDAnalysis.Universe(self.topology.get_topologyfile(),
                                                        self.trjpath, continuous=False)
                    if ext.upper() == ".DCD":
                        m = "\tFor a DCD from LAMMPS is needed to correct timestep: DONE"
                        # DCD uses the AKMA units (https://www.charmm.org/archive/charmm/documentation/basicusage/)
                        # In this unit system the unit of time is 20 AKMA time units is .978 picoseconds.
                        # It is neccesary to correct the timestep due to the MDAnalysis.Universe asssumes the AKMA system
                        dt_new = self.universe.trajectory.dt * 20.45482949774598  # ps
                        self.universe = MDAnalysis.Universe(self.topology.get_topologyfile(),
                                                            self.trjpath, continuous=False, dt=dt_new)
                        print(msg) if self.logger is None else self.logger.info(msg)
                self.nframes = len(self.universe.trajectory)
                self.dt = self.universe.trajectory.dt

                msg = "\tNumber of frames: {}".format(self.nframes)
                print(msg) if self.logger is None else self.logger.info(msg)

                end_time = datetime.datetime.now()
                elapsed_time = end_time - start_time
                msg = "\tTIME(Adding topology): {0:.2f} seconds".format(elapsed_time.total_seconds())
                print(msg) if self.logger is None else self.logger.info(msg)
            else:
                msg = "\tEmpty topology in trajectory"
                print(msg) if self.logger is None else self.logger.info(msg)

        msg = "\t*** End Adding topology\n"
        print(msg) if self.logger is None else self.logger.info(msg)

    # ***********************************************************************************
    def write_trajectory(self, filename, pbc=True, nojump=False, format_trj=None):

        """
        Write trajectory

        Parameters:
            * ``filename`` (string) : Name of the trayectory file
            * ``pbc`` (boolean, default True) : if True remove PBC. Ignored if ``nojump`` is True
            * ``nojump`` (boolean, default False) : if True remove jumps between frames
            * ``format`` (string): Format of the trajectory. The format must be supported by MDAnalysis package
            (https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html)

        Return:
            * ``None``

        """

        msg = "\n\tWriting trajectory... {0:s}".format(filename)
        print(msg) if self.logger is None else self.logger.info(msg)

        natoms = self.topology.natoms
        n_image = np.zeros((natoms, 3), dtype=np.int32)
        nmols_array, l_neigh_array = self.topology.get_array_mols_neigh()

        # Check format
        if format_trj is None:
            W = MDAnalysis.Writer(filename, n_atoms=natoms)
        else:
            _filename = filename + "." + format_trj
            W = MDAnalysis.Writer(_filename, n_atoms=natoms, format=format_trj)

        start_time = datetime.datetime.now()
        indx_iframe = 0
        c_unwrap_0 = None

        # Universe with a virtual contiguous trajectory
        # ########u = MDAnalysis.Universe(self.topology._topologyfile, self.trjpath, continuous=True)
        s_atoms = self.universe.select_atoms("all")

        for ts in self.universe.trajectory:

            # Using unwrapping function from MDAnalysis library
            # if pbc:
            #     s_atoms.unwrap()

            if indx_iframe % 100 == 0:
                msg = "\tFrame {0:9d} of {1:9d}".format(indx_iframe, self.get_numframes())
                mid_time = datetime.datetime.now()
                elapsed_time = mid_time - start_time
                msg += "\ttime: {0:.2f} seconds".format(elapsed_time.total_seconds())
                print(msg) if self.logger is None else self.logger.info(msg)

            # Write coordinates without jump in the molecules (chains)
            if nojump:
                # 1. Get box dimensions
                box = ts.dimensions
                # 2. Unwrap coordinates
                if indx_iframe == 0:
                    c_unwrap_0 = top.unwrap(ts.positions, nmols_array, l_neigh_array, box)
                    # 4. Update coordinates in the trajectory file
                    # for iatom in range(len(c_unwrap_0)):
                    #    ts.positions[iatom] = c_unwrap_0[iatom,:]
                    ts.positions = c_unwrap_0
                else:
                    c_jump = top.unwrap_nojump(ts.positions, c_unwrap_0, n_image, nmols_array, l_neigh_array, box)
                    c_unwrap_0 = top.unwrap(ts.positions, nmols_array, l_neigh_array, box)
                    # 4. Update coordinates in the trajectory file
                    # for iatom in range(len(c_jump)):
                    ts.positions = c_jump
                # 5. Write coordinates
                W.write(s_atoms)
                indx_iframe += 1
            # Remove PBC
            elif pbc and not nojump:
                # 1. Get box dimensions
                box = ts.dimensions
                # 2. Unwrap coordinates
                # c_unwrap = unwrap(ts.positions, nmols_array, l_neigh_array, box)
                # # 3. Update coordinates in the trajectory file
                # for item in range(len(c_unwrap)):
                #     ts.positions[item] = c_unwrap[item]
                ts.positions = top.unwrap(ts.positions, nmols_array, l_neigh_array, box)
                indx_iframe += 1
                # 5. Write coordinates
                W.write(s_atoms)
            # Write trj as is.
            else:
                indx_iframe += 1
                W.write(s_atoms)
