"""
This is a derived class of MolecularGraph. It specializes the molecular graph to a topology.
"""
import numpy as np
import MDAnalysis as mda
import datetime
import os
import utils
from topology.MolecularGraph import MolecularGraph
import topology as top
from collections import defaultdict
import periodictable as PT

"""
Reference 2:    "Automatic Perception of Organic Molecules Based on Essential Structural Information"
                 Yuan Zhao, Tiejun Cheng, and Renxiao Wang*
                 J. Chem. Inf. Model. 2007, 47, 1379-1385

Reference 3:    "A New Algorithm for Exhaustive Ring Perception in a Molecular Graph"
                Th. Hanser, Ph. Jauffret, and G. Kaufmann
                J. Chem. Inf. Comput. Sci. 1996, 36, 1146-1152

ToDo: self._orderbonds can be improved. It is a bad idea to have a natomsxnatoms array

"""


class Topology(MolecularGraph):

    __slots__ = ['_orderbonds', '_nringsCauchy', 'elements', 'charge', 'mass',
                 '_topologyfile', "logger", '_isbackbone', '_totalcharge_mol',
                 '_index_bond_global_dict', '_orderbonds_zero', '_totalcharge',
                 '_types', '_names', '_canonical_order', '_resname']

    # #########################################################################
    def __init__(self, natoms=-1, listbonds=None, undirected=True, logger=None):

        """
        Topology constructor, Topology objects rarely are direcly build by users.
        Segment class initializes the Topology of the system

        Parameters:
            * ``natoms`` (integer , default = -1):  Number of atoms in the topology
            * ``listbonds`` (list, default = None): List of bonds betweem atoms
            * ``undirected`` (boolean, default = True): False for directed graph or True for undirected graph.
            (see warning in `polyanagro.MolecularGraph.MolecularGraph`)
            * ``logger`` (Logger instance): Log the results

        Attributes:
            * ``self._orderbonds`` --> A numpy array (natomsxnatoms) to set up the bond order between two atoms
            * ``self._elements`` --> A list
            * ``self._charge`` --> Atomic charge list
            * ``self._mass`` --> Atomic mass list
            * ``self._nringsCauchy`` -->
            * ``self._topologyfile`` --> Topology file
            * ``self.logger`` (Logger instance): Log the results
            * ``self._isbackbone`` (list of booleans, lenth: natoms): If self._isbackbone[iat]=True, it means
                that iat is a backbone atom


        .. image:: ../../figures/topo01_test06.png

        """

        super().__init__(nvert=natoms, listbonds=listbonds, undirected=undirected)

        if logger is not None:
            self.logger = logger
        else:
            self.logger = None

        self._orderbonds = None
        self._orderbonds_zero = None
        self._index_bond_global_dict = defaultdict(list)
        self._totalcharge_mol = None
        self._totalcharge = None
        self._names = []
        self._resname = []
        self._types = []
        self._nringsCauchy = 0
        self._topologyfile = ""
        self._isbackbone = []
        self._canonical_order = None
        self.elements = []
        self.charge = []
        self.mass = []

    # #########################################################################
    def __copy__(self):

        t = Topology()
        t.natoms = self.natoms
        t._nringsCauchy = self._nringsCauchy
        t._undirected = self._undirected
        t._bonds = self._bonds[:]
        t._cycles = self._cycles[:]
        t._nmols = self._nmols[:]
        t._graphdict = self._graphdict.copy()
        if self._orderbonds is not None:
            t._orderbonds = self._orderbonds.copy()
        else:
            t._orderbonds = None

        t.mass = self.mass[:]
        t.elements = self.elements[:]
        t.charge = self.charge[:]
        t._topologyfile = self._topologyfile
        t._iatch = self._iatch.copy()
        t.logger = self.logger
        t._totalcharge = self._totalcharge
        t._totalcharge_mol = self._totalcharge_mol

        return t

    # #########################################################################
    def __eq__(self, other):

        """
        Overrides equal method

        ``Parameters``:
            * **other** (type: Topology) -->
        """

        if other is None:
            return None

        res = True
        # Get both attributtes from the super and sub clasess because the use of __slots__
        keys = super().__slots__+self.__slots__
        for key in keys:

            #print(key, "self."+key)
            #print(getattr(self,key))

            if isinstance(getattr(self,key), np.ndarray):
                par = np.array_equal(getattr(self,key), getattr(other,key))
                res = res and par
            elif isinstance(getattr(self,key), Topology):
                par = self.__dict__[key] == other.__dict__[key]
                res = res and par
            elif key == "_bonds":
                l1 = getattr(self,key)
                l2 = getattr(other,key)
                for item in l1:
                    par = item in l2
                    res = res and par
            else:
                par = getattr(self,key) == getattr(other,key)
                res = res and par

        return res

    # #########################################################################
    def guess_bonds_topology(self, coords, elements):

        """
        Given a set of coordinate atoms, it guess if a bond exists between two atoms based on the atom distance

        Parameters:
            * ``coords`` (ndarray-float64 (natoms, 3)): Coordinates of the atoms to create bonds
            * ``elements`` (ndarray-string (natoms)): Element for each atom

        Return:
            * ``None``

        """

        natoms = self.natoms

        if np.shape(coords)[0] != natoms:
            raise ValueError('Coord must have same natoms rows. Natoms: {0:d}, Coords: {1:d}'
                             .format(natoms, np.shape(coords)[0]))
        if np.shape(elements)[0] != natoms:
            raise ValueError('Element must have same natoms rows. Natoms: {0:d}, Elements: {1:d}'
                             .format(natoms, np.shape(elements)[0]))

        # Calculate the atom distance matrix
        dist, tmp1, tmp2, tmp3 = top.distance_array(coords, coords)
        # Set up the connectivity of the molecule
        self.detect_connectivity(dist, elements)

        # Elements
        self.elements = elements.tolist()

        for iatom in range(natoms):
            e = self.elements[iatom]
            m = top.atomic_data.atomic_mass[e]
            self.mass.append(m)
            self.charge.append(0.0)

        # Cauchy formula to detect the number of rings in the molecule
        nsegments = len(self.get_forest())
        nbonds = len(self.get_allbonds())
        self._nringsCauchy = nbonds - self.natoms + nsegments
        self.perception_rings()

        self._set_forest()

        self._totalcharge_mol = []
        for imol in range(len(self._nmols)):
            s = 0.0
            for iat in self._nmols[imol]:
                s += self.charge[iat]
            self._totalcharge_mol.append(int(s))

    # #########################################################################
    def guess_nringsCauchy(self):
        """
        Cauchy formula to detect the number of rings in the molecule

        Parameters:
            * ``None``

        Return:
            * ``Number of Cauchy rings``
        """

        nsegments = len(self.get_forest())
        nbonds = len(self.get_allbonds())
        self._nringsCauchy = nbonds - self.natoms + nsegments
        return self._nringsCauchy

    # #########################################################################
    def get_bonds_topologyCONNECTPDB(self, filenamePDB, assign_bo=False, is_unitedatom=False):

        """
        Get the bonds from the PDB CONECT records

        Parameters:
            * ``filenamePDB`` (string): The name of a PDB file containing CONECT records
            * ``assign_bo`` (boolean): If True the bonds are assigned

        Return:
            * ``None``

        Example
        -------
        Topoogy molecule n-hexane

        >>> t = Topology()
        ... fnamePDB = "../data/n-hexane.pdb"
        ... t.get_bonds_topologyCONNECTPDB(filenamePDB=fnamePDB, assign_bo=True)
        ... t.draw_graph_networkx(title="graphs/topo01_test08")

        .. image:: ../../figures/topo01_test08.png

        """

        msg = "\tBuilding bonds from PDB topology...\n"
        if assign_bo:
            msg += "\tSetting bond orders: TRUE"
        else:
            msg += "\tSetting bond orders: FALSE"
        print(msg) if self.logger is None else self.logger.info(msg)
        start_time = datetime.datetime.now()

        # NATOMS
        with open(filenamePDB,'r') as filePDB:

            filePDB.seek(0)
            self._topologyfile = filePDB

            natoms = 0
            while True:
                iline = filePDB.readline()
                if not iline:
                    break
                elif iline.find("ATOM") != -1 or iline.find("HETATM") != -1:
                    natoms += 1
                    name = iline[12:16].replace(" ", "")
                    resname = iline[17:20]
                    e = iline[76:78].replace(" ", "")
                    try:
                        tempfactor = float(iline[60:66])
                    except ValueError:
                        tempfactor = 0.0
                    self.elements.append(e)
                    self._names.append(name)
                    self._resname.append(resname)
                    if tempfactor < 1.0:
                        self._isbackbone.append(0)
                    else:
                        self._isbackbone.append(1)
                    self.charge.append(0.0)
                    m = top.atomic_data.atomic_mass[e]
                    self.mass.append(m)

            for ivert in range(natoms):
                self.add_vertex(ivert)

        # NBONDS
        isthereconnect = False
        with open(filenamePDB, 'r') as filePDB:
            # Go to the begin of the file
            filePDB.seek(0)

            for line in filePDB:
                if not line.startswith('CONECT'):
                    continue
                else:
                    isthereconnect = True
                # The lines containing only CONNECT label are not take into account
                # CONNECT (without numbers)
                if line.split()[1]:
                    iatom = int(line.split()[1])
                    for jatom in line.split()[2:]:
                        self.add_edge([int(iatom)-1, int(jatom)-1])

        self._set_forest()

        self._totalcharge_mol = []
        for imol in range(len(self._nmols)):
            s = 0.0
            for iat in self._nmols[imol]:
                s += self.charge[iat]
            self._totalcharge_mol.append(int(s))

        if assign_bo and not is_unitedatom:
            self.assign_bond_orders()
        self.guess_nringsCauchy()

        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        msg = "\tTIME(Building_bonds_from_PDB_topology): {0:s} seconds".format(str(elapsed_time.total_seconds()))
        print(msg) if self.logger is None else self.logger.info(msg)

    # #########################################################################
    def detect_connectivity(self, distances, elements, test_max_valence=True):

        """
        Identification of bonded atoms using the method proposed in Reference 1
        based on the distances of atoms

        Reference 1:    "A rule-based algorithm for automatic bond type perception"
                Qian Zhang, Wei Zhang, Youyong Li, Junmei Wang, Liling Zhang and Tingjun Hou
                Journal of Cheminformatics 2012, 4:26
                https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-4-26

        Parameters:
            * ``distances`` (ndarray (float64, float64)): Symmetric matrix of distances between atoms
            * ``elements`` (ndarray (str)) : Array of the name of elements
            * ``test_max_valence`` (boolean) : Make the test of maximum valence for each atom in the system

        Return:
            * ``None``


        """
        isbonded = lambda dij,ri, rj : 0.8 < dij < ri+rj+0.4

        for iatom in range(self.natoms):
            for jatom in range(iatom+1, self.natoms):
                d  = distances[iatom, jatom]
                r1 = top.element_cov_radius[elements[iatom]]
                r2 = top.element_cov_radius[elements[jatom]]
                if isbonded(d, r1, r2):
                    self.add_edge([iatom, jatom])
                    if not self._orderbonds is None:
                        self._orderbonds[iatom, jatom] = 0
                        self._orderbonds[jatom, iatom] = 0

        if test_max_valence:
            self.check_atom_max_valence(distances, elements)

    # #########################################################################
    def check_atom_max_valence(self, distances, elements):

        """
        Check number of covalently connected neighbors.
        If the number of neighbors is greater that the value given in maximal valence dictionary
        (`polyanagro.atomic_data.maximal_valences`),
        then remove the edges with the longest distances until match the maximum valence of the atom.

        Parameters:
            * ``distances`` (ndarray (float64, float64)) : Symmetric matrix of distances between atoms
            * ``elements`` (ndarray (str)): Array of the name of elements

        Return:
            * ``None``

        """
        for iatom in range(0, self.natoms):
            e = elements[iatom]
            if e in top.maximal_valences.keys():
                neigh =  self.get_neighbours(iatom)
                n_neigh = len(neigh)

                # For each neighbour
                while n_neigh > top.maximal_valences[e]:
                    neigh =  self.get_neighbours(iatom)
                    max_dist = 0.0
                    iatom_max = -1
                    jatom_max = -1
                    for jatom in neigh:
                        dij = distances[iatom, jatom]
                        if dij > max_dist:
                            max_dist = dij
                            iatom_max = iatom
                            jatom_max = jatom
                    self.remove_edge([iatom_max, jatom_max])
                    n_neigh -= 1

    # #########################################################################
    def get_bonds_topologyPSF(self, filenamePSF, assign_bo=False):

        """
        Get the bonds using a PSF file for the topology

        Parameters:
            * ``filenamePSF`` (string) : The name of a PSF file
                (https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html)
            * ``assign_bo`` (boolean): If True the bonds are assigned

        Return
        ------
            * ``None``

        Example
        -------
        Topoogy 3 chains SCB polymers

        >>> filepsf = "../data/0003Ch-C020-002br04/namd.psf"
        ... t = Topology()
        ... t.get_bonds_topologyPSF(filepsf)

        """

        msg ="\n\tBuilding bonds from PSF topology... \n"
        if assign_bo:
            msg += "\tSetting bond orders: TRUE"
        else:
            msg += "\tSetting bond orders: FALSE"
        print(msg) if self.logger is None else self.logger.info(msg)
        start_time = datetime.datetime.now()

        self._topologyfile = filenamePSF

        # NATOMS
        with open(filenamePSF, 'r') as filePSF:
            # Go to the begin of the file
            filePSF.seek(0)
            while True:
                iline = filePSF.readline()
                if not iline:
                    break
                elif iline.find('!NATOM') != -1:
                    natoms = int(iline.split()[0])
                    for _ in range(natoms):
                        iline = filePSF.readline()
                        name = iline.split()[4]
                        atomtype = iline.split()[5]
                        self._names.append(name)
                        self._types.append(atomtype)
                        self.charge.append(float(iline.split()[6]))
                        mass = float(iline.split()[7])
                        self.mass.append(float("{0:.2f}".format(mass)))
                        # Find atom by molecular mass
                        dd = [k for k, v in top.atomic_mass.items() if v == round(mass,2)]
                        # For united atoms
                        if len(dd) == 0:
                            self.elements.append(atomtype)
                        else:
                            self.elements.append(dd[0])

        for ivert in range(natoms):
            self.add_vertex(ivert)

        # NBONDS
        with open(filenamePSF, 'r') as filePSF:
            # Go to the begin of the file
            filePSF.seek(0)

            while True:
                iline = filePSF.readline()
                if not iline:
                    break
                elif iline.find('!NBOND') != -1:
                    nbonds = int(iline.split()[0])
                    nlines = int(nbonds/4)
                    rest = nbonds % 4
                    if rest > 0: nlines += 1
                    for _ in range(nlines-1):
                        atoms = filePSF.readline().split()
                        for i in range(0,7,2):
                            iat = int(atoms[i]) - 1
                            jat = int(atoms[i+1]) - 1
                            self.add_edge([iat, jat], setforest=False)
                    atoms = filePSF.readline().split()
                    for i in range(0, len(atoms), 2):
                        iat = int(atoms[i]) - 1
                        jat = int(atoms[i+1]) - 1
                        self.add_edge([iat, jat], setforest=False)

        if self._topologyfile.split(".")[-1] == "psf":
            self.assign_backbone_atoms_psf()

        self._set_forest()

        self._totalcharge_mol = []
        for imol in range(len(self._nmols)):
            s = 0.0
            for iat in self._nmols[imol]:
                s += self.charge[iat]
            self._totalcharge_mol.append(int(s))

        if assign_bo:
            self.assign_bond_orders()
        self.guess_nringsCauchy()

        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        msg = "\tTIME(Building_bonds_from_PSF_topology): {0:s} seconds".format(str(elapsed_time.total_seconds()))
        print(msg) if self.logger is None else self.logger.info(msg)

    # #########################################################################
    def get_bonds_topologyMDAnalysis(self, filenametopo, filecoord=None, filetemplate=None,
                                     assign_bo=False, is_unitedatom=False):

        """
        Get the bonds using the python package MDAnalysis (https://www.mdanalysis.org/)

        Parameters:
            * ``filenameTopo`` (string): The name of a topology format available in MDAnalysis
                (https://userguide.mdanalysis.org/1.0.0/formats/format_reference.html)
            * ``filecoord`` (string) : The name of a file containing compatible coordinates in MDAnalysis
            * ``assign_bo`` (boolean): If True the bonds are assigned

        Return:
            * ``None``

        Example
        -------
        Topoogy 3 chains UA SCB polymers

        >>> filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"
        ... top = Topology()
        ... top.get_bonds_topologyMDAnalysis(filename_tpr, assign_bo=False)
        ... top.draw_graph_pygraphviz(title="graphs/topo01_test09")

        .. image:: ../../figures/topo01_test09.png

        """

        msg ="\t*** Building bonds from topology... \n"
        if assign_bo:
            msg += "\tSetting bond orders: TRUE"
        else:
            msg += "\tSetting bond orders: FALSE"
        print(msg) if self.logger is None else self.logger.info(msg)
        start_time = datetime.datetime.now()

        if filecoord is None:
            u = mda.Universe(filenametopo)
        else:
            u = mda.Universe(filenametopo, filecoord)
            if os.path.splitext(filecoord)[-1] == ".dcd":
                m =  "\tFor a DCD from LAMMPS is needed to correct timestep: DONE"
                # DCD uses the AKMA units (https://www.charmm.org/archive/charmm/documentation/basicusage/)
                # In this unit system the unit of time is 20 AKMA time units is .978 picoseconds.
                # It is neccesary to correct the timestep due to the MDAnalysis.Universe asssumes the AKMA system
                dt_new = u.trajectory.dt * 20.45482949774598  # ps
                u = mda.Universe(filenametopo, filecoord, dt = dt_new)

        self._topologyfile = filenametopo

        print("cJJJJ: Topology.py", len(u.atoms))

        # NATOMS
        for iatom in u.atoms:
            self.add_vertex(iatom.index)
            try:
                self.charge.append(iatom.charge)
            except mda.exceptions.NoDataError as e:
                self.charge.append(0.0)
            try:
                element = iatom.element
            except mda.exceptions.NoDataError as e:
                element = iatom.type
            typeatm = iatom.type
            try:
                name = iatom.name
            except mda.exceptions.NoDataError as e:
                name = iatom.type
            try:
                resname = iatom.resname
            except mda.exceptions.NoDataError as e:
                resname = 'UNK'
            mass = float("{0:.2f}".format(iatom.mass))
            self.elements.append(element)
            self.mass.append(mass)
            self._types.append(typeatm)
            self._names.append(name)
            self._resname.append(resname)

        for ibond in u.bonds:
            iat = ibond[0].index
            jat = ibond[1].index
            self.add_edge([iat, jat], setforest=False)

        self._set_forest()

        self._totalcharge_mol = []
        for imol in range(len(self._nmols)):
            s = 0.0
            for iat in self._nmols[imol]:
                s += self.charge[iat]
            self._totalcharge_mol.append(int(s))

        if assign_bo and not is_unitedatom:
            self.assign_bond_orders()
        self.guess_nringsCauchy()
        del u

        if filetemplate is None:
            if self._topologyfile.split(".")[-1] == "psf":
                self.assign_backbone_atoms_psf()
        else:
            self.assign_backbone_atoms_template(filetemplate)

        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        msg = "\tTIME(Building_bonds_from_topology): {0:s} seconds\n".format(str(elapsed_time.total_seconds()))
        msg += "\t*** End Building bonds from topology... \n"
        print(msg) if self.logger is None else self.logger.info(msg)

    # #########################################################################
    def assign_bond_orders(self, select_resonance=0 ):

        """
        This method assigns bond orders to the bonds according to the algorithm reported in
        "Automated simultaneous assignment of bond orders and formal charges"
        Ivan D. Welsh and Jane R. Allison
        J. Cheminform (2019) 11:18

        https://doi.org/10.1186/s13321-019-0340-0

        The function uses the external software **indigo-bondorders** (located in thirdparty/indigo-bondorder).
        This code is compiled and installed in thirdparty/indigox

        The structure to assign bonds needs to have all hydrogen bonds. Thus, united atom models do not work
        with this function.

        Parameters:
            * ``None``

        Return:
            * ``None``

        ..warning:: This only works with all-atom models, for united-atoms the order bond is not assigned correctly.

        """
        try:
            import indigox as ix
        except ModuleNotFoundError:
            msg = "Bond orders cannot be calculated. The indigox module is not installed."
            print(msg) if self.logger is None else self.logger.info(msg)
            return None

        if self._orderbonds is None:
            try:
                self._orderbonds = defaultdict(float)
                self._orderbonds_zero = defaultdict(float)
            except (ValueError, Exception):
                msg = "Bond orders cannot be assigned. There are many atoms in the system."
                print(msg) if self.logger is None else self.logger.info(msg)
                return False

        if self._totalcharge_mol is None:
            msg = "Bond orders cannot be assigned. Total charge is not setup."
            print(msg) if self.logger is None else self.logger.info(msg)
            return False

        # Periodic Table data from indigox
        PTa = ix.PeriodicTable()

        # Convert graph to networkx forest
        nx_list = self.convert_graph_to_networkx_forest()

        # Apply indigo to each molecule
        imol = 0
        bonds_account = 0
        nresonances = 0
        order_bonds_resonance_dict = defaultdict(list)
        #order_bonds = defaultdict(float)
        for inx in nx_list:

            # Build a molecule in the indigox framework
            mol_indigo = ix.Molecule()
            mol_indigo.SetTotalCharge(self._totalcharge_mol[imol])
            # Add all atoms in a dictionary and get the bonds in the
            # framework of indigox program
            all_atoms = dict()

            for iedge in inx.edges():
                i, j = iedge
                if i not in all_atoms:
                    # Element of i
                    e = self.elements[i]
                    all_atoms[i] = mol_indigo.NewAtom(PTa[e])
                    index = all_atoms[i].SetIndex(i)
                    name = e + str(index)
                    all_atoms[i].SetName(name)

                if j not in all_atoms:
                    # Element of j
                    e = self.elements[j]
                    all_atoms[j] = mol_indigo.NewAtom(PTa[e])
                    index = all_atoms[j].SetIndex(j)
                    name = e+str(index)
                    all_atoms[j].SetName(name)

                if all_atoms[i].GetIndex() < all_atoms[j].GetIndex():
                    mol_indigo.NewBond(all_atoms[i], all_atoms[j])
                else:
                    mol_indigo.NewBond(all_atoms[j], all_atoms[i])

            # Setup to use the FPT algorithm with single electrons without preplacing
            # to calculate bond orders and formal charges
            opts = ix.Options.AssignElectrons
            opts.ALGORITHM = opts.Algorithm.FPT
            opts.FPT.ADD_EDGES_TO_TD = False
            opts.FPT.MINIMUM_PROPAGATION_DEPTH = 1
            opts.USE_ELECTRON_PAIRS = False

            # Calculate bond orders and formal charges.
            # Count have the total number of resonance structures
            nresonances = mol_indigo.AssignElectrons()
            nbonds_inx = len(inx.edges())
            tmp_bo_resonance = defaultdict(list)

            if nresonances != 0:
                for iresonance in range(nresonances):
                    tmp_bo_resonance[iresonance] = [0.0] * nbonds_inx
                    mol_indigo.ApplyElectronAssignment(iresonance)
                    index_bond_local = 0
                    for ibond in mol_indigo.GetBonds():
                        i = ibond.GetSourceAtom().GetIndex()
                        j = ibond.GetTargetAtom().GetIndex()
                        if {i, j} in self._bonds:
                            data = sorted([i,j])
                            index_bond_global = int(self._bonds.index({i,j}))
                            #index_bond_local = int(self._bonds.index({i,j})) - bonds_account
                            self._index_bond_global_dict[index_bond_global] = data
                        else:
                            m = "Bond {}-{} does not exist in the topology. Bond orders should be incorrectly assigned\n"
                            print(m) if self.logger is None else self.logger.warning(m)
                            index_bond_global = -1
                            index_bond_local = -1
                        bo = ibond.GetOrder()
                        if bo == bo.SINGLE_BOND:
                            self._orderbonds[index_bond_global] += 1.
                            tmp_bo_resonance[iresonance][index_bond_local] = 1.
                        elif bo == bo.DOUBLE_BOND:
                            self._orderbonds[index_bond_global] += 2.
                            tmp_bo_resonance[iresonance][index_bond_local] = 2.
                        elif bo == bo.TRIPLE_BOND:
                            self._orderbonds[index_bond_global] += 3.
                            tmp_bo_resonance[iresonance][index_bond_local] = 3.
                        else:
                            m = "Bond order cannot be assigned between {} and {} atoms".format(i, j)
                            m + "Bond order: {}".format(bo)
                            print(m) if self._logger is None else self._logger.warning(m)
                        index_bond_local += 1
                    order_bonds_resonance_dict[iresonance] += tmp_bo_resonance[iresonance]

            else:  # Likely is a united atom system (nresonances == 0)

                nresonances = 1
                tmp_bo_resonance[0] = [0.0] * nbonds_inx
                for ibond in mol_indigo.GetBonds():
                    i = ibond.GetSourceAtom().GetIndex()
                    j = ibond.GetTargetAtom().GetIndex()
                    # Remove None label write by indigox
                    i_name = ibond.GetSourceAtom().GetName()[0:-4]
                    j_name = ibond.GetTargetAtom().GetName()[0:-4]
                    if {i, j} in self._bonds:
                        data = sorted([i, j])
                        index_bond_global = int(self._bonds.index({i, j}))
                        index_bond_local = int(self._bonds.index({i, j})) - bonds_account
                        self._index_bond_global_dict[index_bond_global] = data
                    else:
                        m = "Bond {}-{} does not exist in the topology. Bond orders should be incorrectly assigned\n"
                        print(m) if self.logger is None else self.logger.warn(m)
                    # Find equivalences
                    label1 = i_name+"-"+j_name
                    label2 = j_name+"-"+i_name
                    if label1 in top.united_bonds_equivalence.keys():
                        v = top.united_bonds_equivalence[label1]
                        self._orderbonds[index_bond_global] += v
                        tmp_bo_resonance[0][index_bond_local] = v
                    elif label2 in top.united_bonds_equivalence.keys():
                        v = top.united_bonds_equivalence[label2]
                        self._orderbonds[index_bond_global] += v
                        tmp_bo_resonance[0][index_bond_local] = v
                    else:
                        m = "Cannot assing bond order to the united atom system.\n"
                        m += "\tCheck the equivalence tables in atomic_data.py file.\n"
                        m += "\t{} or {} must be defined in that file.\n".format(label1, label2)
                        m += "\ti={}  j={}.\n".format(i, j)
                        m += "\tThis stop the program!!!!\n"
                        print(m) if self.logger is None else self.logger.error(m)
                        exit()

                order_bonds_resonance_dict[0] += tmp_bo_resonance[0]

            istart = bonds_account
            bonds_account += nbonds_inx
            iend = bonds_account
            # This is for systems with different molecules
            for i in range(istart, iend):
                m = self._orderbonds[i] % nresonances
                if m != 0:
                    self._orderbonds[i] = 1.5
                else:
                    self._orderbonds[i] /= nresonances
            imol += 1

        # Correct for aromaticity
        nbonds = len(self._bonds)
        # Return firts resonance structure
        for ibond in range(nbonds):
            self._orderbonds_zero[ibond]  = order_bonds_resonance_dict[0][ibond]

        # Debug =================
        # for iresonance in range(nresonances):
        #     print(order_bonds_resonance_dict[iresonance])
        # print ("cJ", order_bonds)
        # print("============")
        # Debug =================

        return self._orderbonds_zero

    # #########################################################################
    def get_array_mols_neigh(self):

        """
        Return:
            * ``nmols_array``: nmols array with padding (using -1 value). Each list contains
                the vertices include in the molecule
            * ``l_neigh_array``: The element in position 0 contains the neighbors of the node 0 and so on.
                The array is padded with -1.

        Examples
        --------

        >>> t = Topology(natoms=6, listbonds=[(0,1), (1,5), (1,6), (2,3)])

        .. image:: ../../figures/topo02_test06.png

        >>> nmols_array, l_neigh_array = t.get_array_mols_neigh()
        ... nmols_array = [[0, 1, 5], [2, 3, -1], [4, -1, -1]]
        ... l_neigh_array = [[1, -1], [0, 5], [3, -1], [2, -1], [-1, -1], [1, -1]]

       """

        nmols = self._nmols
        # This padding is needed in the case of non-equal molecules.
        nmols_array = utils.padding_list(nmols, fillval=-1)
        nmols_array = np.array(nmols_array, dtype=np.int32)

        # Array of neighbours for each atom
        l_neigh = []
        for item in self._graphdict.values():
            l_neigh.append(item)
        l_neigh_array = utils.padding_list(l_neigh,fillval=-1)
        l_neigh_array = np.array(l_neigh_array, dtype=np.int32)

        return nmols_array, l_neigh_array

    # #########################################################################
    def assign_backbone_atoms_psf(self):

        """
        Assign backbone atoms using the psf file.
        The unused field (field #9) in the atom section.
        (https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html)

        If the value is 0 the atom is backbone otherwise is a branch atom

        ``Parameters``:
            * ``None``

        ``Return``:
            * ``True if the backbone atoms are assigned otherwise False``
        """

        ext_topology = self._topologyfile.split(".")[-1]
        self._isbackbone = self.natoms * [False]
        if ext_topology != "psf":
            m = "\tBackbone atoms cannot be assigned. A psf file is needed\n"
            m+= "\tTopology file: {}".format(self._topologyfile)
            print(m) if self.logger is None else self.logger.info(m)
            return False

        with open(self._topologyfile, 'r') as f:
            lines = f.readlines()
            iat = 0
            first = True
            for iline in lines:
                if first and iline.find("NATOM") == -1:
                    continue
                elif iline.find("NATOM") != -1:
                    nat_read = int(iline.split()[0])

                if iat > nat_read-1:
                    break
                elif first:
                    first = False
                    continue
                else:
                    iat +=1
                    isbb = int(iline.split()[8])
                    if isbb == 0:
                        self._isbackbone[iat-1] = True
                    else:
                        self._isbackbone[iat-1] = False

            m = "\t*** Backbone atoms have been assigned to the topology.\n"
            m+= "\tTopology file: {}".format(self._topologyfile)
            print(m) if self.logger is None else self.logger.info(m)

        return True

    # #########################################################################
    def assign_backbone_atoms_template(self, template_file):

        """
        ToDo --> This function does not work, at this moment
        Args:
            template_file:

        Returns:

        """

        self._isbackbone = self.natoms * [False]
        if not os.path.isfile(template_file):
            m = "\tBackbone atoms cannot be assigned.The template file does not exist.\n"
            m+= "\t\tTemplate file: {}\n".format(template_file)
            m += "\t\tBackbone and brach atoms are not correctly assigned\n"
            print(m) if self.logger is None else self.logger.warning(m)
            return False

    # #########################################################################
    def get_all_bb_bonds(self):

        """
        Get all backbone atoms in a list

        :return:
        """

        ext_topology = self._topologyfile.split(".")[-1]
        if ext_topology != "psf":
            m = "\tBackbone atoms cannot be assigned. A psf file is needed\n"
            m+= "\tTopology file: {}".format(self._topologyfile)
            print(m) if self.logger is None else self.logger.info(m)
            all_bb_bonds = [[0,0]]
            all_bb_bonds = np.array(all_bb_bonds, dtype=np.int32)
            return all_bb_bonds, {}

        all_bonds = self.get_allbonds()
        nbonds_bb_perch = defaultdict(int)
        all_bb_bonds = []
        for item in all_bonds:
            at1 = item[0]
            at2 = item[1]
            ich1 = self._iatch[at1]
            ich2 = self._iatch[at2]
            if not self._isbackbone[at1]: continue
            if not self._isbackbone[at2]: continue
            nbonds_bb_perch[ich1] += 1
            all_bb_bonds.append([at1, at2])

        all_bb_bonds_array = np.asarray(all_bb_bonds, dtype=np.int32)

        return all_bb_bonds_array, nbonds_bb_perch

    # #########################################################################
    def get_nmols(self):

        """
        Get number of molecules in the topology

        Returns:
            A list containing a sublist for each molecule with the atomic indexes for that molecule
        """

        return self._nmols

    # #########################################################################
    def get_topologyfile(self):

        """
        Get the topology file

        Returns:
            Path to the topology file
        """

        return self._topologyfile

    # #########################################################################
    def set_elements(self, elements):

        if not isinstance(elements, list):
            return False

        self.elements = []
        for item in elements:
            self.elements.append(item)

        return True

    # #########################################################################
    def set_charge(self, charge):

        if not isinstance(charge, list):
            return False

        self.charge = []
        for item in charge:
            self.charge.append(item)

        self._totalcharge = sum(self.charge)

        return True

    # #########################################################################
    def set_charge_mol(self, charge):

        if not isinstance(charge, list):
            return False

        self._totalcharge_mol = []

        self.charge = []
        for item in charge:
            self.charge.append(item)

        for imol in self._nmols:
            s = 0
            for idx in imol:
                s += charge[idx]
            self._totalcharge_mol.append(int(round(s)))

        self._totalcharge = 0
        for imol in range(len(self._nmols)):
            self._totalcharge += self._totalcharge_mol[imol]

        return True

    # #########################################################################
    def set_mass(self, mass):

        if not isinstance(mass, list):
            return False

        self.mass = []
        for item in mass:
            self.mass.append(item)

        return True

    # #########################################################################
    def set_resname(self, resname):

        if not isinstance(resname, list):
            return False

        self._resname = []
        for item in resname:
            self._resname.append(item)

        return True

    # #########################################################################
    def set_isbackbone(self, isbackbone):

        if not isinstance(isbackbone, list):
            return False

        self._isbackbone = []
        for item in isbackbone:
            self._isbackbone.append(item)

        return True

    # #########################################################################
    def set_type(self, type_list):

        if not isinstance(type_list, list):
            return False

        self._types = []
        for item in type_list:
            self._types.append(item)

        return True

    # #########################################################################
    def set_name(self, name_list):

        if not isinstance(name_list, list):
            return False

        self._names = []
        for item in name_list:
            self._names.append(item)

        return True

    # #########################################################################
    def canonize_atom_ordering(self):

        """
        Create a canonical order of atoms independent of input following the algorithm in:

        Schneider N. et al. "Get Your Atoms in Order—An Open-Source Implementation of
        a Novel and Robust Molecular Canonicalization Algorithm", J. Chem. Inf. Model. 2015, 55, 10, 2111–2120

        https://pubs.acs.org/doi/full/10.1021/acs.jcim.5b00543
        """

        try:
            import rdkit
        except ModuleNotFoundError:
            msg = "Atom orderding cannot be canonizated. Rdkit library is not available"
            print(msg) if self.logger is None else self.logger.info(msg)
            return None

        rdkit_mol = self.__graph_to_molrdkit()

        #oldorder = tuple(zip(*sorted([(j.GetIdx(), i) for i,j in enumerate(rdkit_mol.GetAtoms())])))[1]
        neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(rdkit.Chem.CanonicalRankAtoms(rdkit_mol))])))[1]

        new_rdkit_mol = rdkit.Chem.RenumberAtoms(rdkit_mol, neworder)

        self._canonical_order = neworder

        return neworder

    # #########################################################################
    def __graph_to_molrdkit(self):

        """
        Get a rdkit mol from the topology graph.

        Returns:
            A rdkit mol instance.
            (https://www.rdkit.org/docs/cppapi/classRDKit_1_1RWMol.html)

        """

        try:
            import rdkit
        except ModuleNotFoundError:
            msg = "Graph to molrdkit conversion cannot be performed. Rdkit library is not available"
            print(msg) if self.logger is None else self.logger.info(msg)
            return None

        mol_rdkit = rdkit.Chem.RWMol()
        nx_graph = self.convert_graph_to_networkx()

        # Add nodes to the rdkit molecule
        nvert = len(nx_graph.nodes())
        for inode in range(nvert):
            nexplicit_H = 0
            for ineigh in self.get_neighbours(inode):
                if self.elements[ineigh] == 'H':
                    nexplicit_H += 1
            element = self.elements[inode]
            atomic_number = getattr(PT, element).number
            a = rdkit.Chem.Atom(atomic_number)
            idx = mol_rdkit.AddAtom(a)

        ibond = 0
        for iedge in nx_graph.edges():
            iat, jat = iedge
            # iiat = node_to_idx[iat]
            # jjat = node_to_idx[jat]
            iiat = iat
            jjat = jat
            if self._orderbonds is not None:
                if self._orderbonds[ibond] == 1.0:
                    mol_rdkit.AddBond(iiat, jjat, order=rdkit.Chem.BondType.SINGLE)
                elif self._orderbonds[ibond] == 1.5:
                    mol_rdkit.AddBond(iiat, jjat, order=rdkit.Chem.BondType.AROMATIC)
                elif self._orderbonds[ibond] == 2.0:
                    mol_rdkit.AddBond(iiat, jjat, order=rdkit.Chem.BondType.DOUBLE)
                elif self._orderbonds[ibond] == 3.0:
                    mol_rdkit.AddBond(iiat, jjat, order=rdkit.Chem.BondType.TRIPLE)
            else:
                mol_rdkit.AddBond(iiat, jjat, order=rdkit.Chem.BondType.SINGLE)
            ibond += 1

        # print("============")
        #print(rdkit.Chem.MolToMolBlock(mol_rdkit))
        # problems = rdkit.Chem.DetectChemistryProblems(mol_rdkit)
        # print(len(problems))
        # for i in range(len(problems)):
        #     print(problems[i].GetType(), problems[i].Message())
        # for iatom in mol_rdkit.GetAtoms():
        #     print(iatom.GetAtomicNum())
        # for ibond in mol_rdkit.GetBonds():
        #     print(ibond.GetEndAtomIdx(), ibond.GetBeginAtomIdx())
        # print("============")

        # Assign bond property --> '_MolFileBondType'
        for ibond in mol_rdkit.GetBonds():
            t = int(ibond.GetBondTypeAsDouble())
            ibond.SetIntProp('_MolFileBondType', t)
            #print(ibond.GetBondType())
            #print(iatom.GetIsAromatic())

        mol_rdkit.ClearComputedProps()
        flag1 = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
        rdkit.Chem.SanitizeMol(mol_rdkit, sanitizeOps=flag1)
        flag2 = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES
        rdkit.Chem.SanitizeMol(mol_rdkit, sanitizeOps=flag2)
        flag3 = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY
        rdkit.Chem.SanitizeMol(mol_rdkit, sanitizeOps=flag3)

        # Ethane can be --> CC or [CH3][CH3] or [H]C([H])([H])C([H])([H])[H]
        # smiles = rdkit.Chem.MolToSmiles(mol_rdkit)
        # print(smiles)

        return mol_rdkit

    # #########################################################################
    def writepsf(self, filename_psf="test.psf", improper_angles=None):

        """
        Write a PSF file using the topology object.
        Improper angles must be given explicitly as a list:
            [[i1, j1, k1, l1], [i2, j2, k2, l2], ...]
        """

        with open(filename_psf, 'w') as f:

            # Title section
            f.write("PSF\n")
            f.write("\n")
            line = str("%8d %s\n" % (2, '!NTITLE'))
            f.write(line)
            # TODO
            line =  "* File created by Topology code *\n"
            line +=  "* by Javier Ramos (IEM-CSIC) *\n"
            now = datetime.datetime.now()
            line += "* DATE: {}\n".format(now)
            f.write(line)
            f.write("\n")

            # Natoms section
            line = str("%8d %s\n" % (self.natoms, '!NATOM'))
            f.write(line)
            for i in range(self.natoms):
                # If self._isbackbone is empty all atoms are labelled as backbones
                try:
                    if self._isbackbone[i]:
                        tmp = '           0'
                    else:
                        tmp = '           1'
                except IndexError:
                    tmp = '           0'
                if len(self._types) == 0:
                    line = str("{0:8d} {1:>4d}  {2:<3d} {3:>3s}  {4:<3s}  {5:<3s}  {6:10.6f}      {7:8.4f}{8:12s}\n"
                           .format(i + 1,
                            self._iatch[i]+1,
                            self._iatch[i]+1,
                            self._resname[i][0:3],
                            self.elements[i],
                            self.elements[i],
                            self.charge[i],
                            self.mass[i],
                            tmp))
                else:
                    line = str("{0:8d} {1:>4d}  {2:<3d} {3:>3s}  {4:<3s}  {5:<3s}  {6:10.6f}      {7:8.4f}{8:12s}\n"
                           .format(i + 1,
                            self._iatch[i]+1,
                            self._iatch[i]+1,
                            self._resname[i][0:3],
                            self.elements[i],
                            self._types[i],
                            self.charge[i],
                            self.mass[i],
                            tmp))
                f.write(line)
            f.write("\n")

            # NBonds section ======================================================
            nbonds = len(self.get_edges())
            line = str("{0:8d} {1:s}\n".format(nbonds, '!NBOND: bonds'))
            f.write(line)
            iline = 1
            for i in range(self.natoms):
                ibond = self.find_all_paths_length(i, 1)
                for j in range(len(ibond)):
                    iat0 = ibond[j][0]
                    iat1 = ibond[j][1]
                    if iat0 < iat1:
                        continue
                    line = str("{0:8d}{1:8d}").format(iat0 + 1, iat1 + 1)
                    if iline % 4 == 0:
                        f.write(line+"\n")
                    else:
                        f.write(line)
                    iline += 1
            if nbonds % 4 != 0:
                f.write("\n\n")
            else:
                f.write("\n")

            # NAngles section =====================================================
            iline = 1
            line = ""
            ntheta = 0
            for i in range(self.natoms):
                itheta = self.find_all_paths_length(i, 2)
                for j in range(len(itheta)):
                    iat0 = itheta[j][0]
                    iat1 = itheta[j][1]
                    iat2 = itheta[j][2]
                    # Exchange extremes
                    if iat0 < iat2:
                        continue
                    ntheta += 1
                    line += str("{0:8d}{1:8d}{2:8d}").format(iat0 + 1, iat1 + 1, iat2 + 1)
                    if iline % 3 == 0:
                        line += "\n"
                    iline += 1

            line0 = str("{0:8d} {1:s}\n".format(ntheta, '!NTHETA: angles'))
            f.write(line0)
            f.write(line)
            if ntheta % 3 != 0:
                f.write("\n\n")
            else:
                f.write("\n")

            # NPhi section ========================================================
            nphi = 0
            iphilist = []
            iphich_tmp = []
            ich_prev = -1
            for i in range(self.natoms):
                iphi = self.find_all_paths_length(i, 3)
                for j in range(len(iphi)):

                    iat0 = iphi[j][0]
                    iat1 = iphi[j][1]
                    iat2 = iphi[j][2]
                    iat3 = iphi[j][3]

                    # Remove repeat dihedrals ===================
                    if ich_prev == -1:
                        ich_prev = self._iatch[iat0]
                        ich_curr = ich_prev
                    else:
                        ich_curr = self._iatch[iat0]

                    if ich_curr != ich_prev:
                        iphich_tmp = []
                        ich_prev = ich_curr
                    iphich_tmp.append([iat0, iat1, iat2, iat3])
                    if iphich_tmp.count([iat3, iat2, iat1, iat0]):
                        continue
                    # Remove repeat dihedrals ====================
                    iphilist.append([iat0, iat1, iat2, iat3])
                    nphi += 1

            iline = 1
            line = ""
            nphi = 0
            for idx in iphilist:
                line += str("{0:8d}{1:8d}{2:8d}{3:8d}").format(idx[3] + 1, idx[2] + 1, idx[1] + 1, idx[0] + 1)
                nphi += 1
                if iline % 2 == 0:
                    line += "\n"
                iline += 1

            line0 = str("{0:8d} {1:s}\n".format(nphi, '!NPHI: dihedrals'))
            f.write(line0)
            f.write(line)
            if nphi % 2 != 0:
                f.write("\n\n")
            else:
                f.write("\n")

            # NImp section ========================================================
            if improper_angles is not None:
                nimpr = len(improper_angles)
                line0 = str("{0:8d} {1:s}\n".format(nimpr, '!NIMPHI: impropers'))
                iline = 1
                line = ""
                for iimpr in improper_angles:
                    line += str("{0:8d}{1:8d}{2:8d}{3:8d}").format(iimpr[0], iimpr[1], iimpr[2], iimpr[3])
                    if iline % 2 == 0:
                        line += "\n"
                    iline += 1
                f.write(line0)
                f.write(line)
                if nimpr % 2 != 0:
                    f.write("\n\n")
                else:
                    f.write("\n")

            else:
                line0 = str("{0:8d} {1:s}\n".format(0, '!NIMPHI: impropers'))
                f.write(line0)
                f.write("\n")

            line0 = str("{0:8d} {1:s}\n".format(0, '!NDON: donors'))
            f.write(line0)
            f.write("\n")
            line0 = str("{0:8d} {1:s}\n".format(0, '!NACC: acceptors'))
            f.write(line0)
            f.write("\n")
            line0 = str("{0:8d} {1:s}\n".format(0, '!NNB'))
            f.write(line0)
            iline = 1
            line1 = ""
