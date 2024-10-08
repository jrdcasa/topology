from MDAnalysis import Universe, exceptions
from MDAnalysis.topology.core import guess_atom_element
import topology as top
import datetime
from topology.readmol.readbaseformat import ReadBaseFormat
from topology.internal_coordinates import unwrap
from collections import defaultdict
from topology.atomic_data import atomic_number
import copy
import os
import numpy as np
import warnings

warnings.simplefilter("ignore", Warning)


class ReadPdbFormat(ReadBaseFormat):

    # *************************************************************************
    def __init__(self, filenamepath, filenametopo=None, inputbondlist=None,
                 assign_bondorders=False, isconect=True, ister=False, logger=None):

        super().__init__(filenamepath=filenamepath, assign_bondorders=assign_bondorders, logger=logger)

        self._atom3d_bfactor = defaultdict()
        self._atom3d_occupancy = defaultdict()
        self._atom3d_kindmolecule = defaultdict()
        self._dimensions = list()
        self._bond_list = list()

        # TER labels can be problematic, remove them and modify atoms and connect matrix
        if ister:
            filenamepath = self._removeterm_and_savepdb(filenamepath)
            self._fnamepath = filenamepath

        if filenamepath is not None:
            if filenametopo is not None:
                self._fnametopopath = filenametopo
                self.read_pdb(guess_bonds=False)
            else:
                self._fnametopopath = None
                if isconect:
                    self.read_pdb(guess_bonds=False)
                else:
                    self.read_pdb(guess_bonds=True)

    # *************************************************************************
    def __copy__(self):

        t = ReadPdbFormat(filenamepath=self._fnamepath, assign_bondorders=self._assign_bo)

        return t

    # *************************************************************************
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
        keys = super().__slots__
        for key in keys:
            print(key)
            if isinstance(getattr(self, key), np.ndarray):
                par = np.array_equal(getattr(self, key), getattr(other, key))
                res = res and par
            elif isinstance(getattr(self, key), top.Topology):
                par = getattr(self, key) == getattr(other, key)
                res = res and par
            elif isinstance(getattr(self, key), Universe):
                # Assuming that both universe are identical.
                # In my best knownledge, Universe in MDAnalysis cannot be compared
                res = res and True
            else:
                par = getattr(self, key) == getattr(other, key)
                res = res and par
            if not res:
                break
        return res

    # *************************************************************************
    def read_pdb(self, guess_bonds):

        """
        Base function to read the PDB format
        Returns:

        """

        degtorad = np.pi / 180.0

        # MDAnalysis Universe
        if guess_bonds:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\t\t\t Start Guessing bonds using MDAnalysis library. ({})".format(now)
            print(m) if self._logger is None else self._logger.info(m)
        self._universe = Universe(self._fnamepath, guess_bonds=guess_bonds, fudge_factor=0.6)
        if guess_bonds:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            m = "\t\t\t End Guessing bonds using MDAnalysis library. ({})".format(now)
            print(m) if self._logger is None else self._logger.info(m)

        self._natoms = self._universe.atoms.n_atoms
        self._dimensions = self._universe.dimensions
        # Check for bond information
        try:
            self._nbonds = len(self._universe.bonds)
        except exceptions.NoDataError:
            if self._fnametopopath is not None:
                u_tmp = Universe(self._fnametopopath)
                self._nbonds = len(u_tmp.bonds)
                bondlist = []
                for ibond in u_tmp.bonds:
                    i = ibond.indices[0]
                    j = ibond.indices[1]
                    bondlist.append([i, j])
                self._universe.add_TopologyAttr('bonds', bondlist)
                del u_tmp
            else:
                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                m = "\t\t\t Start Guessing bonds using MDAnalysis library. ({})".format(now)
                print(m) if self._logger is None else self._logger.info(m)

                self._universe = Universe(self._fnamepath, guess_bonds=True)

                now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
                m = "\t\t\t End Guessing bonds using MDAnalysis library. ({})".format(now)
                print(m) if self._logger is None else self._logger.info(m)

                self._nbonds = len(self._universe.bonds)

        self._nres = len(self._universe.residues)

        # Loop over atoms
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\t\t Start Loop over atoms. ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)
        listelements = list()
        for iatom in self._universe.atoms:
            try:
                e = iatom.element
            except exceptions.NoDataError:
                # Avoid mistakes when the atom name matches some element name
                # Example: CS is not Cs

                if len(iatom.name) > 1 and not iatom.name[1].islower():
                    #    iatom.name = iatom.name[0]
                    element = iatom.name[0]
                e = str(guess_atom_element(element))
                listelements.append(e)
            # cJ idx = iatom.id - 1
            idx = iatom.index
            self._atom3d_element[idx] = e
            self._atom3d_xyz[idx] = list(iatom.position)
            self._atom3d_molname[idx] = iatom.resname
            self._atom3d_isbackbone[idx] = iatom.tempfactor
            self._atom3d_bfactor[idx] = iatom.tempfactor
            self._atom3d_occupancy[idx] = iatom.occupancy
            try:
                self._atom3d_kindmolecule[idx] = iatom.chainID
            except KeyError:
                self._atom3d_kindmolecule[idx] = "X"
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t\t End Loop over atoms. ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)

        # Check if there is information about head-tail atoms
        # for idx, value in self._atom3d_occupancy.items():
        #     if value == 1.0:
        #         self._heads.append(value)
        #     elif value == 2.0:
        #         self._tails.append(value)
        # nheads = len(self._heads)
        # ntails = len(self._tails)

        if len(listelements) > 0:
            self._universe.add_TopologyAttr('element', listelements)

        # Loop over bonds
        for ibond in self._universe.bonds:
            self._bond_list.append([ibond[0].index, ibond[1].index])

        # Dimensions
        with open(self._fnamepath, 'r') as f:
            self._unitcell[0, :] = [0.0, 0.0, 0.0]
            self._unitcell[1, :] = [0.0, 0.0, 0.0]
            self._unitcell[2, :] = [0.0, 0.0, 0.0]
            while True:
                line = f.readline()
                if not line:
                    break
                # if line.find("CRYST1") != -1:
                if line[0:6] == "CRYST1":
                    d = line.split()
                    if float(d[1]) != 0.0 and float(d[2]) != 0.0 and float(d[3]) != 0.0:
                        self._boxlength[0] = float(d[1])  # Angstroms
                        self._boxlength[1] = float(d[2])
                        self._boxlength[2] = float(d[3])
                        self._boxangle[0] = float(d[4]) * degtorad  # radians
                        self._boxangle[1] = float(d[5]) * degtorad
                        self._boxangle[2] = float(d[6]) * degtorad

                        if self._boxlength[0] != 0.0 and self._boxlength[1] != 0.0 and self._boxlength[2] != 0.0:
                            self._unitcell[0, :], self._unitcell[1, :], self._unitcell[2, :] = \
                                self.lengths_and_angles_to_box_vectors(self._boxlength[0], self._boxlength[1],
                                                                       self._boxlength[2], self._boxangle[0],
                                                                       self._boxangle[1], self._boxangle[2])
                    else:
                        self._boxlength = None
                        self._boxangle = None

        # Topology
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\n\t\t\t Start Creating topology. ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)
        self._topology = top.Topology(natoms=self._natoms, listbonds=self._bond_list)

        self._nmols = len(self._topology.get_nmols())
        self._mol_residue_list = list(self._universe.residues.resnames)
        ires = 0
        for item_list in self._universe.residues.indices:
            for idx in item_list:
                self._atom3d_resname[idx] = self._mol_residue_list[ires]
                self._atom3d_residue[idx] = ires + 1
            ires += 1
        # Get atom properties
        imol_idx = 0
        for imol in self._topology.get_nmols():
            for idx in imol:
                self._atom3d_imol[idx] = imol_idx
                self._atom3d_charge[idx] = 0.000
                try:
                    self._atom3d_element[idx] = self._universe.atoms.elements[idx]
                except exceptions.NoDataError:
                    self._atom3d_element[idx] = str(guess_atom_element(self._universe.atoms.names[idx]))
                self._atom3d_mass[idx] = round(self._universe.atoms.masses[idx], 2)

            imol_idx += 1

        # List
        charge_list = []
        elements_list = []
        mass_list = []
        isbackbone_list = []
        name_list = []
        type_list = []
        resname_list = []
        for idx in range(self._natoms):
            charge_list.append(self._atom3d_charge[idx])
            elements_list.append(self._atom3d_element[idx])
            name_list.append(self._universe.atoms.names[idx])
            type_list.append(self._universe.atoms.types[idx])
            mass_list.append(self._atom3d_mass[idx])
            resname_list.append(self._atom3d_resname[idx])
            isbackbone_list.append(self._atom3d_bfactor[idx])

        self._topology.set_elements(elements_list)
        self._topology.set_charge_mol(charge_list)
        self._topology.set_isbackbone(isbackbone_list)
        self._topology.set_mass(mass_list)
        self._topology.set_type(type_list)
        self._topology.set_name(name_list)
        self._topology.set_resname(resname_list)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t\t End Creating topology. ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)

        # Guess box length
        xmax = -1.0E12
        ymax = -1.0E12
        zmax = -1.0E12
        xmin = 1E20
        ymin = 1E20
        zmin = 1E20
        if self._boxlength is None and self._boxangle is None:
            self._boxlength = np.zeros(3, dtype=np.float32)
            self._boxangle = np.zeros(3, dtype=np.float32)
            for key, icoords in self._atom3d_xyz.items():
                if icoords[0] > xmax:
                    xmax = icoords[0]
                if icoords[1] > ymax:
                    ymax = icoords[1]
                if icoords[2] > zmax:
                    zmax = icoords[2]
                if icoords[0] < xmin:
                    xmin = icoords[0]
                if icoords[1] < ymin:
                    ymin = icoords[1]
                if icoords[2] < zmin:
                    zmin = icoords[2]

            self._boxlength[0] = xmax - xmin  # Angstroms
            self._boxlength[1] = ymax - ymin
            self._boxlength[2] = zmax - zmin

            self._boxangle[0] = 90.0 * degtorad  # Degrees
            self._boxangle[1] = 90.0 * degtorad
            self._boxangle[2] = 90.0 * degtorad

            self._unitcell[0, :], self._unitcell[1, :], self._unitcell[2, :] = \
                self.lengths_and_angles_to_box_vectors(self._boxlength[0], self._boxlength[1],
                                                       self._boxlength[2], self._boxangle[0],
                                                       self._boxangle[1], self._boxangle[2])

        # Assign bond orders ==================================================
        if self._assign_bo:
            self._topology.assign_bond_orders()

        # Check if there are hydrogen atoms bonded to the carbon atoms
        # in order to assign the correct
        # Loop over atoms
        listelements = list()
        name_list = []
        for iatom in self._universe.atoms:
            try:
                e = iatom.element
            except exceptions.NoDataError:
                e = str(guess_atom_element(iatom.name))
                listelements.append(e)
            if e == 'C' and \
                    'H' not in iatom.bonded_atoms.elements and \
                    set(iatom.bonded_atoms.elements) == {'C'}:

                if len(iatom.bonded_atoms.elements) == 1:
                    name_list.append("_CH3")
                elif len(iatom.bonded_atoms.elements) == 2:
                    name_list.append("_CH2")
                elif len(iatom.bonded_atoms.elements) == 3:
                    name_list.append("_CH")
                elif len(iatom.bonded_atoms.elements) == 4:
                    name_list.append("C")
                else:
                    print("\t\tWARNING: Index: {} Element: {} Number of bonds: {}".
                          format(iatom.id, e, len(iatom.bonded_atoms.elements)))
            else:
                # name_list.append(e) cJ 5-Mar-2024
                name_list.append(iatom.name)

        self._topology.set_name(name_list)

    # *************************************************************************
    def write_renumber_pdb(self, head_idx_atom, tail_idx_atom,
                           assign_residues_info=None, fnameout=None, isunwrap=False):

        """
        This method assign a new numbering scheme to the molecule. The algorithm is as follows:

            1. For each chain
                1.a. Find the path between the head and the tail atoms.
                     Atoms in the path are labelled as backbones. (tempfactor=1.0)
                1.b.

        :param head_idx_atom: A list with the index of the head atoms for each mol/chain.
        :param tail_idx_atom: A list with the index of the tail atom for each mol/chain.
        :param assign_residues_info
        :param fnameout
        :param isunwrap

        :return:
        """

        # ===========
        def create_subgraph_from_node(inode, isvisitedlist, backbonelist):

            """
            Subgraph to calculate longest distance. As the subgraph is searched using a bfs algorithm, the likely
            nodes for longest distance are the last ones in the list leaves

            :param inode:
            :param isvisitedlist:
            :return:
                A subgraph (dictionary) and the leaves of the subgraph
            """

            queue_sg = [inode]
            subgraph = defaultdict(list)
            put_on_subgraph = []
            leaves_sub = []
            root = inode

            while queue_sg:
                cnode = queue_sg.pop()
                put_on_subgraph.append(cnode)
                if cnode not in isvisitedlist:
                    neighs = topo.get_neighbours(cnode)
                    for ineigh in neighs:
                        if ineigh in isvisitedlist or ineigh in backbonelist:
                            continue
                        else:
                            subgraph[cnode].append(ineigh)
                            if len(neighs) == 1 and cnode != root:
                                leaves_sub.append(cnode)
                            if ineigh not in put_on_subgraph:
                                queue_sg.insert(0, ineigh)

            return subgraph, leaves_sub

        # ===========

        # ===========
        def longest_path(start, end, graph):

            """ Recursive solution """

            # Create a stack for DFS
            stack = list()
            depth = list()
            path = list()
            paths = list()
            stack.append(start)
            depth.append(start)

            while len(stack):
                # Pop a vertex from stack
                cur = stack.pop()
                lvl = depth.pop()

                while len(path) > 0 and path[-1] != lvl:
                    path.pop()

                path.append(cur)
                hasnew = False
                neighs = graph[cur]
                for ineigh in neighs:
                    if ineigh in path:
                        continue
                    if ineigh == end:
                        res = path
                        res.append(end)
                        paths.append(copy.deepcopy(res))
                        continue
                    else:
                        hasnew = True
                        stack.append(ineigh)
                        depth.append(cur)
                if not hasnew:
                    path.pop()

            return paths

        # ===========

        # ===========
        def write_new_pdb(gnew, fnameout2=None):

            if fnameout2 is None:
                fnameout2 = os.path.splitext(self._fname)[0] + "_renumber.pdb"
            else:
                fnameout2 = fnameout2

            try:
                remark_line = self._universe.trajectory.remarks
                remark_line[0] += ". Renumbered with topology"
            except IndexError:
                remark_line = [" Renumbered with topology"]
            remark_line += ["Backbone atoms --> Beta field == 0.00"]
            remark_line += ["Branch   atoms --> Beta field >= 1.00"]
            remark_line += ["Head     atoms --> Occupanccy field == 1.00"]
            remark_line += ["Tail     atoms --> Occupanccy field == 2.00"]
            remark_line += ["Middle   atoms --> Occupanccy field == 0.00"]

            radtodeg = 180 / np.pi
            with open(fnameout2, 'w') as f:
                for ii in remark_line:
                    f.write("REMARK     {0:s}\n".format(ii))

                f.write("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}"
                        "{3:7.2f}{4:7.2f}{5:7.2f} "
                        "{6:<11s}{7:4d}\n".format(self._boxlength[0],
                                                  self._boxlength[1],
                                                  self._boxlength[2],
                                                  self._boxangle[0] * radtodeg,
                                                  self._boxangle[1] * radtodeg,
                                                  self._boxangle[2] * radtodeg,
                                                  self._spacegroup, self._zvalue))

                for aidx_new, aidx_old in gnew.items():
                    # cJ: Bug for systems with more than 100000 atoms
                    aidx_new_normal = (aidx_new + 1) % 100000

                    f.write("HETATM{0:5d} {1:<4s}{2:<1s}{3:<4s}"
                            "{4:1s}{5:4d}{6:1s}   "
                            "{7:8.3f}{8:8.3f}{9:8.3f}{10:6.2f}"
                            "{11:6.2f}      {12:<4s}{13:>2s}\n".format(aidx_new_normal,
                                                                       self._topology._names[aidx_old],
                                                                       '',
                                                                       self._atom3d_resname[aidx_old],
                                                                       '',
                                                                       self._atom3d_residue[aidx_old],
                                                                       '',
                                                                       self._atom3d_xyz[aidx_old][0],
                                                                       self._atom3d_xyz[aidx_old][1],
                                                                       self._atom3d_xyz[aidx_old][2],
                                                                       self._atom3d_occupancy[aidx_old],
                                                                       self._atom3d_bfactor[aidx_old],
                                                                       '',
                                                                       self._atom3d_element[aidx_old]))
                f.write("END\n")
                bondlist = self._topology.get_allbonds()

                dd = defaultdict(list)
                res = dict((v, k) for k, v in gnew.items())
                for ibond in bondlist:
                    iat_old = ibond[0]
                    jat_old = ibond[1]
                    iat = res[iat_old]
                    jat = res[jat_old]
                    dd[iat].append(jat)
                    dd[jat].append(iat)

                # cJ: Bug for systems with more than 100000 atoms
                if self._natoms < 100000:
                    for k in sorted(dd.keys()):
                        cc = [k + 1]
                        for ival in dd[k]:
                            cc.append(ival + 1)
                        conect = ["{0:5d}".format(entry) for entry in cc]
                        conect = "".join(conect)
                        f.write("CONECT{0}\n".format(conect))

        # =====================================
        # Dictionary g_new = {new_idx: old_idx, ...}
        g_new = defaultdict()
        idx_new = 0
        imol = 0

        for iatom, _ in self._atom3d_occupancy.items():
            self._atom3d_occupancy[iatom] = 0.0
            self._atom3d_bfactor[iatom] = 1.0

        for ihead in head_idx_atom:
            now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            print("Mol = {}, Head = {}, Tail = {} ({})".format(imol, ihead, tail_idx_atom[imol], now))
            visitedlist = []
            queue = [ihead]
            topo = self._topology
            itail = tail_idx_atom[imol]

            self._atom3d_occupancy[ihead] = 1.0
            self._atom3d_occupancy[itail] = 2.0

            # backbone_list = self._topology.find_all_paths(ihead, itail)[0]
            backbone_list = self._topology.find_all_paths_iterative(ihead, itail)
            for item in backbone_list:
                self._atom3d_bfactor[item] = 0.0

            # Create dictionary with the new numbers
            while queue:
                idx_old = queue.pop(0)
                if idx_old in visitedlist:
                    continue
                g_new[idx_new] = idx_old
                idx_new += 1
                visitedlist.append(idx_old)
                neigh_old = topo.get_neighbours(idx_old)
                dict_neigh = defaultdict(list)
                for ineigh_old in neigh_old:
                    an = atomic_number[self._atom3d_element[ineigh_old]]
                    if ineigh_old not in visitedlist:
                        dict_neigh[an].append(ineigh_old)

                for ikey, values in sorted(dict_neigh.items()):
                    if len(values) > 0 and ikey != 1:
                        # Get the bb_atom and create a list with the br_atoms (notbackbone)
                        br_atom = list()
                        for item in values:
                            if item in backbone_list:
                                bb_atom = item
                                # Backbone atom is last
                                if bb_atom not in queue:
                                    queue.append(bb_atom)
                            else:
                                br_atom.append(item)
                                # Branch atoms are first ###############
                                if len(br_atom) > 1:
                                    lengths_path = defaultdict()
                                    for item2 in br_atom:
                                        sg, leaves = create_subgraph_from_node(item2, visitedlist, backbone_list)
                                        # Only the  last five nodes are inspected.
                                        max_path = 0
                                        for ileave in leaves[-5:]:
                                            ll = longest_path(item2, ileave, sg)
                                            if len(ll[0]) > max_path:
                                                max_path = len(ll[0])
                                        lengths_path[item2] = max_path
                                    for ikey2, value in sorted(lengths_path.items(), key=lambda item3: item3[1]):
                                        if ikey2 not in queue:
                                            queue.append(ikey2)
                                else:
                                    item = br_atom[0]
                                    queue.append(item)
                                    # Put in the queue all branch atoms before to the backbone atom
                                    sg, leaves = create_subgraph_from_node(item, visitedlist, backbone_list)
                                    for ileave in leaves:
                                        ll = longest_path(item, ileave, sg)
                                        for ipath in ll:
                                            for item2 in ipath:
                                                if item2 not in queue:
                                                    queue.append(item2)

                    # Hydrogen atoms are always assigned as branch atoms
                    else:
                        for item in values:
                            if item not in queue:
                                queue.append(item)
                # Put the backbone atom in the last position of the queue:
                for item in queue:
                    if item in backbone_list:
                        queue.remove(item)
                        queue.append(item)
                        break
            imol += 1

        # Unwrap coordinates
        if isunwrap:
            coords = np.array(list(self._atom3d_xyz.values()))
            tupleout = self._topology.get_array_mols_neigh()
            nmols_array = tupleout[0]
            l_neigh_array = tupleout[1]
            c = unwrap(coords, nmols_array, l_neigh_array, self._boxlength)
            for idx, item in enumerate(c):
                self._atom3d_xyz[idx] = item

        # Write PDB before to change topology
        if assign_residues_info is None:
            write_new_pdb(g_new, fnameout2=fnameout)
            fnamegro = os.path.splitext(fnameout)[0] + ".gro"
            self.write_gro(filename_gro=fnamegro)

        # Change self._topology object for PSF
        elements_list = []
        charge_list = []
        mass_list = []
        isbackbone_list = []
        iatch_list = []
        name_list = []
        type_list = []
        resname_list = []
        xlist = []
        ylist = []
        zlist = []
        occlist = []
        bfactorlist = []
        for aidx_newt, aidx_oldt in g_new.items():
            elements_list.append(self._topology.elements[aidx_oldt])
            charge_list.append(self._topology.charge[aidx_oldt])
            mass_list.append(self._topology.mass[aidx_oldt])
            isbackbone_list.append(self._atom3d_bfactor[aidx_oldt])
            iatch_list.append(self._topology._iatch[aidx_oldt])
            name_list.append(self._universe.atoms.names[aidx_oldt])
            type_list.append(self._universe.atoms.types[aidx_oldt])
            resname_list.append(self._atom3d_resname[aidx_oldt])
            xlist.append(self._atom3d_xyz[aidx_oldt][0])
            ylist.append(self._atom3d_xyz[aidx_oldt][1])
            zlist.append(self._atom3d_xyz[aidx_oldt][2])
            occlist.append(self._atom3d_occupancy[aidx_oldt])
            bfactorlist.append(self._atom3d_bfactor[aidx_oldt])

        self._topology.set_elements(elements_list)
        self._topology.set_charge_mol(charge_list)
        self._topology.set_mass(mass_list)
        self._topology.set_isbackbone(isbackbone_list)
        for idx, ich in enumerate(iatch_list):
            self._topology._iatch[idx] = ich
        self._topology.set_type(type_list)
        self._topology.set_name(name_list)
        self._topology.set_resname(resname_list)

        graph_dict_new = defaultdict(list)
        map_old_to_new = defaultdict()
        for aidx_newt, aidx_oldt in g_new.items():
            map_old_to_new[aidx_oldt] = aidx_newt
        for oldkey, olditem in self._topology._graphdict.items():
            lll = [map_old_to_new[i] for i in olditem]
            graph_dict_new[map_old_to_new[oldkey]] = lll
        self._topology._graphdict = graph_dict_new
        self._topology._bonds = []

        for idx, item in enumerate(elements_list):
            self._atom3d_element[idx] = item
            self._atom3d_xyz[idx][0] = xlist[idx]
            self._atom3d_xyz[idx][1] = ylist[idx]
            self._atom3d_xyz[idx][2] = zlist[idx]
            self._atom3d_occupancy[idx] = occlist[idx]
            self._atom3d_bfactor[idx] = bfactorlist[idx]
            self._atom3d_isbackbone[idx] = int(bfactorlist[idx])

        if assign_residues_info is not None:
            self.assign_residues_chains(assign_residues_info, fout=fnameout)

        return g_new

    # *************************************************************************
    def write_label_pdb(self, head_idx_atom, tail_idx_atom,
                        assign_residues_info=None, fnameout=None, isunwrap=False):

        """
        This method label a pdb file

        :param head_idx_atom: A list with the index of the head atoms for each mol/chain.
        :param tail_idx_atom: A list with the index of the tail atom for each mol/chain.
        :param assign_residues_info
        :param fnameout
        :param isunwrap

        :return:
        """

        # ===========
        def write_new_pdb_norenumber(gnew, fnameout2=None):

            if fnameout2 is None:
                fnameout2 = os.path.splitext(self._fname)[0] + "_labeled.pdb"
            else:
                fnameout2 = fnameout2

            try:
                remark_line = self._universe.trajectory.remarks
                remark_line[0] += ". Labeled with topology"
            except IndexError:
                remark_line = [" Labeled with topology"]
            remark_line += ["Backbone atoms --> Beta field == 0.00"]
            remark_line += ["Branch   atoms --> Beta field >= 1.00"]
            remark_line += ["Head     atoms --> Occupanccy field == 1.00"]
            remark_line += ["Tail     atoms --> Occupanccy field == 2.00"]
            remark_line += ["Middle   atoms --> Occupanccy field == 0.00"]

            radtodeg = 180 / np.pi
            with open(fnameout2, 'w') as f:
                for ii in remark_line:
                    f.write("REMARK     {0:s}\n".format(ii))

                f.write("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}"
                        "{3:7.2f}{4:7.2f}{5:7.2f} "
                        "{6:<11s}{7:4d}\n".format(self._boxlength[0],
                                                  self._boxlength[1],
                                                  self._boxlength[2],
                                                  self._boxangle[0] * radtodeg,
                                                  self._boxangle[1] * radtodeg,
                                                  self._boxangle[2] * radtodeg,
                                                  self._spacegroup, self._zvalue))

                for aidx_new, aidx_old in gnew.items():
                    # cJ: Bug for systems with more than 100000 atoms
                    aidx_new_normal = (aidx_new + 1) % 100000

                    f.write("HETATM{0:5d} {1:<4s}{2:<1s}{3:<4s}"
                            "{4:1s}{5:4d}{6:1s}   "
                            "{7:8.3f}{8:8.3f}{9:8.3f}{10:6.2f}"
                            "{11:6.2f}      {12:<4s}{13:>2s}\n".format(aidx_new_normal,
                                                                       self._topology._names[aidx_old],
                                                                       '',
                                                                       self._atom3d_resname[aidx_old],
                                                                       '',
                                                                       self._atom3d_residue[aidx_old],
                                                                       '',
                                                                       self._atom3d_xyz[aidx_old][0],
                                                                       self._atom3d_xyz[aidx_old][1],
                                                                       self._atom3d_xyz[aidx_old][2],
                                                                       self._atom3d_occupancy[aidx_old],
                                                                       self._atom3d_bfactor[aidx_old],
                                                                       '',
                                                                       self._atom3d_element[aidx_old]))
                f.write("END\n")
                bondlist = self._topology.get_allbonds()

                dd = defaultdict(list)
                res = dict((v, k) for k, v in gnew.items())
                for ibond in bondlist:
                    iat_old = ibond[0]
                    jat_old = ibond[1]
                    iat = res[iat_old]
                    jat = res[jat_old]
                    dd[iat].append(jat)
                    dd[jat].append(iat)

                # cJ: Bug for systems with more than 100000 atoms
                if self._natoms < 100000:
                    for k in sorted(dd.keys()):
                        cc = [k + 1]
                        for ival in dd[k]:
                            cc.append(ival + 1)
                        conect = ["{0:5d}".format(entry) for entry in cc]
                        conect = "".join(conect)
                        f.write("CONECT{0}\n".format(conect))

        # =====================================
        # Dictionary g_new = {new_idx: old_idx, ...}
        g_new = defaultdict()
        for i in range(self._natoms):
            g_new[i] = i
        imol = 0

        for iatom, _ in self._atom3d_occupancy.items():
            self._atom3d_occupancy[iatom] = 0.0
            self._atom3d_bfactor[iatom] = 1.0
            self._atom3d_isbackbone[iatom] = 1.0

        for ihead in head_idx_atom:

            itail = tail_idx_atom[imol]

            self._atom3d_occupancy[ihead] = 1.0
            self._atom3d_occupancy[itail] = 2.0
            try:
                # backbone_list = self._topology.find_all_paths(ihead, itail)[0]
                backbone_list = self._topology.find_all_paths_iterative(ihead, itail)
            except IndexError:
                backbone_list = None
                m = "\n\t\t ERROR: There is a problem with the head {} and tail {} atoms in mol {}". \
                    format(ihead, itail, imol)
                print(m) if self._logger is None else self._logger.error(m)
                exit()
            for item in backbone_list:
                self._atom3d_bfactor[item] = 0.0
                self._atom3d_isbackbone[item] = 0.0
            imol += 1

        # Unwrap coordinates
        if isunwrap:
            coords = np.array(list(self._atom3d_xyz.values()))
            tupleout = self._topology.get_array_mols_neigh()
            nmols_array = tupleout[0]
            l_neigh_array = tupleout[1]
            c = unwrap(coords, nmols_array, l_neigh_array, self._boxlength)
            for idx, item in enumerate(c):
                self._atom3d_xyz[idx] = item

        # Write PDB before to change topology
        if assign_residues_info is None:
            write_new_pdb_norenumber(g_new, fnameout2=fnameout)
            fnamegro = os.path.splitext(fnameout)[0] + ".gro"
            self.write_gro(filename_gro=fnamegro)

        if assign_residues_info is not None:
            self.assign_residues_chains(assign_residues_info, fout=fnameout)

        return g_new

    # *************************************************************************
    def read_head_tail_info(self, fpathdat):

        """
        Read a file text contianing the information about head and tails in the molecular
        system. The format is as follows

        nmols: <Integer number of molecules/chains>
        # ich idx.head idx.tail (indexes begins at zero)
        0  4  7
        1  20 45
        (... as many mols as in the molecular system)

        :param fpathdat:
        :return:
        """

        fname_path, fname = os.path.split(fpathdat)[0:2]

        with open(fpathdat, 'r') as f:
            # Read all lines
            lines = f.readlines()
            # Remove lines with comments and return carriage
            lines = [iline.strip("\n") for iline in lines if not iline.strip().startswith("#")]
            # Get line with nmols keyword
            nmols = [iline for iline in lines if iline.strip().startswith("nmols")]
            if len(nmols) != 1:
                print("ERROR: nmols keyword must appear only once in the file {}".format(fpathdat))
                exit()
            try:
                nmols = int(nmols[0].split()[-1])
                if nmols != self._nmols:
                    m = "\n\t\t ERROR: The 'nmols' keyword in {}\n".format(fpathdat)
                    m += "\t\t must match the number of molecules\n"
                    m += "\t\t in the PDB file ({})\n".format(self._fnamepath)
                    m += "\t\t nmols = {}\n".format(nmols)
                    m += "\t\t Number of molecules = {}\n".format(self._nmols)
                    m += "\t\t This is probably caused by the overlapping of atoms from different molecules.\n"
                    m += "\t\t Aborting!!!!".format(fpathdat, self._fnamepath)
                    print(m) if self._logger is None else self._logger.error(m)
                    # Summarize chains in a file
                    natch_dict = defaultdict(list)
                    with open("error_summarize_chains.dat", 'w') as fchain:
                        line1 = "# Expected number of chains: {}\n".format(nmols)
                        line1 += "# Number_of_Atoms Number_of_chains\n"
                        for imol in self._topology.get_nmols():
                            natch_dict[len(imol)].append(imol)
                        for key, values in natch_dict.items():
                            line1 += "{0:6d} {1:6d}\n".format(key, len(values))
                        fchain.writelines(line1)
                    exit()
                elif len(lines) != nmols + 1:
                    print("ERROR: The number of lines must be equal to the to the number of mols\n "
                          "in the DAT file ({})".format(fpathdat))
                    exit()
            except ValueError as e:
                print("ERROR: The number of mols must be an integer in {}".format(fpathdat))
                print("EXCEPTION as {}".format(e))
                exit()

            try:
                for item in lines[1:]:
                    ich, ihead, itail = item.split()
                    self._heads.append(int(ihead))
                    self._tails.append(int(itail))
            except ValueError:
                print("ERROR: There are problems with the lines containing head and tail atoms ({})".format(fpathdat))
                exit()
                # return [], [] #cJ 6-Mar-2024

        return self._heads, self._tails

    # *************************************************************************
    def assign_residues_chains(self, fpathdat, fout=None):

        """
        This method assign the name and number to the diferent residues in the
        system. The information for the residues are given in an external file.
        The PDB file must be numbering using the criterium implemented in write_renumber_pdb
        ChainID is used to label the type of molecules

        The format of the file is the following:

        All-mols-equals = [True/False]
        Nmols = {integer; 1 if all-mols-equals is True}
        Mol {from, int} {to, int}  # (Example Mol 0 9 for the first ten chains) (Mol 0 0 just for the first chain)
        RES {from, int} {to, int} {name_of_residue, str} {number of atoms in the residue}

        :return:
        """

        fname_path, fname = os.path.split(fpathdat)[0:2]
        dict_residues_read = defaultdict()

        atom_kind_molecule_label = defaultdict()

        with open(fpathdat, 'r') as f:
            # Read all lines
            lines = f.readlines()
            # Remove lines with comments and return carriage
            lines = [iline.strip("\n") for iline in lines if not iline.strip().startswith("#")]
            # print(lines)

            # Setup dictionary residues. Check data n the file
            try:
                start_idx = lines.index("<RESIDUES>")
                end_idx = lines.index("</RESIDUES>")
            except ValueError:
                print("ERROR!!!!. <RESIDUES> ... </RESIDUES> labels must exist in the file {}".format(fpathdat))
                return False
            try:
                _ = int(lines[start_idx + 1])
            except ValueError:
                print("ERROR!!!!. Number of residues must be an integer in the file {}".format(fpathdat))
                return False

            for idx in range(start_idx + 2, end_idx):
                iline = lines[idx].split()
                dict_residues_read[iline[1]] = int(iline[2])

            # Read composition
            try:
                start_res_idx = lines.index("<COMPOSITION>")
                end_res_idx = lines.index("</COMPOSITION>")
            except ValueError:
                print("ERROR!!!!. <COMPOSITION> ... </COMPOSITION> labels must exist in the file {}".format(fpathdat))
                return False

            # Check number of residues in composition
            nkinds = int(lines[start_res_idx + 1])
            idx = start_res_idx + 2
            while idx < end_res_idx:
                for ikind in range(nkinds):
                    count = 0
                    count_lines = 0
                    _, _, nr, _ = [ii for ii in lines[idx].split()]
                    nr = int(nr)
                    idx += 1
                    for ires in range(nr):
                        try:
                            aa, start, end = [ii for ii in lines[idx].split()]
                            count += (int(end) - int(start)) + 1
                            count_lines += 1
                            idx += 1
                        except ValueError:
                            break
                    if count_lines != nr:
                        print("ERROR!!!!. <COMPOSITION> ... </COMPOSITION> "
                              "number of residues incorrect in chain kind {}. File {}".format(
                            ikind, fpathdat))
                        print("ERROR!!!!. Counted residue lines {} of {}".format(count_lines, nr))
                        exit()

            idx_line = start_res_idx + 1
            n_type_mols = int(lines[idx_line])
            idx_line += 1
            iat = 0
            jjres = 0
            itype = 0
            for itype in range(0, n_type_mols):
                start_mol_idx, end_mol_idx, nres_mol, kind_molecule_label = lines[idx_line].split()
                # print(start_mol_idx, end_mol_idx, nres_mol, kind_molecule_label)
                for imol in range(int(start_mol_idx), int(end_mol_idx) + 1):
                    for ires in range(1, int(nres_mol) + 1):
                        name_res = str(lines[idx_line + ires].split()[0])
                        start_res_idx = int(lines[idx_line + ires].split()[1])
                        end_res_idx = int(lines[idx_line + ires].split()[2])
                        for jres in range(start_res_idx, end_res_idx + 1):
                            # print(iat + dict_residues_read[name_res])
                            try:
                                for j in range(iat, iat + dict_residues_read[name_res]):
                                    self._atom3d_resname[j] = name_res
                                    self._atom3d_residue[j] = jjres + 1
                                    atom_kind_molecule_label[j] = kind_molecule_label
                                #   print(idx_line+ires, lines[idx_line+ires], name_res, start_res_idx, end_res_idx)
                                iat += dict_residues_read[name_res]
                            except KeyError:
                                print("ERROR!!!!. Name of residue {} is not defined in residue dat file "
                                      "{}".format(name_res, fpathdat))
                                return False
                            jjres += 1
                idx_line += ires + 1

            # Check if all atoms are assigned to a molecule
            if len(atom_kind_molecule_label) != self._natoms:
                print("\nERROR!!!!. <COMPOSITION> ... </COMPOSITION> "
                      "Not all atoms have been correctly assigned to a chain.\n "
                      "Usually, this is due to an incorrect number of residues in file {}.\n"
                      "Check that file.\n"
                      "natoms         = {}\n"
                      "atoms assigned = {}\n".
                      format(fpathdat, self._natoms, len(atom_kind_molecule_label)))
                exit()

            resname_list = []
            for ikey, item in self._atom3d_resname.items():
                resname_list.append(item)
            self._topology.set_resname(resname_list)

            fnamegro = os.path.splitext(fout)[0] + ".gro"
            fnamepdb = os.path.splitext(fout)[0] + ".pdb"
            self.write_gro(filename_gro=fnamegro)
            self.write_pdb(filename_pdb=fnamepdb, atom_kind_molecule_label=atom_kind_molecule_label)

            return True

    # *************************************************************************
    def get_backbone_headtail_files_for_analysis(self):

        d_bb = defaultdict(list)
        d_br = defaultdict(list)

        # for iatom, ivalue in self._atom3d_bfactor.items():
        #     ich = self._topology._iatch[iatom]
        #     if ivalue == 0.0:
        #         d_bb[ich].append(iatom)
        #     else:
        #         d_br[ich].append(iatom)

        with open('backbone_idx.dat', 'w') as f:
            ich = 0
            for ihead in self._heads:
                iatom = ihead
                f.writelines("[mol{}]\n".format(ich))
                f.writelines("{}\n".format(iatom))
                isvisited = list()
                isvisited.append(iatom)
                d_bb[ich].append(iatom)

                while iatom != self._tails[ich]:
                    neigh_bb = [i for i in self._topology._graphdict[iatom] if self._atom3d_bfactor[i] == 0.0]
                    for jatom in neigh_bb:
                        if jatom not in isvisited:
                            f.writelines("{}\n".format(jatom))
                            isvisited.append(jatom)
                            d_bb[ich].append(jatom)
                            iatom = jatom
                            continue

                ich += 1

        # # FOR DEBUG PURPOSES
        # with open('backbone_idx_vmd.dat', 'w') as f:
        #     for ich in d_bb:
        #         f.writelines("[mol{}]\n".format(ich))
        #         for ivalue in d_bb[ich]:
        #             f.writelines("{} ".format(ivalue))
        #         f.writelines("\n")

        with open('branch_idx.dat', 'w') as f:
            for ich in range(len(self._topology._nmols)):
                f.writelines("[mol{}]\n".format(ich))
                for ivalue in self._topology._nmols[ich]:
                    if ivalue not in d_bb[ich]:
                        f.writelines("{}\n".format(ivalue))
                        d_br[ich].append(ivalue)

        with open('listend2end.dat', 'w') as f:
            f.writelines("# ich head tail\n".format(ich))
            for ich in range(self._nmols):
                f.writelines("{} {} {}\n".format(ich, self._heads[ich], self._tails[ich]))

        # # FOR DEBUG PURPOSES
        # with open('listend2end_heads_vmd.dat', 'w') as f:
        #     for ich in range(self._nmols):
        #         f.writelines("{} ".format(self._heads[ich]))
        #     f.writelines("\n")
        #
        # with open('listend2end_tails_vmd.dat', 'w') as f:
        #     for ich in range(self._nmols):
        #         f.writelines("{} ".format(self._tails[ich]))
        #     f.writelines("\n")

    # *************************************************************************
    def get_bond_list(self):

        return self._bond_list

    # *************************************************************************
    def _removeterm_and_savepdb(self, fnamepdb):

        idx_old = 1
        idx_new = -1
        idx_old_to_new = defaultdict()
        nter = 0
        nnter = -2
        lines_new_pdb = list()
        with open(fnamepdb, 'r') as fpdb:
            lines_old_pdb = fpdb.readlines()
            for iline in lines_old_pdb:
                if iline.count("ATOM") == 0 and \
                        iline.count("HETATM") == 0 and \
                        iline.count("CONECT") == 0 and \
                        iline.count("TER") == 0:
                    lines_new_pdb.append(iline)
                elif iline.count("ATOM") != 0 or iline.count("HETATM") != 0:
                    idx_old = int(iline[6:11])
                    idx_new = idx_old - nter
                    idx_old_to_new[idx_old] = idx_new
                    inewline = iline[0:6] + "{0:>5s}".format(str(idx_new)) + iline[11:22] \
                        + "{0:>4s}".format(str(nter + 1)) + iline[26:]
                    lines_new_pdb.append(inewline)
                elif iline.count("TER") != 0:
                    idx_old = int(iline[6:11])
                    idx_new = idx_old - nter
                    idx_old_to_new[idx_old] = idx_new
                    nter += 1
                elif iline.count("CONECT"):
                    tokens = []
                    idx_token = 6
                    while True:
                        tt = iline[idx_token:idx_token+5]
                        if tt.count("\n") != 0:
                            break
                        tokens.append(int(iline[idx_token:idx_token+5]))
                        idx_token = idx_token+5
                    inewline = "CONECT"
                    for item_old in tokens:
                        try:
                            inewline += "{0:>5d}".format(idx_old_to_new[item_old])
                        except KeyError:
                            nnter += 1
                            inewline += "{0:>5d}".format(item_old-nnter)

                    inewline += "\n"
                    lines_new_pdb.append(inewline)
                else:
                    lines_new_pdb.append(iline)

        basefile = os.path.split(fnamepdb)[-1]
        basename = os.path.splitext(basefile)[0]

        fnamepdbmod = basename + "_mod.pdb"
        with open(fnamepdbmod, 'w') as foutpdb:
            foutpdb.writelines(lines_new_pdb)

        self._fnamepath = fnamepdbmod

        return fnamepdbmod

    # # *************************************************************************
    # def check_read_structure(self):
    #
    #     """
    #     Review certain aspects of the structure obtained from the PDB file.
    #     """
    #
    #     # Check for number of molecules and head/tail
    #     if len(self._heads) != 0 and len(self._heads) != len(self._tails):
    #         m = "\n\t\t ERROR: It appears that information about 'heads and tails' " \
    #             "is available in the PDB file (occupancy column)\n"
    #         m += "\t\t ERROR: The count of heads ({}) does not match the count of tail numbers ({})\n".\
    #             format(len(self._heads), len(self._tails))
    #         m += "\t\t Aborting!!!!!"
    #         print(m) if self._logger is None else self._logger.error(m)
    #         exit()
    #
    #     if len(self._heads) != 0 and len(self._heads) != self._nmols:
    #         m = "\n\t\t ERROR: It appears that information about 'heads and tails' " \
    #             "is available in the PDB file (occupancy column)\n"
    #         m += "\t\t ERROR: The count of heads ({}) does not match the count of " \
    #              "molecules after guessing the topology ({})\n".\
    #             format(len(self._heads), self._nmols)
    #         m += "\t\t ERROR: Number of bonds ({})\n".\
    #             format(self._nbonds)
    #         m += "\t\t ERROR: Probable overlap from atoms.\n"
    #         m += "\t\t ERROR: Use a topology file to get the bond list or \n" \
    #              "try to minimize/equilibrate before to use topology.\n"
    #         m += "\t\t Aborting!!!!!"
    #         print(m) if self._logger is None else self._logger.error(m)
    #         # Summarize chains in a file
    #         natch_dict = defaultdict(list)
    #         with open("error_summarize_chains.dar", 'w') as fchain:
    #             line1 = "# Expected number of chains: {}\n".format(self._nmols)
    #             line1 += "# Number_of_Atoms Number_of_chains\n"
    #             for imol in self._topology.get_nmols():
    #                 natch_dict[len(imol)].append(imol)
    #             for key, values in natch_dict.items():
    #                 line1 += "{0:6d} {1:6d}\n".format(key, len(values))
    #             fchain.writelines(line1)
    #
    #         exit()
    #
    #     if len(self._universe.atoms.names) != len(self._topology._names):
    #         m = "\n\t\tERROR!!!!!: The number of atoms ({}) is not equal to the dimension of names ({})".\
    #               format(len(self._universe.atoms.names), len(self._topology._names))
    #         m += "\t\tERROR!!!!!: Probably the system is overlaped. Try minimize or equilibrate before.\n"
    #         m += "\t\t Aborting!!!!!"
    #         print(m) if self._logger is None else self._logger.error(m)
    #         exit()
