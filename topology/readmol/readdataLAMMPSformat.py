import numpy as np
import topology as top
from collections import defaultdict
from MDAnalysis import Universe, exceptions
from MDAnalysis.topology.core import guess_atom_element

from topology.readmol.readbaseformat import ReadBaseFormat


class ReadDataLammpsFormat(ReadBaseFormat):

    # *************************************************************************
    def __init__(self, filenamepath, filemap, assign_bondorders=False):

        """

        :param filenamepath:
        :param filenames: A text file mapping the LAMMPS type with the name and element
        :param assign_bondorders:
        """

        super().__init__(filenamepath=filenamepath, assign_bondorders=assign_bondorders)

        self._filemap = filemap
        self._map_type_name = defaultdict()
        self._map_type_element = defaultdict()
        with open(self._filemap) as f:
            lines = f.readlines()
            ntypes = int(lines[0])
            for idx in range(1, len(lines)):
                iline = lines[idx].strip()
                if iline.startswith("#"):
                    continue
                itype, name, element = iline.split()
                self._map_type_name[int(itype)] = name
                self._map_type_element[int(itype)] = element

        if filenamepath is not None:
            self.read_data()

    # *************************************************************************
    def __copy__(self):

        t = ReadDataLammpsFormat(filenamepath=self._fnamepath, assign_bondorders=self._assign_bo)

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
            if isinstance(getattr(self,key), np.ndarray):
                par = np.array_equal(getattr(self,key), getattr(other,key))
                res = res and par
            elif isinstance(getattr(self,key), top.Topology):
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
    def read_data(self):

        """
        Base function to read the DATA format from LAMMPS
        Returns:

        """

        degtorad = np.pi/180.0

        # MDAnalysis Universe
        self._universe = Universe(self._fnamepath)
        self._natoms = self._universe.atoms.n_atoms
        self._dimensions = self._universe.dimensions
        try:
            self._nbonds = len(self._universe.bonds)
        except exceptions.NoDataError:
            self._universe = Universe(self._fnamepath, guess_bonds=True)
            self._nbonds = len(self._universe.bonds)

        self._nres = len(self._universe.residues)

        # As data format has not information about the atom element and name
        # is needed to provide a file to correspond the atom type (integer) of LAMMPS
        # with the name of the atom
        # Example:
        #   4
        #   1 C
        #   2 H
        #   3 C
        #   4 H
        # Assign the atom name in the universe
        name_list = []
        element_list = []
        for iatom in self._universe.atoms:
            name_list.append(self._map_type_name[int(iatom.type)])
            element_list.append(self._map_type_element[int(iatom.type)])

        # Resname 'UNK'
        resname_list = ['UNK']*self._nres

        self._universe.add_TopologyAttr("name", name_list)
        self._universe.add_TopologyAttr("element", element_list)
        self._universe.add_TopologyAttr("resname", resname_list)

        # Loop over atoms
        listelements = list()
        for iatom in self._universe.atoms:
            try:
                e = iatom.element
            except exceptions.NoDataError:
                e = str(guess_atom_element(iatom.name))
                listelements.append(e)
            idx = iatom.id - 1
            self._atom3d_element[idx] = e
            self._atom3d_xyz[idx] = list(iatom.position)
            self._atom3d_molname[idx] = iatom.resname
            self._atom3d_isbackbone[idx] = 0.0

        # Loop over bonds
        bond_list = []
        for ibond in self._universe.bonds:
            bond_list.append([ibond[0].index, ibond[1].index])

        # Dimensions
        with open(self._fnamepath, 'r') as f:
            self._unitcell[0, :] = [self._universe.dimensions[0], 0.0, 0.0]
            self._unitcell[1, :] = [0.0, self._universe.dimensions[1], 0.0]
            self._unitcell[2, :] = [0.0, 0.0, self._universe.dimensions[2]]
            while True:
                line = f.readline()
                if not line:
                    break
                if line.find("xy") != -1:
                    d = line.split()
                    self._unitcell[0, 1] = float(d[0])
                    self._unitcell[1, 0] = float(d[0])
                    self._unitcell[0, 2] = float(d[1])
                    self._unitcell[2, 0] = float(d[1])
                    self._unitcell[1, 2] = float(d[2])
                    self._unitcell[2, 1] = float(d[2])

            self._boxlength[0] = self._universe.dimensions[0]      # Angstroms
            self._boxlength[1] = self._universe.dimensions[1]
            self._boxlength[2] = self._universe.dimensions[2]
            self._boxangle[0] = self._universe.dimensions[3]*degtorad      # radians
            self._boxangle[1] = self._universe.dimensions[4]*degtorad
            self._boxangle[2] = self._universe.dimensions[5]*degtorad

        # Topology
        self._topology = top.Topology(natoms=self._natoms, listbonds=bond_list)
        self._nmols = len(self._topology.get_nmols())
        self._mol_residue_list = list(self._universe.residues.resnames)
        ires = 0
        for item_list in self._universe.residues.indices:
            for idx in item_list:
                self._atom3d_resname[idx] = self._mol_residue_list[ires]
                self._atom3d_residue[idx] = ires+1
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
            isbackbone_list.append(0)

        self._topology.set_elements(elements_list)
        self._topology.set_charge_mol(charge_list)
        self._topology.set_isbackbone(isbackbone_list)
        self._topology.set_mass(mass_list)
        self._topology.set_type(type_list)
        self._topology.set_name(name_list)
        self._topology.set_resname(resname_list)

        # Assign bond orders ==================================================
        if self._assign_bo:
            self._topology.assign_bond_orders()

        return None




