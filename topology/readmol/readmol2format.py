from MDAnalysis import Universe
import topology as top
from collections import defaultdict
from topology.readmol.readbaseformat import ReadBaseFormat


class ReadMol2Format(ReadBaseFormat):

    # *************************************************************************
    def __init__(self, filenamepath, assign_bondorders=False):

        super().__init__(filenamepath=filenamepath, assign_bondorders=assign_bondorders)

        self._atom3d_bfactor = defaultdict()
        self._atom3d_occupancy = defaultdict()
        self._heads = list()
        self._tails = list()

        self.__remove_comments
        self.read_mol2()

    # *************************************************************************
    def read_mol2(self):

        """
        Base function to read the mol2 format
        Returns:

        """

        # MDAnalysis Universe
        self._universe = Universe(self._fnamepath)
        self._natoms = self._universe.atoms.n_atoms
        self._nbonds = len(self._universe.bonds)
        self._nres = len(self._universe.residues)

        # Loop over atoms
        for iatom in self._universe.atoms:
            e = iatom.element
            idx = iatom.id - 1
            self._atom3d_element[idx] = e
            self._atom3d_xyz[idx] = list(iatom.position)
            self._atom3d_molname[idx] = "MOL"
            self._atom3d_isbackbone[idx] = 0
            self._atom3d_occupancy[idx] = 0

        # Loop over bonds
        bond_list = []
        for ibond in self._universe.bonds:
            bond_list.append([ibond[0].index, ibond[1].index])

        # Dimensions
        with open(self._fnamepath, 'r') as f:
            self._unitcell[0, :] = [0.0, 0.0, 0.0]
            self._unitcell[1, :] = [0.0, 0.0, 0.0]
            self._unitcell[2, :] = [0.0, 0.0, 0.0]
            while True:
                line = f.readline()
                if not line:
                    break
                if line.find("@<TRIPOS>CRYSIN") != -1:
                    d = f.readline().split()
                    self._unitcell[0, 0] = d[0]
                    self._unitcell[1, 1] = d[1]
                    self._unitcell[2, 2] = d[2]
                    self._boxlength[0] = d[0]
                    self._boxlength[1] = d[1]
                    self._boxlength[2] = d[2]
                    self._boxangle[0] = d[3]
                    self._boxangle[1] = d[4]
                    self._boxangle[2] = d[5]
        # Topology
        self._topology = top.Topology(natoms=self._natoms, listbonds=bond_list)
        self._nmols = len(self._topology.get_nmols())
        self._mol_residue_list = list(self._universe.residues.resnames)

        # Get atom properties
        imol_idx = 0
        for imol in self._topology.get_nmols():
            for idx in imol:
                self._atom3d_imol[idx] = imol_idx
                self._atom3d_charge[idx] = round(self._universe.atoms.charges[idx], 4)
                self._atom3d_element[idx] = self._universe.atoms.elements[idx]
                self._atom3d_mass[idx] = round(self._universe.atoms.masses[idx], 2)
                self._atom3d_residue[idx] = self._universe.atoms.resids[idx]

            imol_idx += 1

        # List
        charge_list = []
        elements_list = []
        mass_list = []
        isbackbone_list = []
        for idx in range(self._natoms):
            charge_list.append(self._atom3d_charge[idx])
            elements_list.append(self._atom3d_element[idx])
            mass_list.append(self._atom3d_mass[idx])
            isbackbone_list.append(0)

        self._topology.set_elements(elements_list)
        self._topology.set_charge_mol(charge_list)
        self._topology.set_isbackbone(isbackbone_list)
        self._topology.set_mass(mass_list)

        # Assign bond orders ==================================================
        if self._assign_bo:
            self._topology.assign_bond_orders()

        return None

    # *************************************************************************
    def __remove_comments(self):

        """
        This function removes all lines starting with "#"
        from the mol2 file in self.fnamepath.
        The self.fnamepath file is overwritten without comments.

        :return:
        """

        lines = []
        with open(self._fnamepath, 'r') as f:
            while True:
                line = f.readline()
                if not line: break
                if line.startswith('#'):
                    continue
                lines += line

        with open(self._fnamepath, 'w') as f:
            f.writelines(lines)

        return None
