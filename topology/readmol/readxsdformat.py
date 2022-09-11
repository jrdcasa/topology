import xml.etree.ElementTree as etree
from collections import defaultdict
import numpy as np
from periodictable import *
import topology as top
import os
from topology.readmol.readbaseformat import ReadBaseFormat


class ReadXsdFormat(ReadBaseFormat):

    # *************************************************************************
    def __init__(self, filenamepath, assign_bondorders=False):

        super().__init__(filenamepath=filenamepath, assign_bondorders=assign_bondorders)

        self._atom3d_bfactor = defaultdict()
        self._atom3d_occupancy = defaultdict()
        self._heads = list()
        self._tails = list()

        self.read_xsd()

    # *************************************************************************
    def read_xsd(self):

        """
        Base function to read the xml format
        Returns:

        """

        tree = etree.parse(self._fnamepath)
        root = tree.getroot()

        # Find labels in the XML file ==================================================
        atomistic_tree_root_xml = root.find('AtomisticTreeRoot')

        # Setup topology elements ==================================================
        self._isthere_boxdimension = self._setup_box_dimensions(tree)

        self.setup_atoms(tree)
        bondlist = self.setup_bonds(tree)

        self.setup_topology(tree, bondlist)

        molecule_map = self.setup_nmols(tree)
        self.setup_residues(tree, molecule_map)

        # Assign bond orders ==================================================
        if self._assign_bo:
            self._topology.assign_bond_orders()

    # *************************************************************************
    def setup_atoms(self, tree):

        # <Atom3d> label ======================================================
        tmp_atom3d_mapping_xml = tree.findall(".//Atom3d")
        # Transverse the list and remove those Atom elements which contain
        # "ImageOf" attributes. This is something which saves Materials Studio
        # for Supercells
        atom3d_mapping_xml = []
        for item in tmp_atom3d_mapping_xml:
            if item.get("ImageOf"):
                continue
            atom3d_mapping_xml.append(item)

        idx = 0  # Index in the topology
        self._natoms = 0
        for child in atom3d_mapping_xml:
            # index in the XML file (id_xml)
            id_xml = int(child.get("ID"))
            self._atom3d_map[id_xml] = idx
            self._atom3d_idxsd[idx] = int(child.get("ID"))
            try:
                self._atom3d_parent[idx] = int(child.get("Parent"))
            except TypeError:
                self._atom3d_parent[idx] = -1
            self._atom3d_element[idx] = child.get("Components")
            self._atom3d_mass[idx] = formula(self._atom3d_element[idx]).mass
            self._atom3d_occupancy[idx] = 0.0

            # self._atom3d_isbackbone[idx] = 0  --> Backbone
            # self._atom3d_isbackbone[idx] >= 1 --> Branch
            try:
                if child.find("Properties").get("IsBackboneAtom") is not None:
                    if child.find("Properties").get("IsBackboneAtom") == '1':
                        self._atom3d_isbackbone[idx] = 0
                    else:
                        self._atom3d_isbackbone[idx] = 1

                else:
                    self._atom3d_isbackbone[idx] = 1
            except AttributeError:
                self._atom3d_isbackbone[idx] = 0

            try:
                self._atom3d_charge[idx] = float(child.get("Charge"))
            except TypeError:
                self._atom3d_charge[idx] = 0.0

            # from fractional to cartesian coordinates
            if self._isthere_boxdimension and self._nx == 1 and self._ny == 1 and self._nz == 1:
                try:
                    UVW = [float(i) for i in child.get("XYZ").split(",")]  # Fractional coordinates
                except AttributeError:
                    # Supercell sometimes has the following entries:
                    # 	<Atom3d ID="377" Mapping="458" ImageOf="94" Visible="0"/>
                    # These entries must be skipped
                    continue
                omega = np.dot(self._unitcell[0, :], np.cross(self._unitcell[1, :], self._unitcell[2, :]))
                X = self._boxlength[0]*UVW[0] + \
                    self._boxlength[1]*np.cos(self._boxangle[2])*UVW[1] + \
                    self._boxlength[2]*np.cos(self._boxangle[1])*UVW[2]

                Y = self._boxlength[1]*np.sin(self._boxangle[2])*UVW[1] + \
                    self._boxlength[2]*((np.cos(self._boxangle[0]) -
                                         np.cos(self._boxangle[1])*np.cos(self._boxangle[2])) /
                                        np.sin(self._boxangle[1]))*UVW[1]

                Z = omega/(self._boxlength[0]*self._boxlength[1]*np.sin(self._boxangle[2]))*UVW[2]
                XYZ = [X, Y, Z]
                self._atom3d_xyz[idx] = XYZ
            else:
                #for molecules and supercell
                UVW = [float(i) for i in child.get("XYZ").split(",")]
                XYZ = [UVW[0], UVW[1], UVW[2]]
                self._atom3d_xyz[idx] = XYZ
            idx += 1
            self._natoms += 1

        return None
    # *************************************************************************
    def setup_bonds(self, tree):

        # BOND INFORMATION <Bond> label =======================================
        # Get bond properties
        tmp_bond_atom3d_mapping_xml = tree.findall(".//Bond")
        # Transverse the list and remove those Bond elements which contain
        # "ImageOf" attributes. This is something which saves Materials Studio
        # for Supercells
        bond_mapping_xml = []
        for item in tmp_bond_atom3d_mapping_xml:
            if item.get("ImageOf"):
                continue
            bond_mapping_xml.append(item)

        bond_list = []
        self._nbonds = 0
        for child in bond_mapping_xml:
            try:
                iat_id_xml, jat_id_xml = [int(i) for i in child.get("Connects").split(",")]
            except AttributeError:
            # Supercell sometimes has the following entries:
            # 	<Atom3d ID="377" Mapping="458" ImageOf="94" Visible="0"/>
            # These entries must be skipped
                continue

            iat_idx = self._atom3d_map[iat_id_xml]
            jat_idx = self._atom3d_map[jat_id_xml]
            if iat_idx > jat_idx:
                bond_list.append((jat_idx, iat_idx))
            else:
                bond_list.append((iat_idx, jat_idx))
            self._nbonds += 1

        return bond_list

    # *************************************************************************
    def setup_nmols(self, tree):

        tmp_molecule_mapping_xml = tree.findall(".//Molecule")
        # Transverse the list and remove those Atom elements which contain
        # "ImageOf" attributes. This is something which saves Materials Studio
        # for Supercells
        molecule_mapping_xml = []
        for item in tmp_molecule_mapping_xml:
            if item.get("ImageOf"):
                continue
            molecule_mapping_xml.append(item)

        molecule_map = defaultdict()
        imol = 0
        if len(molecule_mapping_xml) !=0 :
            for child in molecule_mapping_xml:
                id_xml = int(child.get("ID"))
                molecule_map[id_xml] = imol
                imol += 1
            self._nmols = imol
        # There is not info about molecules in the xsd
        else:
            self._nmols = len(self._topology.get_nmols())
            for imol in range(0, self._nmols):
                molecule_map[imol] = imol

        nmols = self._topology.get_nmols()
        idx_mol = 0
        for imol in nmols:
            for iatom in imol:
                self._atom3d_imol[iatom] = idx_mol
                self._atom3d_residue[iatom] = idx_mol + 1
            idx_mol += 1

        return molecule_map

    # *************************************************************************
    def setup_residues(self, tree, moleculemap):

        tmp_repeatunit_mapping_xml = tree.findall(".//RepeatUnit")
        # Transverse the list and remove those RepeatUnit elements which contain
        # "ImageOf" attributes. This is something which saves Materials Studio
        # for Supercells
        repeatunit_mapping_xml = []
        for item in tmp_repeatunit_mapping_xml:
            if item.get("ImageOf"):
                continue
            repeatunit_mapping_xml.append(item)

        # <Residue> label ====================================================
        residue_map = defaultdict()
        ires = 0
        # Number of residues equals to repeat units
        if len(repeatunit_mapping_xml) != 0:
            for child in repeatunit_mapping_xml:
                id_xml = int(child.get("ID"))
                residue_map[id_xml] = ires
                self._mol_residue_list.append(child.get("Name"))
                ires += 1
        # Otherwise nunber of residues = nmols
        else:
            ires = len(moleculemap)
            id_res = 0
            for _ in range(ires):
                self._mol_residue_list.append("UNK")
                residue_map[id_res] = id_res

        self._nres = ires

        # <Atom3d> label ======================================================
        tmp_atom3d_mapping_xml = tree.findall(".//Atom3d")
        # Transverse the list and remove those Atom elements which contain
        # "ImageOf" attributes. This is something which saves Materials Studio
        # for Supercells
        atom3d_mapping_xml = []
        for item in tmp_atom3d_mapping_xml:
            if item.get("ImageOf"):
                continue
            atom3d_mapping_xml.append(item)

        idx = 0  # Index in the topology
        ires = 0
        for child in atom3d_mapping_xml:
            # index of the parent in the XML file (id_parent)
            try:
                id_parent = int(child.get("Parent"))
                ires = residue_map[id_parent]
                self._atom3d_molname[idx] = str("{0:03d}".format(ires+1))
            except (KeyError, TypeError):
                self._atom3d_molname[idx] = str("{0:03d}".format(self._atom3d_imol[idx]+1))
            idx += 1  # Index in the topology

    # *************************************************************************
    def setup_topology(self, tree, bondlist=None):

        # Setup topology ======================================================
        self._topology = top.Topology(natoms=self._natoms, listbonds=bondlist)

        elements_list = []
        charge_list = []
        mass_list = []
        isbackbone_list = []
        for idx in self._atom3d_idxsd:
            elements_list.append(self._atom3d_element[idx])
            charge_list.append(self._atom3d_charge[idx])
            mass_list.append(self._atom3d_mass[idx])
            isbackbone_list.append(self._atom3d_isbackbone[idx])
        self._topology.set_elements(elements_list)
        self._topology.set_charge_mol(charge_list)
        self._topology.set_mass(mass_list)
        self._topology.set_isbackbone(isbackbone_list)

    # *************************************************************************
    def _setup_box_dimensions(self, tree):

        # Get the root and main element (AtomisticTreeRoot)
        root = tree.getroot()
        atomistic_tree_root_xml = root.find('AtomisticTreeRoot')

        symmetry_system_xml = atomistic_tree_root_xml.find("SymmetrySystem")
        # No info for box dimensions
        if symmetry_system_xml is None:
            return False

        mapping_set_xml = symmetry_system_xml.find("MappingSet")
        mapping_family_xml = mapping_set_xml.find("MappingFamily")
        identity_mapping_xml = mapping_family_xml.find("IdentityMapping")

        # Simulation Box --> Calculate the dimensions and angles of the simulation box
        for child in identity_mapping_xml:
            if child.tag == "SpaceGroup":
                ll = child.get("DisplayRange").split(",")
                self._nx = int(ll[1]) - int(ll[0])
                self._ny = int(ll[3]) - int(ll[2])
                self._nz = int(ll[5]) - int(ll[4])
                self._unitcell[0, :] = [float(i) for i in child.get("AVector").split(",")]
                self._unitcell[1, :] = [float(i) for i in child.get("BVector").split(",")]
                self._unitcell[2, :] = [float(i) for i in child.get("CVector").split(",")]
                self._boxlength[0] = np.linalg.norm(self._unitcell[0, :], ord=2)
                self._boxlength[1] = np.linalg.norm(self._unitcell[1, :], ord=2)
                self._boxlength[2] = np.linalg.norm(self._unitcell[2, :], ord=2)
                ab = np.dot(self._unitcell[0, :], self._unitcell[1, :])
                ac = np.dot(self._unitcell[0, :], self._unitcell[2, :])
                bc = np.dot(self._unitcell[1, :], self._unitcell[2, :])
                self._boxangle[0] = np.arccos(bc/(self._boxlength[1]*self._boxlength[2]))
                self._boxangle[1] = np.arccos(ac/(self._boxlength[0]*self._boxlength[2]))
                self._boxangle[2] = np.arccos(ab/(self._boxlength[0]*self._boxlength[1]))

        return True
