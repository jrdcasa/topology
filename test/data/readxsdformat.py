import xml.etree.ElementTree as etree
from collections import defaultdict
import sys
import numpy as np
from periodictable import *
import topology as top

class ReadXsdFormat(object):

    __slots__ = ['_fnamepath', '_fname', '_nmols', '_mol_residue_list',
                 '_atomgroup_list', '_units', '_unitcell', '_boxlength',
                 '_boxangle', '_nres', '_map_id_to_idx_dict', '_natoms',
                 '_topology', '_assign_bo', '_namemols_dict']

    # *************************************************************************
    def __init__(self, filenamepath, assign_bondorders=False):

        """
        Read XSD format files creating some useful information for the
        MolecularGraph

        Args:
            filenamepath (str): Path and name of the file for the xsd file

        Members:
            | ``self._filenamepath`` -- Path and name of the file for the xsd file
            | ``self._file`` -- Name of the xsd file
            | ``self._nmols`` -- Number of molecules in the system
            | ``self._nres`` -- Number of residues in the system
            | ``self._natoms`` -- Number of atoms in the system
            | ``self._molreslist`` -- "molreslist  : [ [mol1 name residues] [mol2 name residues], ...]
                                    for molecular graphs
            | ``self._atomgrouplist`` -- "molreslist  :
                     self._atomgrouplist[index_atom (int)]:
                                          [
                                              Atom index,
                                              Local index in the residue,
                                              Name of the atom in the system,
                                              Type of atom in the FF,
                                              Mass atom,
                                              Atom charge,
                                              imol number,
                                              Global number of the residue,
                                              Local number of the residue,
                                              Is or not backbone atom,
                                          ]
            | ``self.map_idx_to_id_dict`` -- dict  : Map the index atomic number to the ID label in the xsd file.
                This map can be used to get the bonds

            | ``self._unitcell`` -- array [3,3] : Contains the box vectors in angstroms
                ([[ a_xx a_xy a_xz ]
                  [ 0.0  b_yy b_yz ]
                  [ 0.0  0.0  c_zz ]]
            | ``self._boxlengths`` -- array [3] : Lengths of the simulation box in angstroms
            | ``self._boxangles`` -- array [3] : Angles of the simulation box in radians
            | ``self._topology`` -- Topology : Topology of the system

        """

        self._fnamepath = filenamepath
        self._fname = filenamepath.split("/")[-1]

        self._nmols = 0
        self._nres = 0
        self._natoms = 0
        self._mol_residue_list = []
        self._atomgroup_list = []
        self._map_id_to_idx_dict = defaultdict()
        self._namemols_dict = defaultdict(list)

        self._units = {'time': None, 'length': 'ang'}

        self._unitcell = np.zeros([3,3],dtype=np.float32)
        self._boxlength = np.zeros(3 ,dtype=np.float32)
        self._boxangle = np.zeros(3, dtype=np.float32)

        self._assign_bo = assign_bondorders

        self.read_xsd()

    # *************************************************************************
    def read_xsd(self):

        """
        Base function to read the xml format
        Returns:

        """

        tree = etree.parse(self._fnamepath)
        root = tree.getroot()

        atomistic_tree_root_xml = root.find('AtomisticTreeRoot')

        # CHANGE: self.molrestlist and self.nmols
        if atomistic_tree_root_xml.find("SymmetrySystem"):
            self.__setup_box_dimensions(tree)
            self.__read_xsd_periodic(tree)
        else:
            self.__read_xsd_nonperiodic(tree)

    # *************************************************************************
    def __setup_box_dimensions(self, tree):

        # Get the root and main element (AtomisticTreeRoot)
        root = tree.getroot()
        atomistic_tree_root_xml = root.find('AtomisticTreeRoot')

        symmetry_system_xml = atomistic_tree_root_xml.find("SymmetrySystem")
        mapping_set_xml = symmetry_system_xml.find("MappingSet")
        mapping_family_xml = mapping_set_xml.find("MappingFamily")
        identity_mapping_xml = mapping_family_xml.find("IdentityMapping")

        # Simulation Box --> Calculate the dimensions and angles of the simulation box
        for child in identity_mapping_xml:
            if child.tag == "SpaceGroup":
                self._unitcell[0,:] = [float(i) for i in child.get("AVector").split(",")]
                self._unitcell[1,:] = [float(i) for i in child.get("BVector").split(",")]
                self._unitcell[2,:] = [float(i) for i in child.get("CVector").split(",")]
                self._boxlength[0] = np.linalg.norm(self._unitcell[0,:], ord=2)
                self._boxlength[1] = np.linalg.norm(self._unitcell[1,:], ord=2)
                self._boxlength[2] = np.linalg.norm(self._unitcell[2,:], ord=2)
                ab = np.dot(self._unitcell[0,:], self._unitcell[1,:])
                ac = np.dot(self._unitcell[0,:], self._unitcell[2,:])
                bc = np.dot(self._unitcell[1,:], self._unitcell[2,:])
                self._boxangle[0] = np.arccos(bc/(self._boxlength[1]*self._boxlength[2]))
                self._boxangle[1] = np.arccos(ac/(self._boxlength[0]*self._boxlength[2]))
                self._boxangle[2] = np.arccos(ab/(self._boxlength[0]*self._boxlength[1]))

    # *************************************************************************
    def __read_xsd_periodic(self, tree):

        # Get the root and main element (AtomisticTreeRoot)
        root = tree.getroot()
        atomistic_tree_root_xml = root.find('AtomisticTreeRoot')

        symmetry_system_xml = atomistic_tree_root_xml.find("SymmetrySystem")
        mapping_set_xml = symmetry_system_xml.find("MappingSet")
        mapping_family_xml = mapping_set_xml.find("MappingFamily")
        identity_mapping_xml = mapping_family_xml.find("IdentityMapping")

        # ATOMIC INFORMATION
        # Create moltypedict -->   {<key ID molecule>: <value type from 0...nmoltypes-1>}
        # Create molresdict -->    {<key ID Linear Chain> : [ <name of the repeat unit>, <ID repeat unit>]
        # Create molsmalldict -->  {<key ID Linear Chain> : <key ID molecule>}
        mol_atomtype_dict = defaultdict(list)
        mol_residue_dict = defaultdict(list)
        mol_small_dict = defaultdict()
        itypemol = 0

        ispolymer = False
        for child in identity_mapping_xml:
            if child.tag == "RepeatUnit":
                ispolymer = True
                break

        if ispolymer:
            imol = 0
            for child in identity_mapping_xml:

                if child.tag == "Molecule":
                    name = child.get("Name")
                    mol_atomtype_dict[name].append(int(itypemol))
                    self._nmols += 1
                    itypemol += 1

                elif child.tag == "RepeatUnit":
                    p = int(child.get('Parent'))
                    n = child.get('Name')
                    ires = int(child.get("ID"))
                    mol_residue_dict[p].append([n,ires])
                    self._nres += 1

                elif child.tag == "LinearChain":
                    p = int(child.get('Parent'))
                    i = int(child.get("ID"))
                    namepol = child.get("Name")
                    self._namemols_dict[namepol].append(imol)
                    imol += 1
                    mol_small_dict[i] = p

                else:
                    continue
        else:
            imol = 0
            for child in identity_mapping_xml:

                if child.tag == "Molecule":
                    namemol = child.get("Name")
                    ires = int(child.get("ID"))
                    mol_atomtype_dict[namemol].append(int(itypemol))
                    mol_residue_dict[imol].append([namemol,ires])
                    self._nmols += 1
                    itypemol += 1
                    self._nres += 1
                    self._namemols_dict[namemol].append(imol)
                    imol += 1

        # Map ID to imol
        imol = 0
        map_imol = defaultdict()
        for key in sorted(mol_residue_dict):
            map_imol[key] = imol
            self._mol_residue_list.append([])
            imol += 1
        # self.molreslist --> Setup
        for key in sorted(mol_residue_dict):
            imol = map_imol[key]
            for item in mol_residue_dict[key]:
                self._mol_residue_list[imol].append(item[0])

        # Map ID to ires
        ires = 0
        map_ires = defaultdict(list)
        for item1 in sorted(mol_residue_dict):
            ireslocal = 0
            for item2 in sorted(mol_residue_dict[item1]):
                key = item2[1]
                imol = map_imol[item1]
                map_ires[key] =[ires, ireslocal, imol]
                ires += 1
                ireslocal += 1

        # Get atomic properties --> idatoms, coordinates, charges , ....
        imol_c = -1
        idxatom_local = 0
        for child in identity_mapping_xml:
            if child.tag == "Atom3d":
                idms_atoms = int(child.get("ID"))
                idxatom = int(child.get("UserID"))-1
                imol = map_ires[int(child.get("Parent"))][2]
                if imol_c == -1 or imol_c != imol:
                    imol_c = imol
                    idxatom_local = 0
                else:
                    idxatom_local += 1
                element = child.get("Components")
                nameatom = child.get("Name")
                typeatom = child.get("ForcefieldType")
                ires = map_ires[int(child.get("Parent"))][0]
                ireslocal = map_ires[int(child.get("Parent"))][1]
                massatom = formula(element).mass
                if child.find("Properties").get("IsBackboneAtom") is not None:
                    isbackbone = int(child.find("Properties").get("IsBackboneAtom"))
                else:
                    isbackbone = 0

                # from fractional to cartesian coordinates
                UVW = [float(i) for i in child.get("XYZ").split(",")]  # Fractional coordinates
                omega = np.dot(self._unitcell[0,:], np.cross(self._unitcell[1,:], self._unitcell[2,:]))
                XYZ = [self._boxlength[0]*UVW[0] +
                       self._boxlength[1]*np.cos(self._boxangle[2])*UVW[1] +
                       self._boxlength[2]*np.cos(self._boxangle[1])*UVW[2],
                       self._boxlength[1]*np.sin(self._boxangle[2])*UVW[1] +
                       self._boxlength[2]*((np.cos(self._boxangle[0]) -
                                            np.cos(self._boxangle[1])*np.cos(self._boxangle[2]))/
                                           np.sin(self._boxangle[1]))*UVW[1],
                       omega/(self._boxlength[0]*self._boxlength[1]*np.sin(self._boxangle[2]))*UVW[2]]

                try:
                    charge_list = float(child.get("Charge"))
                except TypeError:
                    charge_list = 0.0

                self._atomgroup_list.append([idxatom, idxatom_local, element, nameatom,
                                             typeatom, massatom, charge_list, imol, ires,
                                             ireslocal, isbackbone, XYZ])
                self._map_id_to_idx_dict[idms_atoms] = idxatom

        self._natoms = len(self._atomgroup_list)

        # BOND INFORMATION
        # Get bond properties
        bond_list = []
        for child in identity_mapping_xml:
            if child.tag == "Bond":
                iat_id, jat_id = [int(i) for i in child.get("Connects").split(",")]
                iat_idx = self._map_id_to_idx_dict[iat_id]
                jat_idx = self._map_id_to_idx_dict[jat_id]
                if iat_idx > jat_idx:
                    bond_list.append((jat_idx, iat_idx))
                else:
                    bond_list.append((iat_idx, jat_idx))

        # Setup topology
        self._topology = top.Topology(natoms=self._natoms, listbonds=bond_list)

        # Fill up elements, charge, masses, rings,
        elements_list = []
        charge_list = []
        mass_list = []
        isbackbone_list = []
        type_list = []
        name_list = []
        for item in self._atomgroup_list:
            elements_list.append(item[2])
            charge_list.append(item[6])
            mass_list.append(item[5])
            isbackbone_list.append(item[10])
            type_list.append(item[4])
            name_list.append(item[3])
        self._topology.set_elements(elements_list)
        self._topology.set_charge_mol(charge_list)
        self._topology.set_mass(mass_list)
        self._topology.set_isbackbone(isbackbone_list)
        self._topology.set_type(type_list)
        self._topology.set_name(name_list)
        self._topology.guess_nringsCauchy()
        self._topology._topologyfile = self._fname

        # Assign bond orders
        if self._assign_bo:
            self._topology.assign_bond_orders()

        return None

    # *************************************************************************
    def __read_xsd_nonperiodic(self, tree):

        # Get the root and main element (AtomisticTreeRoot)
        root = tree.getroot()
        atomistic_tree_root_xml = root.find('AtomisticTreeRoot')
        molecule_xml = atomistic_tree_root_xml.find("Molecule")

        if molecule_xml is None:
            print("No hierarchy")
            print("Need to be implemented")
        else:
            print("hierarchy")
            self.__read_xsd_nonperiodic_with_hierarchy(molecule_xml)




    # *************************************************************************
    def __read_xsd_nonperiodic_with_hierarchy(self, molecule_xml):

        # Atom3d is the first molecule
        atom_second_xml = molecule_xml.findall("Atom3d")
        idxatom = -1
        idxatom_local= -1

        imol = 0
        ires = 0
        ireslocal = 0
        self._map_id_to_idx_dict = defaultdict()
        for item in atom_second_xml:
            idms_atoms = int(item.get("ID"))
            idxatom += 1
            self._map_id_to_idx_dict[idms_atoms] = idxatom
            idxatom_local += 1
            element = item.get("Components")
            nameatom = item.get("Name")
            typeatom = item.get("ForcefieldType")
            try:
                charge_list = float(item.get("Charge"))
            except TypeError:
                charge_list = 0.0
            massatom = formula(element).mass

            isbackbone = 0
            if item.find("Properties") is not None:
                if item.find("Properties").get("IsBackboneAtom") is not None:
                    isbackbone = int(item.find("Properties").get("IsBackboneAtom"))

            XYZ = [float(i) for i in item.get("XYZ").split(",")]  # Cartesian coordinates
            self._atomgroup_list.append([idxatom, idxatom_local, element, nameatom,
                                         typeatom, massatom, charge_list, imol, ires,
                                         ireslocal, isbackbone, XYZ])

        # The rest of molecules are under "Molecule" label
        molecule_second_xml = molecule_xml.findall("Molecule")
        for item in molecule_second_xml:
            imol += 1
            ires += 1
            ireslocal = 0
            for child in item:
                if child.tag == "Atom3d":
                    idms_atoms = int(child.get("ID"))
                    idxatom += 1
                    self._map_id_to_idx_dict[idms_atoms] = idxatom
                    idxatom_local += 1
                    element = child.get("Components")
                    nameatom = child.get("Name")
                    typeatom = child.get("ForcefieldType")
                    try:
                        charge_list = float(child.get("Charge"))
                    except TypeError:
                        charge_list = 0.0
                    massatom = formula(element).mass

                    isbackbone = 0
                    if child.find("Properties") is not None:
                        if child.find("Properties").get("IsBackboneAtom") is not None:
                            isbackbone = int(child.find("Properties").get("IsBackboneAtom"))

                    XYZ = [float(i) for i in child.get("XYZ").split(",")]  # Cartesian coordinates
                    self._atomgroup_list.append([idxatom, idxatom_local, element, nameatom,
                                                 typeatom, massatom, charge_list, imol, ires,
                                                 ireslocal, isbackbone, XYZ])

        self._natoms = len(self._atomgroup_list)

        # BOND INFORMATION
        # Get bond properties
        bond_list = []
        bond_second_xml = molecule_xml.findall("Bond")
        for item in bond_second_xml:
            iat_id, jat_id = [int(i) for i in item.get("Connects").split(",")]
            iat_idx = self._map_id_to_idx_dict[iat_id]
            jat_idx = self._map_id_to_idx_dict[jat_id]
            if iat_idx > jat_idx:
                bond_list.append((jat_idx, iat_idx))
            else:
                bond_list.append((iat_idx, jat_idx))

        molecule_second_xml = molecule_xml.findall("Molecule")
        for item in molecule_second_xml:
            for child in item:
                if child.tag == "Bond":
                    iat_id, jat_id = [int(i) for i in child.get("Connects").split(",")]
                    iat_idx = self._map_id_to_idx_dict[iat_id]
                    jat_idx = self._map_id_to_idx_dict[jat_id]
                    if iat_idx > jat_idx:
                        bond_list.append((jat_idx, iat_idx))
                    else:
                        bond_list.append((iat_idx, jat_idx))
        # Setup topology
        self._topology = top.Topology(natoms=self._natoms, listbonds=bond_list)

        # Fill up elements, charge, masses, rings,
        elements_list = []
        charge_list = []
        mass_list = []
        isbackbone_list = []
        type_list = []
        name_list = []
        for item in self._atomgroup_list:
            elements_list.append(item[2])
            charge_list.append(item[6])
            mass_list.append(item[5])
            isbackbone_list.append(item[10])
            type_list.append(item[4])
            name_list.append(item[3])
        self._topology.set_elements(elements_list)
        self._topology.set_charge_mol(charge_list)
        self._topology.set_mass(mass_list)
        self._topology.set_isbackbone(isbackbone_list)
        self._topology.set_type(type_list)
        self._topology.set_name(name_list)
        self._topology.guess_nringsCauchy()
        self._topology._topologyfile = self._fname

        # Assign names to the molecule
        imol = 0
        atom_second_xml = molecule_xml.findall("Atom3d")
        namemol = molecule_xml.get("Name")
        self._namemols_dict[namemol].append(imol)

        molecule_second_xml = molecule_xml.findall("Molecule")
        for item in molecule_second_xml:
            imol += 1
            namemol = item.get("Name")
            self._namemols_dict[namemol].append(imol)

        # Assign bond orders
        if self._assign_bo:
            self._topology.assign_bond_orders()

        return None

    # *************************************************************************
    def get_natoms(self):

        return self._natoms

    # *************************************************************************
    def get_nmols(self):

        return self._nmols

    # *************************************************************************
    def get_molreslist(self):

        return self._mol_residue_list

    # *************************************************************************
    def get_atomgrouplist(self):

        return self._atomgroup_list

    # *************************************************************************
    def write_xyz(self, filename_xyz="test.xyz"):

        """
        Write a xyz file to check the structure
        """

        with open(filename_xyz, 'w') as fxyz:

            fxyz.write("{}\n".format(int(self.get_natoms())))
            fxyz.write("MOL\n")

            for item in self._atomgroup_list:
                line = "{0:s} {1:.5f} {2:.5f} {3:.5f}\n".format(item[2],
                                                                item[-1][0],
                                                                item[-1][1],
                                                                item[-1][2])
                fxyz.write(line)

    # *************************************************************************
    def write_gro(self, filename_gro="test.gro"):

        """
        Write a xyz file to check the structure
        """

        fmt = {
            'n_atoms': "{0:5d}\n",  # number of atoms
            # coordinates output format, see http://chembytes.wikidot.com/g-grofile
            'xyz': "{resid:>5d}{resname:<5.5s}{name:>5.5s}{index:>5d}{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}\n",
            # unitcell
            'box_orthorhombic': "{0:10.5f} {1:9.5f} {2:9.5f}\n",
            'box_triclinic': "{0:10.5f} {1:9.5f} {2:9.5f} "
                             "{3:9.5f} {4:9.5f} {5:9.5f} "
                             "{6:9.5f} {7:9.5f} {8:9.5f}\n"
        }

        with open(filename_gro, 'w') as fgro:

            fgro.write('Written by MStoLAGRO (J.Ramos)\n')
            fgro.write(fmt['n_atoms'].format(int(self.get_natoms())))

            for item in self._atomgroup_list:
                line = "{0:s} {1:.5f} {2:.5f} {3:.5f}\n".format(item[2],
                                                                item[-1][0],
                                                                item[-1][1],
                                                                item[-1][2])

                imol = int(item[7])
                dd = [k for k, v in self._namemols_dict.items() if imol in v]

                fgro.write(fmt['xyz'].format(
                    resid=item[7]+1,
                    resname="{}".format(dd[0][0:3]),
                    index = item[0]+1,
                    name=item[2],
                    pos=[0.1*i for i in item[-1]]
                ))

            # self._unitcell =                  GROMACS box vectors
            #       ([[ a_xx a_xy a_xz ]  Transpose   ([[ v1(x) 0.0   0.0   ]
            #         [ 0.0  b_yy b_yz ]   =======>     [ v2(x) v2(y) 0.0   ]
            #         [ 0.0  0.0  c_zz ]]               [ v3(x) v3(y) v3[z] ]]
            #
            #    v1(x) --> [0,0]; v2(y)--> [1,1]; v3(z)--> [2,2]
            #    v1(y) --> [1,0]; v1(z)--> [2,0]; v2(x)--> [0,1]
            #    v2(z) --> [2,1]; v3(x)--> [0,3]; v3(z)--> [1,2]
            #   Gromacs v1(y) = v1(z) = v2(z) = 0.0
            if self._unitcell[0,1] != 0.0 or self._unitcell[0,2] or self._unitcell[1,2]:
                fgro.write(fmt['box_triclinic'].format(self._unitcell[0,0]*0.1,
                                                       self._unitcell[1,1]*0.1,
                                                       self._unitcell[2,2]*0.1,
                                                       self._unitcell[1,0]*0.1,
                                                       self._unitcell[2,0]*0.1,
                                                       self._unitcell[0,1]*0.1,
                                                       self._unitcell[2,1]*0.1,
                                                       self._unitcell[0,3]*0.1,
                                                       self._unitcell[1,2])*0.1)
            else:
                fgro.write(fmt['box_orthorhombic'].format(self._unitcell[0,0]*0.1,
                                                          self._unitcell[1,1]*0.1,
                                                          self._unitcell[2,2]*0.1))
