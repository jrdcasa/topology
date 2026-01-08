import re
import json

from topology.BasicFormat import BasicFormat
from topology.chirality_from_coordinates import chirality_from_coords
from topology.force_fields.ForceField import ForceField
from topology.atomic_data import element_cov_radius


class RaspaFormat(BasicFormat):

    """
    RaspaFormat Class
    """

    # ====================================================================
    def __init__(self, mol_mda, logger=None):

        # Parent class
        super().__init__(mol_mda, logger)

        self._tc = 0.000  # Critical temperature
        self._pc = 0.000  # Critical pressure
        self._omega = 0.000  # Acentric factor

        self._graph_mol, self._graph_mol_components, self._graph_unique_mol, \
            self._graph_component_nodes, self._nmols, \
            self._nmols_unique = self._number_molecules_graph()

        self._molecular_names()

    # ====================================================================
    def _setup_bonds(self, igraph, offset=0):

        """

        Args:
            igraph:
            offset:

        Returns:

        """

        nbonds_corrected = 0
        lines_bond = "# Bond stretch: atom i-j, type, parameters\n"
        list_bonds_json = []
        list_connectivity = []
        for ibond in self._mol_mda.bonds:
            iat = ibond.indices[1]
            jat = ibond.indices[0]
            if iat not in igraph.nodes():
                continue
            itype = ibond.atoms[0].name
            jtype = ibond.atoms[1].name
            # Do not write bond for Na, Li, ...
            if itype.upper() in self._ions_list or jtype.upper() in self._ions_list:
                continue
            blabel_d = itype + "-" + jtype
            blabel_i = jtype + "-" + itype
            if blabel_d not in self._ff._bondtypes_dict and \
                    blabel_i not in self._ff._bondtypes_dict:
                m = "\t\t Bond type {0:s} or {1:s} was not found in the force field parameters\n".\
                     format(blabel_d, blabel_i)
                m1 = "\t\t" + len(m) * "*" + "\n"
                m += "\t\t Suggestion: Insert parameters bond_dict_types in the {0:s} force field\n".\
                     format(self._ff._ffname)
                m += "\t\t Changes should be done in ForceField.py class\n"
                print("\n" + m1 + m + m1) if self._logger is None else self._logger.error("\n"+m1+m+m1)
                exit()
            if blabel_i in self._ff._bondtypes_dict:
                blabel_d = blabel_i
                tmp = iat
                iat = jat
                jat = iat
            if self._ff._bondtypes_dict[blabel_d][0] == 1:
                type_bond_potential = "HARMONIC_BOND"
                type_bond_potential_raspa3 = "HARMONIC"
            lines_bond += "{0:>10d} {1:>10d} {2:>15s} {3:.1f} {4:14.3f}\n". \
                format(iat - offset, jat - offset,
                       type_bond_potential,
                       self._ff._bondtypes_dict[blabel_d][2] / (100 * self._kb_kjmol),  #from kj/molnm2 to K/A2
                       self._ff._bondtypes_dict[blabel_d][1] * 10)

            itypeff = self._index_to_type_dict[iat]
            jtypeff = self._index_to_type_dict[jat]

            list_bonds_json.append([
                ["{0:s}".format(itypeff), "{0:s}".format(jtypeff)],
                "{0:s}".format(type_bond_potential_raspa3),
                [round(self._ff._bondtypes_dict[blabel_d][2] / (100 * self._kb_kjmol), 1),
                 round(self._ff._bondtypes_dict[blabel_d][1] * 10, 3)]
                                   ])
            list_connectivity.append([round(iat-offset,0), round(jat-offset,0)])

            nbonds_corrected += 1

        unique_bonds = []
        for b in list_bonds_json:
            if b not in unique_bonds:
                unique_bonds.append(b)

        list_bonds_json = unique_bonds
        return lines_bond, nbonds_corrected, list_bonds_json, list_connectivity

    # ====================================================================
    def _setup_angles(self, igraph, offset=0):

        """

        Args:
            igraph:
            offset:

        Returns:

        """

        lines_angle = "# Bond Bending: atom i-j-k, type, parameters\n"
        nangles_corrected = 0
        list_bends_json = []
        # angles = sorted(list(Mda.topology.guessers.guess_angles(self._mol_mda.bonds)))
        angles = self._mol_mda.angles
        # self._mol_mda.add_TopologyAttr('angles', angles)
        for iangle in self._mol_mda.angles:
            iat = iangle.indices[0]
            jat = iangle.indices[1]
            kat = iangle.indices[2]
            itype = iangle.atoms[0].name
            jtype = iangle.atoms[1].name
            ktype = iangle.atoms[2].name
            if iat not in igraph.nodes():
                continue
            # Do not write angles for Na, Li, ...
            if itype.upper() in self._ions_list or \
                jtype.upper() in self._ions_list or \
                ktype.upper() in self._ions_list:
                continue

            alabel_d = itype+"-"+jtype+"-"+ktype
            alabel_i = ktype+"-"+jtype+"-"+itype
            if alabel_d not in self._ff._angletypes_dict and \
               alabel_i not in self._ff._angletypes_dict:
                m = "\t\t Angle type {0:s} or {1:s} was not found in the force field parameters\n".\
                     format(alabel_d, alabel_i)
                m1 = "\t\t"+len(m)*"*"+"\n"
                m += "\t\t Suggestion: Insert parameters angle_dict_types in the {0:s} force field\n".\
                     format(self._ff._ffname)
                m += "\t\t Changes should be done in ForceField.py class\n"
                print("\n"+m1+m+m1) if self._logger is None else self._logger.error("\n"+m1+m+m1)
                exit()
            if alabel_i in self._ff._angletypes_dict:
                alabel_d = alabel_i
                tmp = itype
                itype = ktype
                ktype = tmp
                tmp = iat
                iat = kat
                kat = tmp

            nangles_corrected += 1

            if self._ff._angletypes_dict[alabel_d][0] == 1:
                type_bend_potential = "HARMONIC_BEND"
                type_bend_potential_raspa3 = "HARMONIC"
            lines_angle += "{0:>10d} {1:>10d} {2:>10d} {3:>15s} {4:.1f} {5:14.3f}\n". \
                format(iat - offset, jat - offset, kat - offset,
                       type_bend_potential,
                       self._ff._angletypes_dict[alabel_d][2] / self._kb_kjmol,  #from kj/molnm2 to K/A2
                       self._ff._angletypes_dict[alabel_d][1])

            itypeff = self._index_to_type_dict[iat]
            jtypeff = self._index_to_type_dict[jat]
            ktypeff = self._index_to_type_dict[kat]

            list_bends_json.append([
                ["{0:s}".format(itypeff), "{0:s}".format(jtypeff), "{0:s}".format(ktypeff) ],
                "{0:s}".format(type_bend_potential_raspa3),
                [round(self._ff._angletypes_dict[alabel_d][2] / self._kb_kjmol, 1),
                 round(self._ff._angletypes_dict[alabel_d][1], 3)]
            ])

        unique_bends = []
        for b in list_bends_json:
            if b not in unique_bends:
                unique_bends.append(b)

        list_bends_json = unique_bends

        return lines_angle, nangles_corrected, list_bends_json

    # ====================================================================
    def _setup_dihedrals(self, igraph, offset=0):

        """

        Args:
            igraph:
            offset:

        Returns:

        """
        lines_dihedrals = "# Dihedrals: atom i-j-k-l, type, parameters\n"
        ndihedrals_corrected = 0
        list_dihedrals_json = []
        # dihedrals = sorted(list(Mda.topology.guessers.guess_dihedrals(self._mol_mda.angles)))
        # self._mol_mda.add_TopologyAttr('dihedrals', dihedrals)
        # angles = self._mol_mda.angles

        for idih in self._mol_mda.dihedrals:
            iat = idih.indices[0]
            jat = idih.indices[1]
            kat = idih.indices[2]
            lat = idih.indices[3]
            itype = idih.atoms[0].name
            jtype = idih.atoms[1].name
            ktype = idih.atoms[2].name
            ltype = idih.atoms[3].name
            if iat not in igraph.nodes():
                continue
            # Do not write angles for Na, Li, ...
            if itype.upper() in self._ions_list or \
                jtype.upper() in self._ions_list or \
                ktype.upper() in self._ions_list or \
                ltype.upper() in self._ions_list:
                continue

            dlabel_d = itype+"-"+jtype+"-"+ktype+"-"+ltype
            dlabel_i = ltype+"-"+ktype+"-"+jtype+"-"+itype
            if dlabel_d not in self._ff._dihedraltypes_dict and \
               dlabel_i not in self._ff._dihedraltypes_dict:
                m = "\t\t Dihedral type {0:s} or {1:s} was not found in the force field parameters\n".\
                     format(dlabel_d, dlabel_i)
                m1 = "\t\t"+len(m)*"*"+"\n"
                m += "\t\t Suggestion: Insert parameters dihedral_dict_types in the {0:s} force field\n".\
                     format(self._ff._ffname)
                m += "\t\t Changes should be done in ForceField.py class\n"
                print("\n"+m1+m+m1) if self._logger is None else self._logger.error("\n"+m1+m+m1)
                exit()
            if dlabel_i in self._ff._dihedraltypes_dict:
                dlabel_d = dlabel_i
                itype, jtype, ktype, ltype = dlabel_i.split("-")
                str_at = "{0:d}-{1:d}-{2:d}-{3:d}".format(lat, kat, jat, iat)
                iat, jat, kat, lat = [int(i) for i in str_at.split("-")]

            ndihedrals_corrected += 1
            if self._ff._dihedraltypes_dict[dlabel_d][0] == 3:
                type_dihedral_potential = "SIX_COSINE_DIHEDRAL"
                type_dihedral_potential_raspa3 = "SIX_COSINE"
            lines_dihedrals += "{0:>10d} {1:>10d} {2:>10d} {3:>10d} {4:s} ". \
                                format(iat - offset, jat - offset, kat - offset, lat - offset,
                                       type_dihedral_potential)
            lcoeffs = []
            for item in self._ff._dihedraltypes_dict[dlabel_d][1]:
                lines_dihedrals += "{0:10.5f} ".format(item/self._kb_kjmol)
                lcoeffs.append(round(item/self._kb_kjmol, 5))
            lines_dihedrals += "\n"

            itypeff = self._index_to_type_dict[iat]
            jtypeff = self._index_to_type_dict[jat]
            ktypeff = self._index_to_type_dict[kat]
            ltypeff = self._index_to_type_dict[lat]

            list_dihedrals_json.append([
                ["{0:s}".format(itypeff), "{0:s}".format(jtypeff), "{0:s}".format(ktypeff), "{0:s}".format(ltypeff) ],
                "{0:s}".format(type_dihedral_potential_raspa3),
                lcoeffs
            ])

        unique_dihedrals = []
        for b in list_dihedrals_json:
            if b not in unique_dihedrals:
                unique_dihedrals.append(b)

        list_dihedrals_json = unique_dihedrals


        return lines_dihedrals, ndihedrals_corrected, list_dihedrals_json

    # ====================================================================
    def _setup_impropers(self, igraph, offset=0, improper_file=None):

        """

        Args:
            igraph:
            offset:
            improper_file:

        Returns:

        """

        lines_improper = "# Improper torsion atom i-j-k-l\n"
        if improper_file is None:
            improper_list = self._guess_impropers()
            if len(improper_list) != 0:
                m = "\n\t\t Guessed improper angles are found!!!!\n"
                m += "\t\t You must give an improper.ndx with parameters!!!!\n"
                m += "\t\t Use --improper improper.ndx option.\n"
                for item in improper_list:
                    m += "\t\t "
                    m += ", ".join(map(str, item))
                    m += "\n"
                print(m) if self._logger is None else self._logger.error(m)
                exit()
        else:
            with open(improper_file, 'r') as fimp:
                lines = fimp.readlines()
                # nkinds = int(lines[0].split()[1])
                improper_list = []
                nimpropers_list = []
                idx = 1
                while idx < len(lines):
                    if re.match("^#", lines[idx]):
                        idx += 1
                        continue
                    a, b = lines[idx][:-1].split()
                    b = int(b)
                    nimpropers_list.append(b)
                    idx += 1
                    improper_lines = ''
                    for j in range(idx, idx+b):
                        iline = lines[j].split()
                        if int(iline[4]) == 4:
                            type_improper = "CVFF_IMPROPER_DIHEDRAL"
                            lines_improper += "{0:8d} {1:8d} {2:8d} {3:8d}  {4:s}   {5:7.3f} {6:d} {7:6.1f} \n".\
                                              format(int(iline[0]), int(iline[1]), int(iline[2]), int(iline[3]),
                                                     type_improper, float(iline[6])/self._kb_kjmol, int(iline[7]),
                                                     float(iline[5]))
                        elif int(iline[4]) == 2:
                            type_improper = "HARMONIC_IMPROPER_DIHEDRAL"
                            lines_improper += "{0:8d} {1:8d} {2:8d} {3:8d}  {4:s}   {5:7.3f} {6:6.1f}\n".\
                                              format(int(iline[0]), int(iline[1]), int(iline[2]), int(iline[3]),
                                                     type_improper, float(iline[6])/self._kb_kjmol, float(iline[5]) )
                        else:
                            m = "\n\t\tERROR: Improper type {} is not implemented".format(int(iline[4]))
                            print(m) if self._logger is None else self._logger.error(m)
                    improper_list.append(improper_lines)
                    idx += b

        return lines_improper, len(improper_list)

    # ====================================================================
    def _setup_pairs(self, igraph, offset=0):

        """

        Args:
            igraph:
            offset:

        Returns:

        """
        lines_pairs = "# Intra VDW and Coulombic(1-4): i-j\n"

        paths_of_length_3 = self.find_paths_of_length_n(igraph, 4)
        pair_list = []
        for item in paths_of_length_3:
            iat = item[0]
            lat = item[-1]
            if iat > lat:
                if [lat, iat] not in pair_list:
                    pair_list.append([lat, iat])
            else:
                if [iat, lat] not in pair_list:
                    pair_list.append([iat, lat])
        pair_list = sorted(pair_list)

        for ipair in pair_list:
            lines_pairs += "{0:>10d} {1:>10d}\n".format(ipair[0]-offset, ipair[1]-offset)

        return lines_pairs, len(pair_list)

    # ====================================================================
    def write_def_molecules(self, ffname, chargefile, improper_file=None):

        """

        Args:
            ffname:
            chargefile:
            improper_file:

        Returns:

        """

        # Setup force field
        self._ff = ForceField(ffname, log=self._logger)

        self.read_charges_types(chargefile)

        molecules = []
        offset = 0
        for ikey, igraph in self._graph_unique_mol.items():

            G = igraph[0]

            ffname = G.name + '_template.def'
            m = "\t\t Writting {}".format(ffname)
            print(m) if self._logger is None else self._logger.info(m)

            # Section Critical constants =======================================
            iline_def = "# critical constants: Temperature [T], Pressure [Pa], and Acentric factor [omega]\n"
            if self._tc == 0.0 or self._omega == 0.0 or self._pc == 0.0:
                iline_def += "{0:.3f}\n".format(self._tc)
                iline_def += "{0:.3f}\n".format(self._pc)
                iline_def += "{0:.3f}\n".format(self._omega)
                m = "\n\t\t Critical values are not provided for system {}".format(ffname)
                print(m) if self._logger is None else self._logger.warn(m)
            else:
                iline_def += "{0:.2f}\n".format(self._tc)
                iline_def += "{0:.2f}\n".format(self._pc)
                iline_def += "{0:.2f}\n".format(self._omega)

            # Section Number of atoms =======================================
            iline_def += "# Number of atoms\n"
            iline_def += "{}\n".format(len(G))
            iline_def += "# Number of Groups\n"
            iline_def += "{}\n".format(1)
            iline_def += "# {} group\n".format(G.name)
            iline_def += "flexible\n"
            iline_def += "# Number of atoms\n"
            iline_def += "{}\n".format(len(G))
            iline_def += "# Atomic positions\n"
            for iatom_idx in sorted(G):
                iatom = self._mol_mda.atoms[iatom_idx]
                atomtype = self._read_typeatoms[iatom_idx]
                atomname = iatom.name
                # print(iatom_idx, atomname, atomtype)
                iline_def += "{0:<3d} {1:s}\n".format(iatom_idx - offset, atomtype)

            # Topology
            iline_bonds, n_bonds, _, _ = self._setup_bonds(G, offset)
            iline_angles, n_angles, _ = self._setup_angles(G, offset)
            iline_dihedrals, n_dihedrals, _ = self._setup_dihedrals(G, offset)
            iline_impropers, n_impropers = self._setup_impropers(G, offset, improper_file)
            iline_pairs, n_pairs = self._setup_pairs(G, offset)

            iline_def += "# Chiral centers Bond  BondDipoles Bend  UrayBradley InvBend  " \
                         "Torsion Imp. Torsion Bond/Bond Stretch/Bend Bend/Torsion IntraVDW IntraCoulomb\n"

            self._chiral_centers = chirality_from_coords(self._mol_mda, logger=self._logger)

            n_chirals = len(self._chiral_centers)

            iline_def += "{0:2d} ".format(n_chirals)   # Chiral
            iline_def += "{0:3d} ".format(n_bonds)     # Bonds
            iline_def += "{0:3d} ".format(0)           # BondDipoles
            iline_def += "{0:3d} ".format(n_angles)    # Bend
            iline_def += "{0:3d} ".format(0)           # UreyBradley
            iline_def += "{0:3d} ".format(0)           # InvBend
            iline_def += "{0:3d} ".format(n_dihedrals) # Dihedrals
            iline_def += "{0:3d} ".format(n_impropers) # Torsion Impropers
            iline_def += "{0:3d} ".format(0)           # Bond/Bond
            iline_def += "{0:3d} ".format(0)           # Stretch/Bend
            iline_def += "{0:3d} ".format(0)           # Bend/Bend
            iline_def += "{0:3d} ".format(0)           # Stretch/Torsion
            iline_def += "{0:3d} ".format(0)           # Bend/Torsion
            iline_def += "{0:3d} ".format(n_pairs)     # IntraVDW
            iline_def += "{0:3d} ".format(n_pairs)     # IntraCoulomb

            iline_def += "\n"
            offset += len(G.nodes)


            iline_def += iline_bonds
            iline_def += iline_angles
            iline_def += iline_dihedrals
            iline_def += iline_impropers
            iline_def += iline_pairs
            iline_def += "#Number of config moves\n"
            iline_def += "#USER#\n"

            # Write the file
            with (open(ffname, 'w')) as ff:
                ff.writelines(iline_def)

    # ====================================================================
    def write_json_files(self, ffname, chargefile, improper_file=None):

        """

        Args:
            ffmol:
            chargefile:
            improper_file:

        Returns:

        """

        # ===================================================================
        def empty_molecule_json():

            """

            Returns:

            """

            mol_json = {
                    "CriticalTemperature": None,
                    "CriticalPressure": None,
                    "AcentricFactor": None,
                    "Type": "",
                    "pseudoAtoms": [],
                    "Connectivity": [],
                    "Bonds": [],
                    "Bends": [],
                    "Torsions": [],
                    "VanDerWaals": None
            }
            return mol_json

        # ===================================================================
        def empty_forcefield_json():

            """

            Returns:

                PseudoAtoms structure:
                            {
                                "name": None,
                                "framework": False,
                                "print_to_output": False,
                                "element": None,
                                "print_as": None,
                                "mass": 0.0,
                                "charge": 0.0,
                                "source": ""
                            }
                "SelfInteractions": [
                            {
                                "name": None,
                                "type": None,
                                "parameters": [],
                                "source": ""
                            }


            """



            ff_json = {
                    "MixingRule": None,
                        "TruncationMethod": None,
                        "TailCorrections": None,
                        "CutOffVDW": None,
                        "PseudoAtoms": [],
                        "SelfInteractions": [],
            }
            return ff_json

        # ===================================================================
        def empty_simulation_json():

            """

            Returns:

            """

            simulation_json = {
              "SimulationType": "MonteCarlo",
              "NumberOfCycles": 10000,
              "NumberOfInitializationCycles": 5000,
              "NumberOfEquilibrationCycles": 5000,
              "PrintEvery": 1000,
              "Systems": [
                {
                  "Type": "Framework",
                  "Name": "MFI_SI",
                  "NumberOfUnitCells": [
                    2,
                    2,
                    2
                  ],
                  "ExternalTemperature": 353.0,
                  "ExternalPressure": 1.0e5,
                  "ChargeMethod": "Ewald",
                  "HybridMCProbability": 0.01,
                  "HybridMCMoveNumberOfSteps": 100,
                  "TimeStep": 0.005,
                  "Ensemble": "NVE"
                }
              ],
              "Components": [
                {
                  "Name": "CO2",
                  "MoleculeDefinition": "ExampleDefinitions",
                  "FugacityCoefficient": 1.0,
                  "ThermodynamicIntegration": True,
                  "TranslationProbability": 0.5,
                  "RotationProbability": 0.5,
                  "CFCMC_CBMC_SwapProbability": 1.0,
                  "CreateNumberOfMolecules": 0
                }
              ]
            }


            return simulation_json


        # ===================================================================
        def to_compact_json(obj, indent = 2, level = 0, inline=False):
            """
            Pretty-print JSON with objects indented, but all arrays kept compact on one line.
            """
            space = " " * indent * level
            next_space = " " * indent * (level + 1)

            if isinstance(obj, dict):
                items = []
                for k, v in obj.items():
                    items.append(f'{next_space}"{k}": {to_compact_json(v, indent, level + 1)}')
                return "{\n" + ",\n".join(items) + "\n" + space + "}"

            elif isinstance(obj, list):
                if inline:
                    # keep this list compact on one line
                    inner = ", ".join(to_compact_json(v, indent, 0, inline=True) for v in obj)
                    return f"[ {inner} ]"
                else:
                    # pretty-print top-level lists, but force their children inline
                    items = [to_compact_json(v, indent, level + 1, inline=True) for v in obj]
                    return "[\n" + ",\n".join(f"{next_space}{it}" for it in items) + "\n" + space + "]"


            elif isinstance(obj, str):
                return json.dumps(obj)

            else:
                # numbers, booleans, None
                return str(obj).lower()

        # ==================== MAIN - MOLECULES  =======================
        # Setup force field
        self._ff = ForceField(ffname, log=self._logger)

        self.read_charges_types(chargefile)

        molecules = []
        offset = 0
        for ikey, igraph in self._graph_unique_mol.items():

            G = igraph[0]

            ffmol = G.name + '_template.json'
            m = "\t\t Writting {}".format(ffmol)
            print(m) if self._logger is None else self._logger.info(m)

            # Section Critical constants =======================================
            iline_json = empty_molecule_json()
            if self._tc == 0.0 or self._omega == 0.0 or self._pc == 0.0:
                iline_json["CriticalTemperature"] =  round(self._tc, 3)
                iline_json["CriticalPressure"] =  round(self._pc, 3)
                iline_json["AcentricFactor"] =  round(self._omega, 3)
                m = "\n\t\t Critical values are not provided for system {}".format(ffmol)
                print(m) if self._logger is None else self._logger.warn(m)
            else:
                iline_json["CriticalTemperature"] =  round(self._tc, 3)
                iline_json["CriticalPressure"] =  round(self._pc, 3)
                iline_json["AcentricFactor"] =  round(self._omega, 3)

            iline_json["Type"] =  "flexible"

            # Section Number of atoms =======================================
            for iatom_idx in sorted(G):
                iatom = self._mol_mda.atoms[iatom_idx]
                atomtype = self._read_typeatoms[iatom_idx]
                atomname = iatom.name
                # print(iatom_idx, atomname, atomtype)
                iline_json["pseudoAtoms"].append([atomtype, [round(iatom.position[0],3),
                                                             round(iatom.position[1],3),
                                                             round(iatom.position[2],3)]])


            # Topology
            iline_bonds, n_bonds,  list_bonds_json, list_conectivity = self._setup_bonds(G, offset)
            iline_angles, n_angles, list_bends_json = self._setup_angles(G, offset)
            iline_dihedrals, n_dihedrals, list_dihedrals_json = self._setup_dihedrals(G, offset)
            iline_impropers, n_impropers = self._setup_impropers(G, offset, improper_file)

            iline_json["Connectivity"] = list_conectivity
            iline_json["Bonds"] = list_bonds_json
            iline_json["Bends"] = list_bends_json
            iline_json["Torsions"] = list_dihedrals_json
            iline_json["VanDerWaals"] = "auto"

            if n_impropers > 0:
                m = "\n\t\t Improper (out-of-bend) angles are not still implemented in RASPA-3"
                print(m) if self._logger is None else self._logger.warn(m)


            compacted = to_compact_json(iline_json, indent=2, inline=True)
            # Write the file
            with (open(ffmol, 'w')) as ff:
                #json.dump(iline_json, ff, indent=2)
                ff.write(compacted)

        # ==================== MAIN - FORCEFIELD  =======================
        m = "\t\t Writting force_field.json for {}\n".format(ffmol)
        print(m) if self._logger is None else self._logger.info(m)

        iline_ff_json = empty_forcefield_json()

        if ffname == "OPLSAA":
            with open("pseudo_atoms.def", 'w') as ff:

                n_pseudo_atoms = len(self._ff._atomtypes_dict)
                iline_ff_json["MixingRule"] = "Jorgensen"
                iline_ff_json["TruncationMethod"] = "shifted"
                iline_ff_json["TailCorrections"] = False
                iline_ff_json["CutOffVDW"] = 12.0

                for key, item in self._ff._atomtypes_dict.items():
                    if key not in self._read_typeatoms:
                        continue
                    try:
                        charge = self._read_dict_charges[key] # Charge
                        ischargesetup = True
                    except KeyError:
                        charge = item[4]                      # Charge
                        ischargesetup = False
                    iline_ff_json["PseudoAtoms"].append({
                        "name": key,
                        "framework": False,
                        "print_to_output": True,
                        "element": item[0],
                        "print_as": item[0],
                        "mass": item[3],
                        "charge": charge,
                        "source": "OPLS Gromacs Files. Not all charges from FF: {}".format(ischargesetup),
                    })
                    iline_ff_json["SelfInteractions"].append({
                        "name": key,
                        "type": "lennard-jones",
                        "parameters" : [round(item[6]/self._kb_kjmol,7), round(item[5], 6)]
                    })
        else:
            m = "\n\t\t Force field is not implemented in RASPA3.\n"
            m += "\t\t Force field: {}\n".format(ffname)
            m += "\t\t Implement first the force field in the code."
            print(m) if self._logger is None else self._logger.warn(m)

        compacted = to_compact_json(iline_ff_json, indent=2, inline=True)
        # Write the file
        with (open("force_field.json", 'w')) as ff:
            #json.dump(iline_json, ff, indent=2)
            ff.write(compacted)

        iline_sim_json = empty_simulation_json()
        compacted = to_compact_json(iline_sim_json, indent=2, inline=True)
        # Write the file
        with (open("simulation_template.json", 'w')) as ff:
            #json.dump(iline_json, ff, indent=2)
            ff.write(compacted)

    # ====================================================================
    def write_pseudo_atoms(self, ffname):

        """

        Args:
            ffname:

        Returns:

        """

        m = "\t\t Writting pseudo_atoms.def for {}\n".format(ffname)
        print(m) if self._logger is None else self._logger.info(m)

        if ffname == "OPLSAA":
            with open("pseudo_atoms.def", 'w') as ff:

                n_pseudo_atoms = len(self._ff._atomtypes_dict)

                iline = "#number of pseudo_atoms\n"
                iline += "{0:<4d}\n".format(n_pseudo_atoms)
                iline += "#type      print as chem  oxidation   mass         charge "\
                          "polarization B-factor radii  connectivity anisotropic anisotropic-type   tinker-type\n"


                for key, item in self._ff._atomtypes_dict.items():

                    if key not in self._read_typeatoms:
                        continue
                    ii = "{0:<12s} ".format(key)                                    # Type
                    ii += "{0:<3s} ".format("yes")                                  # Print
                    ii += "{0:<3s} ".format(item[0])                                # As
                    ii += "{0:<3s} ".format(item[0])                                # Chem
                    ii += "{0:>2d}     ".format(0)                                  # Oxidation
                    ii += "{0:>12.5f} ".format(item[3])                             # Mass
                    try:
                        ii += "{0:>12.4f}    ".format(self._read_dict_charges[key]) # Charge
                    except KeyError:
                        ii += "{0:>12.4f}    ".format(item[4])                      # Charge
                    ii += "{0:>4.1f}       ".format(0.0)                            # Polarization
                    ii += "{0:>4.1f}     ".format(1.0)                              # B-factor
                    ii += "{0:>4.3f}     ".format(element_cov_radius[item[0]])      # Radii
                    ii += "{0:>2d}          ".format(0)                             # Connectivity
                    ii += "{0:>2d}           ".format(0)                            # Anisotropic
                    ii += "{0:>8s}           ".format("relative")                   # Anisotropic-type
                    ii += "{0:>2d}".format(0)                                       # Tinker-type
                    ii += "\n"
                    iline += ii

                ff.writelines(iline)

        else:
            m = "\n\t\t Force field is not implemented in RASPA.\n"
            m += "\t\t Force field: {}\n".format(ffname)
            m += "\t\t Implement first the force field in the code."
            print(m) if self._logger is None else self._logger.warn(m)

    # ====================================================================
    def write_ff_vdw(self, ffname):

        """

        Args:
            ffname:

        Returns:

        """

        m = "\t\t Writting pseudo_atoms.def for {}\n".format(ffname)
        print(m) if self._logger is None else self._logger.info(m)

        if ffname == "OPLSAA":
            with open("force_field_mixing_rules.def", 'w') as ff:

                iline = "# general rule for shifted vs truncated\n"
                iline += "shifted\n"
                iline += "# general rule tailcorrections\n"
                iline += "no\n"
                iline += "# number of defined interactions\n"
                iline += "{0:<4d}\n".format(len(self._ff._atomtypes_dict))
                iline += "# type interaction, parameters.\n"
                for key, item in self._ff._atomtypes_dict.items():
                    if key not in self._read_typeatoms:
                        continue
                    ii = "{0:<12s}  ".format(key)                        # Type
                    ii += "{0:>12s}  ".format("lennard-jones")           # Interaction
                    ii += "{0:>12.7f}  ".format(item[6]/self._kb_kjmol)  # epsilon (K)
                    ii += "{0:>8.6f}  ".format(item[5])                 # sigma (A)
                    ii += "\n"
                    iline += ii
                iline += "#  general mixing rule for Lennard-Jones\n"
                iline += "Jorgensen\n"
                ff.writelines(iline)

        else:
            m = "\n\t\t Force field is not implemented in RASPA.\n"
            m += "\t\t Force field: {}\n".format(ffname)
            m += "\t\t Implement first the force field in the code."
            print(m) if self._logger is None else self._logger.warn(m)
