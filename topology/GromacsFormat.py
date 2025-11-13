import datetime
import os
import re
import MDAnalysis as Mda
from collections import defaultdict
from MDAnalysis.topology.guessers import guess_improper_dihedrals
from topology.BasicFormat import BasicFormat
from topology.force_fields.ForceField import ForceField


class GromacsFormat(BasicFormat):

    """
    GromacsFormat Class
    """

    # ====================================================================
    def __init__(self, mol_mda, logger=None):

        # Parent class
        super().__init__(mol_mda, logger)

        self._graph_mol, self._graph_mol_components, self._graph_unique_mol, \
            self._graph_component_nodes, self._nmols, \
            self._nmols_unique = self._number_molecules_graph()

        self._molecular_names()


    # ====================================================================
    def write_gro(self):

        """

        Returns:

        """

        ffname_pattern = os.path.splitext(os.path.split(self._mol_mda.trajectory.filename)[-1])[0]
        self._mol_mda.atoms.write(ffname_pattern+'.gro', reindex=False)
        m = "\t\t Writting {0:s}".format(ffname_pattern+'.gro\n')
        print(m) if self._logger is None else self._logger.info(m)

    # ====================================================================
    def write_mdp(self):

        """

        Returns:

        """

        ffname_pattern = os.path.splitext(os.path.split(self._mol_mda.trajectory.filename)[-1])[0]
        ffname = ffname_pattern+'_min.mdp'

        mdp_lines = '; minim.mdp - used as input into grompp to generate em.tpr\n'
        mdp_lines += '; Parameters describing what to do, when to stop and what to save\n'
        mdp_lines += 'integrator  = steep         ; Algorithm (steep = steepest descent minimization)\n'
        mdp_lines += 'emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm\n'
        mdp_lines += 'emstep      = 0.01          ; Minimization step size\n'
        mdp_lines += 'nsteps      = 50000         ; Maximum number of (minimization) steps to perform\n'
        mdp_lines += '\n'
        mdp_lines += '; Parameters describing how to find the neighbors and how to calculate the interactions\n'
        mdp_lines += 'nstlist         = 40         ; Frequency to update the neighbor list and long range forces\n'
        mdp_lines += 'cutoff-scheme   = Verlet    ; Buffered neighbor searching\n'
        mdp_lines += 'ns_type         = grid      ; Method to determine neighbor list (simple, grid)\n'
        mdp_lines += 'coulombtype     = PME       ; Treatment of long range electrostatic interactions\n'
        mdp_lines += 'rcoulomb        = 1.0       ; Short-range electrostatic cut-off\n'
        mdp_lines += 'rvdw            = 1.0       ; Short-range Van der Waals cut-off\n'
        mdp_lines += 'pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions\n'

        with open(ffname, 'w') as ff:
            ff.writelines(mdp_lines)

        m = "\t\t Writting {0:s}".format(ffname)
        print(m) if self._logger is None else self._logger.info(m)

    # ====================================================================
    def write_top(self, ffname, chargefile, improper_file):

        """

        Args:
            ffname:
            chargefile:
            improper_file:

        Returns:

        """

        # Check force field
        self._ff = ForceField(ffname, log=self._logger)
        pattern_top = os.path.splitext(os.path.split(self._mol_mda.trajectory.filename)[-1])[0]

        self.read_charges_types(chargefile)

        # Default section
        iline_top = '; Force field from typing_gases_to_VLE_cmd.py\n'
        iline_top += '; # 1 (Lennard-Jones) or 2 (Buckingham)\n'
        iline_top += '; # 1 and 3 (Geometric rule) and 2 (Arithmetic rule)\n'
        iline_top += '; # Generate pairs (yes/no). If yes generates pairs not present in [pairtypes]\n'
        iline_top += '; # Factor to multiply LJ 1-4 interctions, default 1. Not used if genpairs=no\n'
        iline_top += '; # Factor to multiply electrostatic 1-4 interctions, default 1. Always used.\n'
        iline_top += '[defaults]\n'
        iline_top += '; nbfunc    comb-rule   gen-pairs   fudgeLJ     fudgeQQ\n'
        iline_top += '{0:2d}     {1:2d}     {2:3s}     {3:.1f}   {4:.1f}\n'.\
            format(self._ff._nbfunctiontype, self._ff._combrule, self._ff._genpairs,
                   self._ff._fudge_lj, self._ff._fudge_coulomb)

        # Include section
        iline_top += "\n"
        iline_top += '#include "{0:s}"\n'.format(self._ff._general_itp_name)

        molecules = []
        offset = 0
        for ikey, igraph in self._graph_unique_mol.items():
            itpname = self._write_unique_itp(igraph[0], offset, improper_file)
            iline_top += '#include "{0:s}"\n'.format(itpname+".itp")
            molecules.append([itpname, len(igraph)])
            offset += len(igraph[0].nodes)

        iline_top += "\n"
        iline_top += "[ system ]\n"
        iline_top += "{0:s}\n".format(pattern_top)

        iline_top += "\n"
        iline_top += "[ molecules ]\n"
        iline_top += ";system molecules\n"
        for igraph in self._graph_mol_components:
            iline_top +=  igraph.name + "        {0:s} \n".format("#NMOLS#")

        with open(pattern_top+".top", "w") as fout:
            fout.writelines(iline_top)

    # ====================================================================
    def _write_unique_itp(self, igraph, offset=0, improper_file=None):

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # Neccesary to generate the pairlist
        # components = list(nx.connected_components(self._graph_mol))
        # subgraphs = [self._graph_mol.subgraph(component).copy() for component in components]
        lines_itp = ""
        if igraph.name == "":
            pattern_itp = "MOL_{0:03d}".format(self._imol)
        else:
            pattern_itp = igraph.name
        self._imol += 1

        m = "\t\t Name of the molecule {0:3d}: {1:s}".format(self._imol, pattern_itp)
        print(m) if self._logger is None else self._logger.info(m)

        # header
        lines_itp += "; itp created with typing_gases_to_VLE_cmd.py program from Topology library\n"
        lines_itp += ";                      Javier Ramos (IEM-CSIC)                        \n"
        lines_itp += ";             Created:    {}                                   \n".format(now)

        # [ moleculetype ]
        lines_itp += "\n"
        lines_itp += "[ moleculetype ]\n"
        lines_itp += "; Name            nrexcl\n"
        lines_itp += "{0:s}             {1:d}\n".format(pattern_itp, self._ff._nrexcl)

        # [ atoms ]
        lines_itp += "\n"
        lines_itp += "[ atoms ]\n"
        lines_itp += ";   nr       type  resnr residue  atom   cgnr     charge        mass\n"
        qtot = 0.0
        resid_current = 0
        first_res = True
        qtot_res = defaultdict()
        improper_list = []
        iatom_idx_local = 0
        resid_local = 0

        for iatom_idx in sorted(igraph):

            iatom = self._mol_mda.atoms[iatom_idx]
            resname = iatom.resname
            # atomtype = iatom.name
            atomtype = self._read_typeatoms[iatom_idx]
            atomname = iatom.name
            resid = iatom.resid

            # Check if atom type is defined in the force field
            try:
                self._ff._atomtypes_dict[atomtype]
            except KeyError:
                m = "\t\t Atom type {0:s} (index {1:d}) is not defined in the ForceField class.\n".\
                    format(atomtype, iatom_idx_local)
                m1 = "\n\t\t" + "*" * len(m)+"\n"
                m2 = "\t\t" + "*" * len(m)+"\n"
                m += "\t\t Check the atom type in the PDB file ({})\n".format(self._mol_mda.trajectory.filename)
                print(m1+m+m2) if self._logger is None else self._logger.error(m1+m+m2)
                exit()

            # Residues
            if first_res:
                resid_current = resid
                first_res = False
                qtot_res[resid_local] = 0.0
                lines_itp += "; residue {0:>8d}\n".format(resid_local+1)
            if resid != resid_current:
                resid_local += 1
                resid_current = resid
                qtot_res[resid_local] = 0.0
                lines_itp += "; residue {0:>8d}\n".format(resid_local+1)

            charge = self._read_charges[iatom.index]
            mass = self._ff._atomtypes_dict[atomtype][3]
            qtot += charge
            qtot_res[resid_local] += charge
            lines_itp += "{0:^12d}{1:^7s}{2:>8d}{3:>8s}  {4:>8s}{5:>12d}{6:>13.7f}{7:>11.4f} ; qtot {8:>13.3f}\n" \
                .format(iatom_idx_local + 1, atomtype, resid_local+1, resname,
                        atomname, resid_local+1, charge, mass, qtot)
            iatom_idx_local += 1

        # [ bonds ]
        lines_itp += "\n"
        lines_itp += "[ bonds ]\n"
        lines_itp += "; ai    aj funct bo(nm) kb(kJ/molnm2)\n"
        nbonds_corrected = 0

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

            nbonds_corrected += 1
            lines_itp += "{0:>10d} {1:>10d} {2:>10d} {3:.4f} {4:14.3f}\n".\
                          format(iat + 1 - offset, jat + 1 - offset,
                       self._ff._bondtypes_dict[blabel_d][0],
                       self._ff._bondtypes_dict[blabel_d][1],
                       self._ff._bondtypes_dict[blabel_d][2])

        # [ angles ]
        lines_itp += "\n"
        lines_itp += "[ angles ]\n"
        lines_itp += ";  ai    aj    ak funct\n"
        nangles_corrected = 0
        angles = sorted(list(Mda.topology.guessers.guess_angles(self._mol_mda.bonds)))
        self._mol_mda.add_TopologyAttr('angles', angles)
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

            nangles_corrected += 1
            lines_itp += "{0:>10d} {1:>10d} {2:>10d} {3:>10d} {4:12.3f} {5:12.3f}\n".\
                format(iat+1-offset, jat+1-offset, kat+1-offset,
                        self._ff._angletypes_dict[alabel_d][0],
                        self._ff._angletypes_dict[alabel_d][1],
                        self._ff._angletypes_dict[alabel_d][2])


        # [ dihedrals ]
        lines_itp += "\n"
        lines_itp += "[ dihedrals ]\n"
        lines_itp += "; ai    aj  ak  al  funct\n"
        ndihedrals_corrected = 0
        dihedrals = sorted(list(Mda.topology.guessers.guess_dihedrals(self._mol_mda.angles)))
        self._mol_mda.add_TopologyAttr('dihedrals', dihedrals)

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

            ndihedrals_corrected += 1
            lines_itp += ("{0:>10d} {1:>10d} {2:>10d} {3:>10d} {4:>10d} ".
                          format(iat+1-offset, jat+1-offset, kat+1-offset, lat+1-offset,
                                 self._ff._dihedraltypes_dict[dlabel_d][0]))
            for item in self._ff._dihedraltypes_dict[dlabel_d][1]:
                lines_itp += "{0:10.5f}  ".format(item)
            lines_itp += "\n"

        # [ impropers ]
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
                    lines_itp += "; Impropers\n"
                    for j in range(idx, idx+b):
                        iline = lines[j].split()
                        if int(iline[4]) == 4:
                            lines_itp += "{0:8d} {1:8d} {2:8d} {3:8d}  {4:1d}   {5:7.3f} {6:s} {7:d}\n".\
                                              format(int(iline[0]), int(iline[1]), int(iline[2]), int(iline[3]),
                                                     int(iline[4]), float(iline[5]), iline[6], int(iline[7]))
                        elif int(iline[4]) == 2:
                            lines_itp += "{0:8d} {1:8d} {2:8d} {3:8d}  {4:1d}   {5:7.3f} {6:s}\n".\
                                              format(int(iline[0]), int(iline[1]), int(iline[2]), int(iline[3]),
                                                     int(iline[4]), float(iline[5]), iline[6])
                        else:
                            m = "\n\t\tERROR: Improper type {} is not implemented".format(int(iline[4]))
                            print(m) if self._logger is None else self._logger.error(m)
                    improper_list.append(improper_lines)
                    idx += b

        # [ pairs ]
        lines_itp += "\n"
        lines_itp += "[ pairs ]\n"
        lines_itp += "; ai    aj funct\n"
        pair_function = 1
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
            lines_itp += "{0:>10d} {1:>10d} {2:>10d}\n".format(ipair[0]+1-offset, ipair[1]+1-offset, pair_function)

        # Write file
        fname = "{0:4s}.itp".format(pattern_itp)
        m = "\t\t Writting {0:s}".format(fname)
        print(m) if self._logger is None else self._logger.info(m)
        lines_itp += "\n"
        with open(fname, "w") as fout:
            fout.writelines(lines_itp)
            # fout.writelines(lines_rest_labels)
            # fout.writelines(lines_rest_atoms)
            # fout.writelines("#endif\n")

        return pattern_itp


