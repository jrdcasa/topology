import datetime
import os.path
import re
import sys
import warnings
from collections import defaultdict

import MDAnalysis as Mda
import networkx as nx

from utils.Logger import init_logger
from utils.parse_arguments import parse_arguments_typing_carb, print_header_type_cb, command_info

warnings.filterwarnings('ignore')

ions_list = ["NA", "LI"]
set_atoms_CMN_residue = ["Hm1", "Hm2", "Cm1", "Cm2", "Om1", "Om2", "Na"]
improper_CMN_residue = [["Cm2", "Cm1", "Om1", "Om2", 2, 0.0, 401.664]]
# charmm/toppar/top_all36_carb.rtf --> BGLC and PATCHES for the disaccharides
# Format: atomtype, [charge if is not link atom, charge if is link atom, mass]
name_to_type_bglc_charmm_nopatch = {"C1": ["CC3162", 0.3400, 12.0110],
                                    "C2": ["CC3161", 0.1400, 12.0110],
                                    "C3": ["CC3161", 0.1400, 12.0110],
                                    "C4": ["CC3161", 0.1400, 12.0110],
                                    "C5": ["CC3163", 0.1100, 12.0110],
                                    "C6": ["CC321", 0.0500, 12.0110],
                                    "O1": ["OC311", -0.650, 15.9994],
                                    "O2": ["OC311", -0.650, 15.9994],
                                    "O3": ["OC311", -0.650, 15.9994],
                                    "O4": ["OC311", -0.650, 15.9994],
                                    "O5": ["OC3C61", -0.400, 15.9994],
                                    "O6": ["OC311", -0.650, 15.9994],
                                    "O6i": ["OC311", -0.980, 15.9994],  # Inventado
                                    "H1": ["HCA1", 0.090, 1.0080],
                                    "H2": ["HCA1", 0.090, 1.0080],
                                    "H3": ["HCA1", 0.090, 1.0080],
                                    "H4": ["HCA1", 0.090, 1.0080],
                                    "H5": ["HCA1", 0.090, 1.0080],
                                    "H61": ["HCA2", 0.090, 1.0080],
                                    "H62": ["HCA2", 0.090, 1.0080],
                                    "H6A": ["HCA2", 0.090, 1.0080],
                                    "H6B": ["HCA2", 0.090, 1.0080],
                                    "HO1": ["HCP1", 0.420, 1.0080],
                                    "HO2": ["HCP1", 0.420, 1.0080],
                                    "HO3": ["HCP1", 0.420, 1.0080],
                                    "HO4": ["HCP1", 0.420, 1.0080],
                                    "HO6": ["HCP1", 0.420, 1.0080],
                                    "Cm1": ["CC321", 0.000, 12.0110],
                                    "Hm1": ["HCA2", 0.09, 1.0080],
                                    "Hm2": ["HCA2", 0.09, 1.0080],
                                    "Cm2": ["CC2O2", 0.520, 12.0110],
                                    "Om1": ["OC2D2", -0.760, 15.9994],
                                    "Om2": ["OC2D2", -0.760, 15.9994],
                                    "Na": ["SOD", 1.000, 22.98977],
                                    }

name_to_type_bglc_charmm_patch = {"C1": ["CC3162", 0.2900, 12.0110],
                                  "C2": ["CC3161", 0.0900, 12.0110],
                                  "C3": ["CC3161", 0.0900, 12.0110],
                                  "C4": ["CC3161", 0.0900, 12.0110],
                                  "C6": ["CC321", 0.0000, 12.0110],
                                  "O1": ["OC301", -0.360, 15.9994],
                                  "O2": ["OC301", -0.360, 15.9994],
                                  "O3": ["OC301", -0.360, 15.9994],
                                  "O4": ["OC301", -0.360, 15.9994],
                                  "O6": ["OC301", -0.360, 15.9994],
                                  "Cm1": ["CC321", 0.000, 1.0080],
                                  }

atoms_backbone = ["C1", "C2", "C3", "C4", "C5", "O5"]
atoms_sidechain = ["O1", "O2", "O3", "O4", "C6", "O6"]

# =============================================================================
def find_paths_of_length_n(graph, length):


    def find_paths(graphf, start, target_length):

        stack = [(start, [start])]
        pathsf = []

        while stack:
            current_node, current_path = stack.pop()

            if len(current_path) == target_length:
                pathsf.append(current_path)
                continue

            for neighbor in graphf.neighbors(current_node):
                if neighbor not in current_path:
                    stack.append((neighbor, current_path + [neighbor]))

        return pathsf

    # Convert the directed graph to an undirected graph
    undirected_graph = nx.Graph(graph)
    all_paths_of_length = []
    for node in graph.nodes:
        paths = find_paths(undirected_graph, node, length)
        all_paths_of_length.extend(paths)

    return all_paths_of_length


# =============================================================================
def number_polymer_chains(u):
    graph_mol = nx.Graph()
    graph_component_nodes = defaultdict()

    # Add nodes and edges
    for idx_bond, ibond in enumerate(u.bonds):
        idx1 = ibond.indices[0]
        idx2 = ibond.indices[1]
        graph_mol.add_node(idx1)
        graph_mol.add_node(idx2)
        graph_mol.add_edge(idx1, idx2)

    # Set attributtes for each node
    attrs = defaultdict(dict)
    oxygen_links = ["O1", "O2", "O3", "O4", "O6"]

    for inode in graph_mol.nodes:
        attrs[inode] = {"resid": None, "name": None, "link": [], "patch_atom": False}

    for inode in graph_mol.nodes:
        resid = u.atoms[inode].resid
        name = u.atoms[inode].name
        attrs[inode]["resid"] = resid
        attrs[inode]["name"] = name

        if name in oxygen_links:
            neighs = [n for n in graph_mol.neighbors(inode)]
            neighs_element = [u.atoms[i].element for i in graph_mol.neighbors(inode)]
            if 'H' not in neighs_element:
                attrs[inode]["patch_atom"] = True
                for jnode in neighs:
                    attrs[jnode]["patch_atom"] = True

    nx.set_node_attributes(graph_mol, attrs)

    # Number of disconnected graphs
    nmols = nx.number_connected_components(graph_mol)
    for idx, imol in enumerate(nx.connected_components(graph_mol)):
        graph_component_nodes[idx] = sorted(imol)


    return graph_mol, graph_component_nodes, nmols


# =============================================================================
def write_itp(imol, graph_mol, graph_component_nodes, u, logger=None):



    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    lines_itp = ""

    # Neccesary to generate the pairlist
    components = list(nx.connected_components(graph_mol))
    subgraphs = [graph_mol.subgraph(component).copy() for component in components]

    # header
    lines_itp += "; itp created with typing_carb_cmd.py program from Topology library\n"
    lines_itp += ";                   Javier Ramos (IEM-CSIC)                        \n"
    lines_itp += ";           Created:    {}                                   \n".format(now)

    # [ moleculetype ]
    lines_itp += "\n"
    lines_itp += "[ moleculetype ]\n"
    lines_itp += "; Name            nrexcl\n"
    lines_itp += "CARB{0:02d}             3\n".format(imol)

    # [ atoms ]
    lines_itp += "\n"
    lines_itp += "[ atoms ]\n"
    lines_itp += ";   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB\n"
    lines_rest_atoms = ""
    qtot = 0.0
    resid_current = 0
    first_res = True
    qtot_res = defaultdict()
    improper_list = []
    iatom_idx_local = 0
    resid_local = 0
    for iatom_idx in graph_component_nodes[imol]:
        iatom = u.atoms[iatom_idx]

        atomtype = iatom.name
        resid = iatom.resid
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

        atomname = iatom.name
        resname = iatom.resname
        ispatch = graph_mol.nodes[iatom_idx]["patch_atom"]

        # Carboxymethyl group impropers
        improper_sublist = [-1, -1, -1, -1]
        #improper_list.append()
        if atomname == "Cm2":
            neighs = [n for n in graph_mol.neighbors(iatom_idx)]
            improper_sublist[0] = iatom_idx+1
            for ineigh in neighs:
                tt = u.atoms[ineigh].name
                if tt == "Cm1":
                    improper_sublist[1] = ineigh+1
                elif tt == "Om1":
                    improper_sublist[2] = ineigh+1
                elif tt == "Om2":
                    improper_sublist[3] = ineigh+1
                else:
                    continue
        if not any(elem == -1 for elem in improper_sublist):
            improper_sublist.append(2)
            improper_sublist.append(0.0)
            improper_sublist.append(401.664)
            improper_list.append(improper_sublist)

        if resname == "BGL" or resname == "BGLC":
            try:
                if ispatch:
                    atomtype = name_to_type_bglc_charmm_patch[atomname][0]
                    charge = float(name_to_type_bglc_charmm_patch[atomname][1])
                    mass = float(name_to_type_bglc_charmm_patch[atomname][2])
                else:
                    atomtype = name_to_type_bglc_charmm_nopatch[atomname][0]
                    charge = float(name_to_type_bglc_charmm_nopatch[atomname][1])
                    mass = float(name_to_type_bglc_charmm_nopatch[atomname][2])
            except (KeyError, IndexError):
                m = "\t\t\t ERROR!!!! atomtype {} is not found.".format(atomtype)
                charge = -99.9
                mass = 0.0
                print(m) if logger is None else logger.error(m)
                exit()
        else:
            charge = -99.9
            mass = 0.0
            m = "\t\t\t ERROR!!!! resname {} in PDB is unknown.".format(resname)
            print(m) if logger is None else logger.error(m)
            exit()
        qtot += charge
        qtot_res[resid_local] += charge
        lines_itp += "{0:^12d}{1:^7s}{2:>8d}{3:>8s}{4:>8s}{5:>12d}{6:>13.7f}{7:>11.4f} ; qtot {8:>13.3f}\n" \
            .format(iatom_idx_local + 1, atomtype, resid_local+1, resname, atomname, resid_local+1, charge, mass, qtot)
        iatom_idx_local += 1

        if atomname in atoms_backbone:
            lines_rest_atoms += "{0:^12d}  1  POSRES_FC_BB  POSRES_FC_BB   POSRES_FC_BB\n".format(iatom_idx_local + 1)
        elif atomname in atoms_sidechain:
            lines_rest_atoms += "{0:^12d}  1  POSRES_FC_SC  POSRES_FC_SC   POSRES_FC_SC\n".format(iatom_idx_local + 1)

    for ires in range(0, resid_local):
        substring_to_replace = "; residue {0:>8d}\n".format(ires+1)
        new_substring = "; residue {0:>8d} q {1:>6.3f} \n".format(ires, qtot_res[ires])
        lines_itp = re.sub(substring_to_replace, new_substring, lines_itp)

    # Offset for different molecules
    offset = 0
    for i in range(0, imol):
        offset += len(graph_component_nodes[i])

    # [ bonds ]
    lines_itp += "\n"
    lines_itp += "[ bonds ]\n"
    lines_itp += "; ai    aj funct\n"
    bond_function = 1
    nbonds_corrected = 0

    for ibond in u.bonds:
        iat = ibond.indices[1]
        jat = ibond.indices[0]
        if iat not in graph_component_nodes[imol]:
            continue
        itype = ibond.type[0]
        jtype = ibond.type[1]
        # Do not write bond for Na, Li, ...
        if itype.upper() in ions_list or jtype.upper() in ions_list:
            continue
        nbonds_corrected += 1
        lines_itp += "{0:>10d} {1:>10d} {2:>10d}\n".format(iat+1-offset, jat+1-offset, bond_function)

    # [ pairs ]
    lines_itp += "\n"
    lines_itp += "[ pairs ]\n"
    lines_itp += "; ai    aj funct\n"
    pair_function = 1
    paths_of_length_3 = find_paths_of_length_n(subgraphs[imol], 4)
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

    # [ angles ]
    lines_itp += "\n"
    lines_itp += "[ angles ]\n"
    lines_itp += ";  ai    aj    ak funct\n"
    angle_function = 5
    nangles_corrected = 0
    angles = sorted(list(Mda.topology.guessers.guess_angles(u.bonds)))
    nangles_guess = len(angles)
    u.add_TopologyAttr('angles', angles)
    for iangle in u.angles:
        iat = iangle.indices[0]
        jat = iangle.indices[1]
        kat = iangle.indices[2]
        itype = iangle.type[0]
        jtype = iangle.type[1]
        ktype = iangle.type[2]
        if iat not in graph_component_nodes[imol]:
            continue
        # Do not write angles for Na, Li, ...
        if itype.upper() in ions_list or \
            jtype.upper() in ions_list or \
            ktype.upper() in ions_list:
            continue
        nangles_corrected += 1
        lines_itp += "{0:>10d} {1:>10d} {2:>10d} {3:>10d}\n".format(iat+1-offset, jat+1-offset,
                                                                    kat+1-offset, angle_function)

    # [ dihedrals ]
    lines_itp += "\n"
    lines_itp += "[ dihedrals ]\n"
    lines_itp += "; ai    aj  ak  al  funct\n"
    dihedral_function = 9
    ndihedrals_corrected = 0
    dihedrals = sorted(list(Mda.topology.guessers.guess_dihedrals(u.angles)))
    ndihedrals_guess = len(dihedrals)
    u.add_TopologyAttr('dihedrals', dihedrals)

    for idih in u.dihedrals:
        iat = idih.indices[0]
        jat = idih.indices[1]
        kat = idih.indices[2]
        lat = idih.indices[3]
        itype = idih.type[0]
        jtype = idih.type[1]
        ktype = idih.type[2]
        ltype = idih.type[3]
        if iat not in graph_component_nodes[imol]:
            continue
        # Do not write angles for Na, Li, ...
        if itype.upper() in ions_list or \
            jtype.upper() in ions_list or \
            ktype.upper() in ions_list or \
            ltype.upper() in ions_list:
            continue
        ndihedrals_corrected += 1
        lines_itp += "{0:>10d} {1:>10d} {2:>10d} {3:>10d} {4:>10d}\n".format(iat+1-offset, jat+1-offset,
                                                                             kat+1-offset, lat+1-offset,
                                                                             dihedral_function)

    nimpropers_corrected = 0
    if len(improper_list) != 0:
        lines_itp += "; impropers\n"
        for imp in improper_list:
            iat = imp[0]
            jat = imp[1]
            kat = imp[2]
            lat = imp[3]
            ff = imp[4]
            angle = imp[5]
            kk = imp[6]
            lines_itp += "{0:>10d} {1:>10d} {2:>10d} {3:>10d} {4:>10d} {5:>4.1f} {6:10.3f}\n".\
                format(iat-offset, jat-offset, kat-offset, lat-offset, ff, angle, kk)
            nimpropers_corrected += 1

    lines_rest_labels = "#ifdef POSRES\n"
    lines_rest_labels +="[ position_restraints ]\n"

    m = "\t\t ****** PDB INFO (corrected) ******\n"
    m = "\t\t ******    imol = {}         ******\n".format(imol+1)
    m += "\t\t; Bonded parameters with ions are removed\n"
    m += "\t\tNumber of bonds (corrected)    : {}\n".format(nbonds_corrected)
    m += "\t\tNumber of angles (guess)       : {} by MDA using all atoms not components\n".format(nangles_guess)
    m += "\t\tNumber of angles (corrected)   : {} angles in each component or molecule\n".format(nangles_corrected)
    m += "\t\tNumber of dihedrals (guess)    : {}\n".format(ndihedrals_guess)
    m += "\t\tNumber of dihedrals (corrected): {}\n".format(ndihedrals_corrected)
    m += "\t\tNumber of impropers (corrected): {}\n".format(nimpropers_corrected)
    print(m) if logger is None else logger.error(m)

    # Write file
    fname = "CARB{0:02d}.itp".format(imol)
    lines_itp += "\n"
    with open(fname, "w") as fout:
        fout.writelines(lines_itp)
        fout.writelines(lines_rest_labels)
        fout.writelines(lines_rest_atoms)
        fout.writelines("#endif\n")

    molname_itp = "CARB{0:02d}".format(imol)
    return fname,  molname_itp


# =============================================================================
def gromacs_mdp_templates():

    # Minimization
    lines_mdp ='define                  = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0 -DDIHRES -DDIHRES_FC=4.0\n'
    lines_mdp +='integrator              = steep\n'
    lines_mdp +='emtol                   = 1000.0\n'
    lines_mdp +='nsteps                  = 5000\n'
    lines_mdp +='nstlist                 = 10\n'
    lines_mdp +='cutoff-scheme           = Verlet\n'
    lines_mdp +='rlist                   = 1.2\n'
    lines_mdp +='vdwtype                 = Cut-off\n'
    lines_mdp +='vdw-modifier            = Force-switch\n'
    lines_mdp +='rvdw_switch             = 1.0\n'
    lines_mdp +='rvdw                    = 1.2\n'
    lines_mdp +='coulombtype             = PME\n'
    lines_mdp +='rcoulomb                = 1.2\n'
    lines_mdp +=';\n'
    lines_mdp +='constraints             = h-bonds\n'
    lines_mdp +='constraint_algorithm    = LINCS\n'

    fname = "minimization_template.mdp"
    with open(fname, "w") as fout:
        fout.writelines(lines_mdp)

    # Equilibration
    lines_mdp ='define                  = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0 -DDIHRES -DDIHRES_FC=4.0\n'
    lines_mdp +='integrator              = md\n'
    lines_mdp +='dt                      = 0.001\n'
    lines_mdp +='nsteps                  = 125000\n'
    lines_mdp +='nstxout-compressed      = 5000\n'
    lines_mdp +='nstxout                 = 0\n'
    lines_mdp +='nstvout                 = 0\n'
    lines_mdp +='nstfout                 = 0\n'
    lines_mdp +='nstcalcenergy           = 100\n'
    lines_mdp +='nstenergy               = 1000\n'
    lines_mdp +='nstlog                  = 1000\n'
    lines_mdp +=';\n'
    lines_mdp +='cutoff-scheme           = Verlet\n'
    lines_mdp +='nstlist                 = 20\n'
    lines_mdp +='rlist                   = 1.2\n'
    lines_mdp +='vdwtype                 = Cut-off\n'
    lines_mdp +='vdw-modifier            = Force-switch\n'
    lines_mdp +='rvdw_switch             = 1.0\n'
    lines_mdp +='rvdw                    = 1.2\n'
    lines_mdp +='coulombtype             = PME\n'
    lines_mdp +='rcoulomb                = 1.2\n'
    lines_mdp +=';\n'
    lines_mdp +='tcoupl                  = v-rescale\n'
    lines_mdp +='tc_grps                 = SOLU SOLV\n'
    lines_mdp +='tau_t                   = 1.0 1.0\n'
    lines_mdp +='ref_t                   = 303.15 303.15\n'
    lines_mdp +=';\n'
    lines_mdp +='constraints             = h-bonds\n'
    lines_mdp +='constraint_algorithm    = LINCS\n'
    lines_mdp +=';\n'
    lines_mdp +='nstcomm                 = 100\n'
    lines_mdp +='comm_mode               = linear\n'
    lines_mdp +='comm_grps               = SOLU SOLV\n'
    lines_mdp +=';\n'
    lines_mdp +='gen-vel                 = yes\n'
    lines_mdp +='gen-temp                = 303.15\n'
    lines_mdp +='gen-seed                = -1\n'

    fname = "equilibration_template.mdp"
    with open(fname, "w") as fout:
        fout.writelines(lines_mdp)

    lines_mdp ='integrator              = md\n'
    lines_mdp +='dt                      = 0.002\n'
    lines_mdp +='nsteps                  = 500000\n'
    lines_mdp +='nstxout-compressed      = 50000\n'
    lines_mdp +='nstxout                 = 0\n'
    lines_mdp +='nstvout                 = 0\n'
    lines_mdp +='nstfout                 = 0\n'
    lines_mdp +='nstcalcenergy           = 100\n'
    lines_mdp +='nstenergy               = 1000\n'
    lines_mdp +='nstlog                  = 1000\n'
    lines_mdp +=';\n'
    lines_mdp +='cutoff-scheme           = Verlet\n'
    lines_mdp +='nstlist                 = 20\n'
    lines_mdp +='vdwtype                 = Cut-off\n'
    lines_mdp +='vdw-modifier            = Force-switch\n'
    lines_mdp +='rvdw_switch             = 1.0\n'
    lines_mdp +='rvdw                    = 1.2\n'
    lines_mdp +='rlist                   = 1.2\n'
    lines_mdp +='rcoulomb                = 1.2\n'
    lines_mdp +='coulombtype             = PME\n'
    lines_mdp +=';\n'
    lines_mdp +='tcoupl                  = v-rescale\n'
    lines_mdp +='tc_grps                 = SOLU SOLV\n'
    lines_mdp +='tau_t                   = 1.0 1.0\n'
    lines_mdp +='ref_t                   = 303.15 303.15\n'
    lines_mdp +=';\n'
    lines_mdp +='pcoupl                  = C-rescale\n'
    lines_mdp +='pcoupltype              = isotropic\n'
    lines_mdp +='tau_p                   = 5.0\n'
    lines_mdp +='compressibility         = 4.5e-5\n'
    lines_mdp +='ref_p                   = 1.0\n'
    lines_mdp +=';\n'
    lines_mdp +='constraints             = h-bonds\n'
    lines_mdp +='constraint_algorithm    = LINCS\n'
    lines_mdp +='continuation            = yes\n'
    lines_mdp +=';\n'
    lines_mdp +='nstcomm                 = 100\n'
    lines_mdp +='comm_mode               = linear\n'
    lines_mdp +='comm_grps               = SOLU SOLV\n'
    lines_mdp +=';\n'

    fname = "production_template.mdp"
    with open(fname, "w") as fout:
        fout.writelines(lines_mdp)


# =============================================================================
def create_gro_top_charmm36(u, boxdim, logger=None):

    resname_dict = defaultdict(list)
    idx_atom_remove_bonds = []

    m = "\t\t ****** PDB INFO ******\n"
    m += "\t\tFilename                    : {}\n".format(u.filename)
    m += "\t\tNumber of atoms             : {}\n".format(len(u.atoms))
    m += "\t\tNumber of residues          : {}\n".format(len(u.residues))
    m += "\t\tNumber of bonds             : {}\n".format(len(u.bonds))
    for ires in u.residues:
        resname_dict[ires.resname].append(ires.resid)
    if boxdim is not None:
        if len(boxdim) == 1:
            u.dimensions = [boxdim[0], boxdim[0], boxdim[0], 90.0, 90.0, 90.0]
        else:
            u.dimensions = [boxdim[0], boxdim[1], boxdim[2], 90.0, 90.0, 90.0]
    m += "\t\tLength box dimensions (Ang) : {} {} {}\n".format(u.dimensions[0], u.dimensions[1], u.dimensions[2])
    m += "\t\tAngle  box dimensions (º)   : {} {} {}".format(u.dimensions[3], u.dimensions[4], u.dimensions[5])

    pattern = os.path.splitext(u.filename)[0]
    u.atoms.write(pattern + '.gro', reindex=False)
    print(m) if logger is None else logger.info(m)

    # Number of different polymer chains
    mol_graph, graph_component_nodes, nmols = number_polymer_chains(u)

    # Setup dictionaries to build itps
    filename_itp = defaultdict()
    molname_itp = defaultdict()
    for imol in range(0, nmols):
        filename_itp[imol], molname_itp[imol] = write_itp(imol, mol_graph, graph_component_nodes, u)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    lines_top = ""
    # header
    lines_top += "; top created with typing_carb_cmd.py program from Topology library\n"
    lines_top += ";                   Javier Ramos (IEM-CSIC)                        \n"
    lines_top += ";           Created:    {}                                   \n".format(now)

    lines_top += '#include "./charmm36_ljpme-jul2022.ff/forcefield.itp"\n'
    lines_top += '#include "./charmm36_ljpme-jul2022.ff/tip3p.itp"\n'
    lines_top += '#include "./charmm36_ljpme-jul2022.ff/ions.itp"\n'
    for i in filename_itp:
        lines_top += '#include "{}"\n'.format(filename_itp[i])

    lines_top += "\n"
    lines_top += "[ system ]\n"
    lines_top += "; Name\n"
    lines_top += "CARB\n"

    lines_top += "\n"
    lines_top += "[ molecules ]\n"
    lines_top += "; Compound        mols\n"
    for i in filename_itp:
        lines_top += '{}             1\n'.format(molname_itp[i])

    with open("topology.top", "w") as ftop:
        ftop.writelines(lines_top)

    gromacs_mdp_templates()


# =============================================================================
def edit_charm36_ff(logger=None):

    """
    Add some terms to the force field charmm36_ljpme-jul2022 by comparison with other terms

    #   Angle UB: CC202 - CC321 - OC301
    #    CC2O2    CC311    OC301     5   103.000000   376.560000   0.00000000         0.00 ; og amop mp2/631gd
         CC2O2    CC301    OC301     5   103.000000   376.560000   0.00000000         0.00 ; from CC2O2 CC311 OC301; og amol ok
         CC2O2   CC3062    OC301     5   103.000000   376.560000   0.00000000         0.00 ; from CC2O2   CC301   OC301


    """
    dir_charmm = "charmm36_ljpme-jul2022.ff"

    # Copy charmm36 directory as parameters
    charmm_path = os.path.split(sys.argv[0])[0]
    charmm_path += "/../data/{}".format(dir_charmm)
    os.system("cp -r {} .".format(charmm_path))

    fname = os.path.join(dir_charmm, "ffbonded.itp")
    with open(fname, 'r') as fileff:
        lines_ffbonded = fileff.readlines()
        # Insert UB Angle
        string_to_replace = "   CC2O2    CC311    OC301     5   " \
                            "103.000000   376.560000   0.00000000         0.00 ; og amop mp2/631gd\n"
        new_substring = "   CC2O2    CC321    OC301     5   " \
                      "103.000000   376.560000   0.00000000         0.00 ; Add cJ\n"
        try:
            index = lines_ffbonded.index(string_to_replace)
            lines_ffbonded.insert(index, new_substring)
        except ValueError:
            m = f"\n\t\t'{string_to_replace}' not found in the list."
            print(m) if logger is None else logger.error(m)
            exit()
        # Insert dihedral CC3163-CC321-OC301-CC321
        string_to_replace3 = '  CC3163    CC321    OC301     CTL2     9     0.000000     2.552240     3 ; \n'
        new_substring1 = '  CC3163    CC321    OC301    CC321     9   180.000000     2.677760     1 ; Add cJ\n'
        new_substring2 = '  CC3163    CC321    OC301    CC321     9   180.000000     0.125520     2 ; Add cJ\n'
        new_substring3 = '  CC3163    CC321    OC301    CC321     9     0.000000     2.552240     3 ; Add cJ\n'
        try:
            index3 = lines_ffbonded.index(string_to_replace3)
            lines_ffbonded.insert(index3+1, new_substring1)
            lines_ffbonded.insert(index3+2, new_substring2)
            lines_ffbonded.insert(index3+3, new_substring3)
        except ValueError:
            m = f"\n\t\t'{string_to_replace3}' not found in the list."
            print(m) if logger is None else logger.error(m)
            exit()

        # Insert dihedral -CC2O2-CC321-OC301-CC321
        string_to_replace = '   CC2O2    CC311    OC301    CC331     9     0.000000     0.836800     3 ' \
                            '; og/sng thp CC321C CC321C OC3C61 CC321C\n'
        new_substring = '   CC2O2    CC321    OC301    CC321     9     0.000000     0.836800     3 ' \
                            '; Add cJ\n'
        try:
            index = lines_ffbonded.index(string_to_replace)
            lines_ffbonded.insert(index+1, new_substring)
        except ValueError:
            m = f"\n\t\t'{string_to_replace}' not found in the list."
            print(m) if logger is None else logger.error(m)
            exit()

        # Insert dihedral -CC2O2-CC321-OC301-CC321
        string_to_replace = '   OC301    CC301    CC2O2    OC2D2     9   180.000000     2.677760     2 ' \
                            '; og amop mp2/ccpvtz\n'
        new_substring1 = '   OC301    CC321    CC2O2    OC2D2     9   180.000000     2.677760     2 ' \
                            '; og amop mp2/ccpvtz\n'
        try:
            index = lines_ffbonded.index(string_to_replace)
            lines_ffbonded.insert(index+1, new_substring1)
        except ValueError:
            m = f"\n\t\t'{string_to_replace}' not found in the list."
            print(m) if logger is None else logger.error(m)
            exit()

    with open(fname, 'w') as fileout:
        fileout.writelines(lines_ffbonded)

    m = "\t\t Some bonded terms related with the methylcarboxilic group" \
        " have been added to ffbonded.itp file by similarity.\n"
    m += "\t\t Urey-Bradley: CC202-CC321-OC301 (Cm2-Cm1-O6) similar to CC202-CC311-OC301\n"
    m += "\t\t Torsion     : CC3163-CC321-OC301-CC321 similar to  CC3163-CC321-OC301-CTL\n "
    m += "\t\t Torsion     : CC2O2-CC321-OC301-CC321  similar to  CC2O2-CC311-OC301-CC331\n "
    m += "\t\t Torsion     : OC301-CC321-CC2O2-OC2D2  similar to  OC301-CC301-CC2O2-OC2D2\n "

    print(m) if logger is None else logger.info(m)




# =============================================================================
def main_app(version=None):

    # Init logger =============================================================
    logger = init_logger("Output", fileoutput="InfoTypeCB.log",
                         append=False, inscreen=True)
    print_header_type_cb(version, logger)
    opts = parse_arguments_typing_carb()
    starttime = datetime.datetime.now()
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    logger.info("\t\tStarting: \t {}\n".format(now))

    # Print command info
    command_info(opts, logger=logger)

    # Copy and add some parameters neccesary for the carboxymethyl group
    edit_charm36_ff(logger=logger)


    u = Mda.Universe(opts.pdbfile, opts.pdbfile)
    create_gro_top_charmm36(u, opts.boxdimensions, logger)

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
