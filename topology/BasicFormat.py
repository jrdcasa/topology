import MDAnalysis as Mda
import networkx as nx
from collections import defaultdict
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


class BasicFormat(object):

    # ====================================================================
    def __init__(self, mol_mda, logger=None):

        self._mol_mda = mol_mda
        self._logger = logger
        self._graph_mol = None
        self._ff = None
        self._imol = 0
        self._read_charges = list()
        self._read_typeatoms = list()
        self._read_dict_charges = dict()
        self._index_to_type_dict = dict()
        self._formulas_name_list = list()
        self._chiral_centers = list()
        self._ions_list = ["NA", "LI"]
        self._kb_kjmol = 0.0083144621  # kJ/molK

        if not isinstance(mol_mda, Mda.core.universe.Universe):
            m = "\n\t\t A Universe object must be passed to GromacsFormat class.\n"
            m1 = "\t\t" + "*" * len(m)
            m += "\t\t Received type: {0}\n".format(type(mol_mda))
            m += "\t\t Expected type: MDAnalysis.core.universe.Universe\n"
            m = "\n" + m1 + m + m1
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        self._graph_mol = None
        self._graph_unique_mol = None
        self._graph_component_nodes = None
        self._nmols = None
        self._nmols_unique = None

    # =============================================================================
    def _number_molecules_graph(self):

        graph_mol = nx.Graph()
        graph_component_nodes = defaultdict()

        # Add nodes and edges
        for idx_bond, ibond in enumerate(self._mol_mda.bonds):
            idx1 = ibond.indices[0]
            idx2 = ibond.indices[1]
            graph_mol.add_node(idx1)
            graph_mol.add_node(idx2)
            graph_mol.add_edge(idx1, idx2)

        # Set attributtes for each node
        attrs = defaultdict(dict)

        for inode in graph_mol.nodes:
            attrs[inode] = {"resid": None, "name": None, "resname": None, "link": [],
                            "element": None, "patch_atom": False}

        for inode in graph_mol.nodes:
            resid = self._mol_mda.atoms[inode].resid
            name = self._mol_mda.atoms[inode].name
            resname = self._mol_mda.atoms[inode].resname
            element = self._mol_mda.atoms[inode].element
            attrs[inode]["resid"] = resid
            attrs[inode]["name"] = name
            attrs[inode]["smiles"] = ""
            attrs[inode]["resname"] = resname
            attrs[inode]["element"] = element
            attrs[inode]["type"] = ""

        nx.set_node_attributes(graph_mol, attrs)

        # Number of disconnected graphs
        nmols = nx.number_connected_components(graph_mol)
        for idx, imol in enumerate(nx.connected_components(graph_mol)):
            graph_component_nodes[idx] = sorted(imol)

        # Isomorphic components
        graph_unique_mol = defaultdict(list)
        components = [graph_mol.subgraph(c).copy() for c in nx.connected_components(graph_mol)]

        for subgraph in components:
            h = weisfeiler_lehman_graph_hash(subgraph, node_attr="name")
            graph_unique_mol[h].append(subgraph)
        nmols_unique = len(graph_unique_mol)

        return graph_mol, components, graph_unique_mol, graph_component_nodes, nmols, nmols_unique

    # ====================================================================
    def _molecular_names(self):

        for ikey, igraph in self._graph_unique_mol.items():
            G_tmp = igraph[0]

            # 1. Convert to RDKit molecule
            mol = Chem.RWMol()
            node_to_idx = {}

            # 2. Add atoms
            for n, data in G_tmp.nodes(data=True):
                atom = Chem.Atom(data["element"])
                idx = mol.AddAtom(atom)
                node_to_idx[n] = idx

            # 3. Add bonds (single bonds by default)
            for i, j in G_tmp.edges():
                mol.AddBond(node_to_idx[i], node_to_idx[j], Chem.BondType.SINGLE)

            # 4. Sanitize molecule
            mol = mol.GetMol()
            Chem.SanitizeMol(mol)

            # Molecular formula
            formula = rdMolDescriptors.CalcMolFormula(mol)
            smiles = Chem.MolToSmiles(mol, canonical=True)
            # print("Formula:", formula)
            # print("Canonical SMILES:", smiles)

            # Insert formula in the list
            tmp_formula = formula
            if formula not in self._formulas_name_list:
                self._formulas_name_list.append(formula)
            else:
                sentinel = 0
                while True:
                    formula = tmp_formula + "_{0:1d}".format(sentinel)
                    if formula not in self._formulas_name_list:
                        self._formulas_name_list.append(formula)
                        break
                    sentinel += 1

            for item in igraph:
                item.name = formula
                item.smiles = smiles

    # =============================================================================
    def read_charges_types(self, chargefile):

        """

        Args:
            chargefile:

        Returns:

        """

        # Read the charge file
        with open(chargefile, 'r') as fchg:
            lines = fchg.readlines()
            for iline in lines:
                if iline.find("#") != -1:
                    continue
                try:
                    tokens = iline.split()
                    self._read_typeatoms.append(str(tokens[2]))
                    self._read_charges.append(float(tokens[3]))
                    self._read_dict_charges[str(tokens[2])] = float(tokens[3])
                    idx = int(tokens[0]) - 1
                    self._index_to_type_dict[idx] = str(tokens[2])
                except (IndexError, ValueError):
                    continue


    # =============================================================================
    def _guess_impropers(self):
        """
        Guess impropers as LigParGen server

        ligpargen/ligpargen/tools/geometry.py
        (https://github.com/Isra3l/ligpargen.git)

        Returns:

        """

        impropers = []
        for iatom in self._mol_mda.atoms:

            neigh = iatom.bonded_atoms

            if len(neigh) == 3:

                atomIndex = iatom.index
                atomNIndex = [atomA.index for atomA in neigh]

                imp = [atomNIndex[2],  atomIndex, atomNIndex[0], atomNIndex[1] ]
                impropers.append(imp)

        return impropers

    # =============================================================================
    @staticmethod
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
