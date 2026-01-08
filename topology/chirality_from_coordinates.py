import numpy as np
import networkx as nx
from collections import deque, defaultdict
from topology.atomic_data import atomic_number
from functools import cmp_to_key


# ===============================================================================================
def max_flat(tree):
    """Recursively find the maximum atomic number in a substituent tree.

    Works even if children are integers (leaf nodes).

    Args:
        tree (list): Nested list representing a substituent tree. Leaves can be integers.

    Returns:
        int: Maximum atomic number found in the tree.
    """
    if not tree:
        return 0
    max_val = tree[0]
    for c in tree[1:]:
        if isinstance(c, list):
            max_val = max(max_val, max_flat(c))
        else:  # c is int (leaf)
            max_val = max(max_val, c)
    return max_val


# ===============================================================================================
def bfs_sequence(G, start, center, max_depth=6):
    """Build a CIP substituent sequence from a starting atom outward, avoiding the central atom.

    Args:
        G (nx.Graph): Molecular graph.
        start (int): Index of the starting atom.
        center (int): Index of the stereocenter atom to avoid.
        max_depth (int, optional): Maximum recursion depth. Defaults to 6.

    Returns:
        list[tuple]: List of tuples, each representing a level of atomic numbers in descending order.
    """
    visited = {center, start}
    q = deque([(start, 0)])
    levels = {}

    while q:
        node, depth = q.popleft()
        if depth >= max_depth:
            continue

        levels.setdefault(depth, []).append(atomic_number[G.nodes[node]['element']])

        for neigh in G.neighbors(node):
            if neigh not in visited:
                visited.add(neigh)
                q.append((neigh, depth + 1))

    # sort each level descending for CIP
    seq = [tuple(sorted(levels[d], reverse=True)) for d in sorted(levels)]
    return seq


# ===============================================================================================
def compare_seq(seq1, seq2):

    """Compare two CIP sequences for priority ranking.

    Args:
        seq1 (list[tuple]): First sequence.
        seq2 (list[tuple]): Second sequence.

    Returns:
        int: 1 if seq1 > seq2, -1 if seq1 < seq2, 0 if equal.
    """



    for a, b in zip(seq1, seq2):
        # extend shorter sequence with zeros if needed
        la, lb = list(a), list(b)
        while len(la) < len(lb): la.append(0)
        while len(lb) < len(la): lb.append(0)

        if la > lb:
            return 1
        elif la < lb:
            return -1
    # if one is longer
    if len(seq1) > len(seq2):
        return 1
    elif len(seq1) < len(seq2):
        return -1
    return 0


# ===============================================================================================
def cip_sort(G, center, max_depth=6):
    """Sort neighbors of a stereocenter by CIP priority using recursive trees.

    Args:
        G (nx.Graph): Molecular graph.
        center (int): Index of the stereocenter atom.
        max_depth (int, optional): Maximum recursion depth for BFS. Defaults to 6.

    Returns:
        list[int] or None: List of neighbor indices sorted by descending CIP priority,
                           or None if priorities are not unique.
    """

    def comp(a, b):
        return -compare_seq(seqs[a], seqs[b])  # negate so higher priority comes first

    neighbors = list(G.neighbors(center))

    seqs = defaultdict()
    for n in neighbors:
        seqs[n] = bfs_sequence(G, n, center, max_depth)

    # Check that all neighbors are differents
    values = list(seqs.values())
    all_unique = len(values) == len(set(map(tuple, values)))

    if all_unique:
        sorted_neighbors = sorted(neighbors, key=cmp_to_key(comp))
        return sorted_neighbors
    else:
        return None


# ===============================================================================================
def find_tetrahedral_centers(u):
    """Identify tetrahedral stereocenters in an MDAnalysis Universe.

    Criteria:
        - Atom has exactly 4 bonded neighbors.
        - Neighbors are not all identical (at least 2 different elements).

    Args:
        u (MDAnalysis.core.universe.Universe): MDAnalysis Universe object.

    Returns:
        list[MDAnalysis.core.groups.Atom]: List of potential stereocenter atoms.
    """

    centers = []
    for atom in u.atoms:
        neighbors = list(atom.bonded_atoms)
        if len(neighbors) != 4:  # tetrahedral
            continue
        elements = [a.element for a in neighbors]
        if len(set(elements)) >= 2:  # at least 2 different substituents
            centers.append(atom)
    return centers


# ===============================================================================================
def mol_to_graph(u):
    """Build a NetworkX graph from an MDAnalysis Universe.

    Nodes are atoms with element attribute. Edges are bonds (default order 1).

    Args:
        u (MDAnalysis.core.universe.Universe): MDAnalysis Universe object.

    Returns:
        nx.Graph: Molecular graph.
    """

    G = nx.Graph()
    for atom in u.atoms:
        G.add_node(atom.index, element=atom.element)
    for bond in u.bonds:
        G.add_edge(bond.atoms[0].index, bond.atoms[1].index, order=1)  # default single bond
    return G


# ===============================================================================================
def minimum_image_vector(pos_i, pos_j, box):
    """Compute displacement vector considering periodic boundary conditions (PBC).

    Args:
        pos_i (np.ndarray): Position of first point (3,).
        pos_j (np.ndarray): Position of second point (3,).
        box (np.ndarray): Simulation box dimensions (3,).

    Returns:
        np.ndarray: PBC-corrected vector from pos_i to pos_j.
    """
    delta = pos_j - pos_i
    for k in range(3):
        if delta[k] > 0.5 * box[k]:
            delta[k] -= box[k]
        elif delta[k] < -0.5 * box[k]:
            delta[k] += box[k]
    return delta


# ===============================================================================================
def chirality_r_s(center, neighbors_coords, priority_labels, box):
    """Compute R or S configuration for a stereocenter using PBC-aware scalar triple product.

    Args:
        center (np.ndarray): Coordinates of stereocenter (3,).
        neighbors_coords (dict[str, np.ndarray]): Mapping of neighbor label to coordinates (3,).
        priority_labels (list[str]): Labels of 4 substituents in descending CIP priority.
        box (np.ndarray): Box dimensions for PBC (3,).

    Returns:
        str: "R" or "S" configuration.
    """

    p1, p2, p3, p4 = priority_labels
    # Compute vectors using minimum image convention
    v1 = minimum_image_vector(neighbors_coords[p1], center, box)
    v2 = minimum_image_vector(neighbors_coords[p2], center, box)
    v3 = minimum_image_vector(neighbors_coords[p3], center, box)
    # Scalar triple product
    s = np.dot(np.cross(v1, v2), v3)

    # print("cJ1:", s, p1, p2, p3, p4)

    return "R" if s > 0 else "S"


# ===============================================================================================
def assign_r_s(u):
    """Assign R/S configuration to all tetrahedral stereocenters in a MDAnalysis Universe.

    Steps:
        - Convert universe to NetworkX graph.
        - Detect tetrahedral stereocenters.
        - Compute CIP priorities with NetworkX.
        - Compute R/S using PBC-aware scalar triple product.

    Args:
        u (MDAnalysis.core.universe.Universe): MDAnalysis Universe object.

    Returns:
        list[tuple]: List of tuples (atom index, "R"/"S", sorted neighbor indices).
    """

    G = mol_to_graph(u)
    box = u.dimensions[:3]
    u.trajectory[0]  # pick first frame
    results = []
    centers = find_tetrahedral_centers(u)
    for center_atom in centers:
        sorted_neighbor_indices = cip_sort(G, center_atom.index)
        if sorted_neighbor_indices is None:
            continue
        # convert indices to Atom objects
        neighbors_atoms = [u.atoms[i] for i in sorted_neighbor_indices]
        neighbors_coords = {a.name: a.position for a in neighbors_atoms}
        priority_labels = [a.name for a in neighbors_atoms]
        r_s = chirality_r_s(center_atom.position, neighbors_coords, priority_labels, box)
        sorted_neighbor_idx = list(map(lambda x: x - 1, sorted_neighbor_indices))
        results.append((center_atom.index, r_s, sorted_neighbor_idx))
    return results


# ===============================================================================================
def chirality_from_coords(u_mda, logger=None):
    """Determine R/S chirality of all stereocenters from an MDAnalysis Universe.

    Args:
        u_mda (MDAnalysis.core.universe.Universe): MDAnalysis Universe object.
        logger (logging.Logger, optional): Logger to print information. Defaults to None.

    Returns:
        list[tuple]: List of tuples (atom index, "R"/"S", neighbor indices).
    """

    chirality_list = assign_r_s(u_mda)

    if len(chirality_list) == 0:
        m = "\n\t\t Not Chiral centers \n"
    else:
        m = "\n\t\t Chiral centers (Index/Chiral): \n"
        m += "\t\t (R: Right-handed and S:Left-handed forms) \n"
        for item in chirality_list:
            m += "\t\t\t {0:d}/{1:s} atoms =".format(item[0], item[1])
            m += "  ".join(map(str, item[2]))
    print(m) if logger is None else logger.info(m)

    return chirality_list



