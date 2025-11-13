import numpy as np
import networkx as nx
from collections import deque, defaultdict
from topology.atomic_data import atomic_number
from functools import cmp_to_key


# ===============================================================================================
def max_flat(tree):
    """
    Recursively find the maximum atomic number in a substituent tree.
    Works even if children are integers (leaf nodes).
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
    """
    Build CIP substituent sequence from 'start' atom outward, avoiding 'center'.
    Returns list of lists: each entry is a level of atomic numbers.
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
    """
    Sort neighbors of a stereocenter by CIP priority using recursive trees
    G : NetworkX graph
    center : int - index of stereocenter atom
    max_depth : int - maximum recursion depth
    Returns: list of neighbor indices sorted by descending CIP priority
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
    """
    Identify tetrahedral stereocenters in an MDAnalysis universe
    Criteria:
        - Atom has exactly 4 bonded neighbors
        - Neighbors are not all identical (at least 2 different elements)
    Returns: list of Atom objects that are potential stereocenters
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
    """
    Build a NetworkX graph of the molecule from MDAnalysis universe
    Nodes: atom indices with 'element' attribute
    Edges: bonds with 'order' attribute (default 1)
    Returns: NetworkX Graph
    """

    G = nx.Graph()
    for atom in u.atoms:
        G.add_node(atom.index, element=atom.element)
    for bond in u.bonds:
        G.add_edge(bond.atoms[0].index, bond.atoms[1].index, order=1)  # default single bond
    return G


# ===============================================================================================
def minimum_image_vector(pos_i, pos_j, box):
    """
    Compute displacement vector between two points considering periodic boundary conditions (PBC)
    pos_i, pos_j : np.array([x, y, z])
    box : np.array([Lx, Ly, Lz])
    Returns: PBC-corrected vector from pos_i to pos_j
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
    """
    Compute R or S configuration for a stereocenter
    center : np.array([3,]) - coordinates of stereocenter
    neighbors_coords : dict {label: np.array([3,])} - coordinates of 4 substituents
    priority_labels : list of 4 labels in descending CIP priority
    box : np.array([Lx, Ly, Lz]) - box dimensions for PBC
    Returns: "R" or "S"
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
    """
    Assign R/S configuration to all tetrahedral stereocenters in a MDAnalysis Universe
    Steps:
        - Convert universe to NetworkX graph
        - Detect tetrahedral stereocenters
        - Compute CIP priorities with NetworkX
        - Compute R/S using PBC-aware scalar triple product
    Returns: list of tuples (atom index, "R"/"S")
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
    """
    Determine R/S chirality of all stereocenters.

    Pipeline for manual R/S assignment from an MDAnalysis Universe,
    using NetworkX for CIP ranking, handling ties and multiple bonds,
    and computing PBC-aware R/S configuration.

    u_mda: MDAnalysis Universe
    logger: logging function

    Returns: "R" or "S"
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



