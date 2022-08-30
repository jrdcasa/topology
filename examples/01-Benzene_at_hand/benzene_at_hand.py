import topology as top


"""

                H6
                |
                C0
              /   \
       H7 - C1    C5 - H11
             |     |
       H8 - C2    C4 - H10
              \   /
               C3
               |
               H9
"""

# Number of atoms of the molecule
natoms = 12  # 6C and 6H
# List of all bonds
bond_list = [[0, 1], [0, 2], [0, 5], [0, 6], [1, 2],
             [1, 7], [2, 3], [2, 8], [3, 4], [3, 9],
             [4, 10], [5, 11]]

# Create a topology
t = top.Topology(natoms=natoms, listbonds=bond_list)

# Assign elements
elements = ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
t.set_elements(elements)

# Assign charges to the atoms
chargelist = [0.0]*natoms
t.set_charge_mol(chargelist)

# Assign bond orders
t.assign_bond_orders()