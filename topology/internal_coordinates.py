"""
Miscellaneous mathematical functions for molecular simulations.
"""

import math
import numpy as np
import random
import psutil
from scipy.spatial import distance_matrix


# #############################################################################
# noinspection PyUnboundLocalVariable
def distance_array(ref, conf, openmp=False):

    """
    Calculate pairwise distances and distance-vector components between
    all points in two coordinate arrays. Periodic boundary conditions
    are not applied.

    Args:
        ref (ndarray): Reference coordinates of shape (npoints1, dim),
            where dim is 2 or 3.
        conf (ndarray): Target coordinates of shape (npoints2, dim),
            where dim is 2 or 3.
        openmp (bool, optional): If True, use the parallel OpenMP C/Cython
            implementation. Defaults to False.

    Returns:
        tuple:
            dist (ndarray): Distances of shape (npoints1, npoints2).
            rijx (ndarray): X-components of distance vectors.
            rijy (ndarray): Y-components of distance vectors.
            rijz (ndarray): Z-components of distance vectors.
    """

    if openmp:
        from topology.ext_libc.c_distances_openmp import calc_distance_array_openmp
        # from c_distances_openmp import calc_distance_array_openmp
    else:
        from topology.ext_libc.c_distances import calc_distance_array
        # from c_distances import calc_distance_array

    rows = ref.shape[0]
    cols = conf.shape[0]

    # Check memory
    vmem_available = psutil.virtual_memory()[1]  # bytes
    # 4 arrays of rows*cols --> float64
    mem_req = 4*rows*cols*8  # bytes
    if 4 * rows * cols * 8 > vmem_available:
        print("There is not sufficient memory to calculate distances...")
        print("Internal_coordinates.py --> distance_array fucntion.")
        print("Memory available: {0:.1f} Mb".format(vmem_available/1024/1024))
        print("Memory required : {0:.1f} Mb".format(mem_req/1024/1024))
        print("Aborting program...")
        exit()

    dist = np.zeros([rows, cols])
    rijx = np.zeros([rows, cols])
    rijy = np.zeros([rows, cols])
    rijz = np.zeros([rows, cols])
    if openmp:
        calc_distance_array_openmp(ref, conf, dist, rijx, rijy, rijz)
    else:
        calc_distance_array(ref, conf, dist, rijx, rijy, rijz)

    return dist, rijx, rijy, rijz


# #############################################################################
# noinspection PyUnboundLocalVariable
def distance_diagonal_array(ref, conf, openmp=False):

    """
    Compute distances between corresponding points (i -> i) in two arrays.
    Only diagonal distances are evaluated. Periodic boundary conditions
    are not applied.

    Args:
        ref (ndarray): Reference coordinates of shape (npoints, dim).
        conf (ndarray): Target coordinates of shape (npoints, dim).
        openmp (bool, optional): If True, use the parallel OpenMP version.

    Returns:
        tuple:
            dist (ndarray): Diagonal distances of shape (npoints,).
            rijx (ndarray): X-components of distance vectors.
            rijy (ndarray): Y-components.
            rijz (ndarray): Z-components.
    """
    if openmp:
        from topology.ext_libc.c_distances_openmp import calc_distance_diagonal_openmp
        # from c_distances_openmp import calc_distance_array_openmp
    else:
        from topology.ext_libc.c_distances import calc_distance_diagonal
        # from c_distances import calc_distance_array

    rows = ref.shape[0]
    dist = np.zeros([rows])
    rijx = np.zeros([rows])
    rijy = np.zeros([rows])
    rijz = np.zeros([rows])
    if openmp:
        calc_distance_diagonal_openmp(ref, conf, dist, rijx, rijy, rijz)
    else:
        calc_distance_diagonal(ref, conf, dist, rijx, rijy, rijz)

    # calc_distance_diagonal_numba(ref, conf, dist, rijx, rijy, rijz)

    return dist, rijx, rijy, rijz


# #############################################################################
def distance_array_numpypython(ref, conf):

    """
    Compute all pairwise distances and distance-vector components using
    pure NumPy operations.

    Args:
        ref (ndarray): Coordinates of shape (npoints1, 3).
        conf (ndarray): Coordinates of shape (npoints2, 3).

    Returns:
        tuple:
            dist (ndarray): Distances of shape (npoints1, npoints2).
            rijx (ndarray): X-components of distance vectors.
            rijy (ndarray): Y-components.
            rijz (ndarray): Z-components.
    """

    rows = ref.shape[0]     # m
    cols = conf.shape[0]    # n
    rijx = np.zeros([rows, cols])
    rijy = np.zeros([rows, cols])
    rijz = np.zeros([rows, cols])

    irow = 0
    for ipoint in ref:
        rijx[irow, :] = conf[:, 0] - ipoint[0]
        rijy[irow, :] = conf[:, 1] - ipoint[1]
        rijz[irow, :] = conf[:, 2] - ipoint[2]
        irow += 1

    # dist = np.zeros([rows,cols])
    dist = distance_matrix(ref, conf)
    return dist, rijx, rijy, rijz


# #############################################################################
def distance_array_purepython(ref, conf):

    """
    Compute all pairwise distances and vector components using pure
    Python loops. This is slower than NumPy or compiled versions.

    Args:
        ref (ndarray): Coordinates of shape (npoints1, 3).
        conf (ndarray): Coordinates of shape (npoints2, 3).

    Returns:
        tuple:
            dist (ndarray): Distances of shape (npoints1, npoints2).
            rijx (ndarray): X-components.
            rijy (ndarray): Y-components.
            rijz (ndarray): Z-components.
    """

    rows = ref.shape[0]
    cols = conf.shape[0]
    dist = np.zeros([rows, cols])
    rijx = np.zeros([rows, cols])
    rijy = np.zeros([rows, cols])
    rijz = np.zeros([rows, cols])
    irow = 0
    for ipoint in ref:
        icol = 0
        for jpoint in conf:
            rijx[irow, icol] = jpoint[0] - ipoint[0]
            rijy[irow, icol] = jpoint[1] - ipoint[1]
            rijz[irow, icol] = jpoint[2] - ipoint[2]
            dist[irow, icol] = math.sqrt(rijx[irow, icol]*rijx[irow, icol] +
                                         rijy[irow, icol]*rijy[irow, icol] +
                                         rijz[irow, icol]*rijz[irow, icol])
            icol += 1
        irow += 1

    return dist, rijx, rijy, rijz

# #############################################################################
# jit(nopython=True, nogil=True, cache=True)
# def calc_distance_array_numba(ref, conf):
#
#     rows = ref.shape[0]     # m
#     cols = conf.shape[0]    # n
#     rijx = np.zeros([rows,cols])
#     rijy = np.zeros([rows,cols])
#     rijz = np.zeros([rows,cols])
#     for irow in numba.prange(rows):
#         for icol in numba.prange(cols):
#             rijx[irow, icol] = conf[icol, 0] - ref[irow, 0]
#             rijy[irow, icol] = conf[icol, 1] - ref[irow, 1]
#             rijz[irow, icol] = conf[icol, 2] - ref[irow, 2]
#     dist = distance_matrix(ref, conf)
#
#
#     return dist, rijx, rijy, rijz

# #############################################################################
# jit(nopython=True, cache=True, parallel=True)
# def calc_distance_diagonal_numba(ref, conf, dist, rijx, rijy, rijz):
#
#     row = ref.shape[0]
#
#     for r in numba.prange(row):
#         rijx[r] = conf[r, 0] - ref[r, 0]
#         rijy[r] = conf[r, 1] - ref[r, 1]
#         rijz[r] = conf[r, 2] - ref[r, 2]
#         dist[r] = math.sqrt((rijx[r] * rijx[r]) +
#                              (rijy[r] * rijy[r]) +
#                              (rijz[r] * rijz[r]) )
#
#     return None


# #############################################################################
def bend_angle_purepython(c1, c2, c3, radians=False):

    """
    Compute the bend angle defined by three atoms (1–2–3).

    Args:
        c1 (ndarray): Coordinates of atom 1, shape (3,).
        c2 (ndarray): Coordinates of atom 2, shape (3,).
        c3 (ndarray): Coordinates of atom 3, shape (3,).
        radians (bool, optional): If True, return angle in radians.
            Otherwise return degrees.

    Returns:
        float: Bend angle.
    """

    p12 = distance_array(c1.reshape(1, 3), c2.reshape(1, 3))
    p23 = distance_array(c2.reshape(1, 3), c3.reshape(1, 3))
    p13 = distance_array(c1.reshape(1, 3), c3.reshape(1, 3))

    theta = math.acos((pow(p12[0], 2.)+pow(p23[0], 2.)-pow(p13[0], 2.))/(2.*p12[0]*p23[0]))
    if not radians:
        theta = theta * 180. / math.pi
    return theta


# #############################################################################
def dihedral_angle_purepython(ci, cj, ck, cl, radians=False):

    """
    Compute the proper or improper dihedral angle defined by four atoms.

    Args:
        ci (ndarray): Atom i coordinates.
        cj (ndarray): Atom j coordinates.
        ck (ndarray): Atom k coordinates.
        cl (ndarray): Atom l coordinates.
        radians (bool, optional): If True, return radians; else degrees.

    Returns:
        float: Dihedral angle.
    """

    rij = ci - cj
    rjk = cj - ck
    rlk = cl - ck
    m = np.cross(rij, rjk)
    n = np.cross(rlk, rjk)
    m_norm = np.linalg.norm(m)
    n_norm = np.linalg.norm(n)
    cos_ijkl = np.dot(m, n)/(m_norm*n_norm)
    sin_ijkl = np.dot(n, rij)*np.linalg.norm(rjk)/(m_norm*n_norm)
    phi = -np.arctan(sin_ijkl/cos_ijkl)
    if not radians:
        phi = phi * 180. / math.pi
    return phi


# #############################################################################
def dihedral_angle_purepython_vmd(ci, cj, ck, cl, radians=False):

    """
    Compute a dihedral or improper angle, using a formulation following
    VMD conventions. Periodic boundary conditions are not applied.

    The function follows the geometry illustrated in:

        ../../figures/dih_imp.png

    Args:
        ci (ndarray): Coordinates of atom i.
        cj (ndarray): Coordinates of atom j.
        ck (ndarray): Coordinates of atom k.
        cl (ndarray): Coordinates of atom l.
        radians (bool): If True, return radians; otherwise degrees.

    Returns:
        float: Dihedral angle.
    """

    # 1st bond
    delx1 = ci[0] - cj[0]
    dely1 = ci[1] - cj[1]
    delz1 = ci[2] - cj[2]

    # 2nd bond
    delx2 = ck[0] - cj[0]
    dely2 = ck[1] - cj[1]
    delz2 = ck[2] - cj[2]

    delx2m = -delx2
    dely2m = -dely2
    delz2m = -delz2

    # 3rd bond
    delx3 = cl[0] - ck[0]
    dely3 = cl[1] - ck[1]
    delz3 = cl[2] - ck[2]

    # c,s calculation
    ax = dely1*delz2m - delz1*dely2m
    ay = delz1*delx2m - delx1*delz2m
    az = delx1*dely2m - dely1*delx2m
    bx = dely3*delz2m - delz3*dely2m
    by = delz3*delx2m - delx3*delz2m
    bz = delx3*dely2m - dely3*delx2m

    rasq = ax*ax + ay*ay + az*az
    rbsq = bx*bx + by*by + bz*bz
    rgsq = delx2m*delx2m + dely2m*dely2m + delz2m*delz2m
    rg = np.sqrt(rgsq)

    rginv = ra2inv = rb2inv = 0.0
    if rg > 0:
        rginv = 1.0/rg
    if rasq > 0:
        ra2inv = 1.0/rasq
    if rbsq > 0:
        rb2inv = 1.0/rbsq;
    rabinv = np.sqrt(ra2inv*rb2inv)

    c = (ax*bx + ay*by + az*bz)*rabinv
    s = rg*rabinv*(ax*delx3+ay*dely3+az*delz3)

    if c > 1.0:
        c = 1.0
    if c < -1.0:
        c = -1.0

    phi = np.arctan2(s, c)
    if not radians:
        phi = phi * 180. / math.pi

    return phi


# #############################################################################
def cos_angle_purepython(ci, cj, ck, cl, radians=False):

    """
    Compute the cosine angle measure used for bond perception in:

        Zhang et al., J. Cheminformatics 2012, 4:26.

    Args:
        ci (ndarray): Atom i coordinates.
        cj (ndarray): Atom j coordinates.
        ck (ndarray): Atom k coordinates.
        cl (ndarray): Atom l coordinates.
        radians (bool): If True, return radians; else degrees.

    Returns:
        float: Angle between particles.
    """

    rji = cj - ci
    rki = ck - ci
    rli = cl - ci
    m = np.cross(rji, rki)

    m_norm = np.linalg.norm(m)
    rli_norm = np.linalg.norm(rli)

    cos_ijkl = np.dot(m, rli)/(m_norm*rli_norm)

    phi = np.arccos(cos_ijkl)
    if not radians:
        phi = phi * 180. / math.pi
    return phi


# #############################################################################
def generate_random_euler_angles(seed=None):

    """
    Generate uniformly distributed Euler angles (α, β, γ).

    Args:
        seed (int, optional): Seed for the random number generator.

    Returns:
        list: Three Euler angles in radians.
    """

    if seed is None:
        random.seed()
    else:
        random.seed(seed)

    alpha = 2.0 * math.pi * random.random()  # Radians
    c_beta = 1.0 - 2.0 * random.random()
    beta = math.acos(c_beta)                # Radians
    gamma = 2.0 * math.pi * random.random()  # Radians

    return [alpha, beta, gamma]


# #############################################################################
def euler_rotation_matrix(euler):

    """
    Construct a rotation matrix from Euler angles using the convention
    described in:

        Arfken, Weber & Harris, "Mathematical Methods for Physicists",
        7th edition, pp. 140–142.

    Args:
        euler (array-like): Euler angles [α, β, γ] in radians.

    Returns:
        ndarray: 3×3 rotation matrix.
    """

    a = euler[0]
    b = euler[1]
    g = euler[2]
    ca = math.cos(a)
    cb = math.cos(b)
    cg = math.cos(g)
    sa = math.sin(a)
    sb = math.sin(b)
    sg = math.sin(g)

    S = np.array([[(cg*cb*ca) - (sg*sa), (cg*cb*sa) + (sg*ca), -cg*sb],
                  [(-sg*cb*ca) - (cg*sa), (-sg*cb*sa) + (cg*ca),  sg*sb],
                  [sb*ca, sb*sa, cb]])

    return S


# #############################################################################
def center_of_geom_purepython(coords):

    """
    Compute the center of geometry (average position). Periodic
    boundary conditions are not considered. Pure Python implementation.

    Args:
        coords (list or ndarray): Atomic coordinates (unwrapped).

    Returns:
        ndarray: Geometry center.
    """

    if isinstance(coords, list):
        _coords = np.asarray(coords)
    elif isinstance(coords, np.ndarray):
        _coords = coords
    else:
        return None

    tmp = np.zeros(3)
    natoms = _coords.shape[0]
    for iatom in range(natoms):
        tmp += _coords[iatom, :]

    cog = tmp/natoms
    return cog


# #############################################################################
def center_of_geom(coords):

    """
    Compute the center of geometry using a compiled C routine.
    Periodic boundary conditions are not applied.

    Args:
        coords (list or ndarray): Atomic coordinates (unwrapped).

    Returns:
        ndarray: Geometry center.
    """

    from topology.ext_libc.c_unwrap_openmp import calc_center_of_geom

    if isinstance(coords, list):
        _coords = np.asarray(coords, dtype=np.float32)
    elif isinstance(coords, np.ndarray):
        _coords = np.float32(coords)
    else:
        return None

    cog = np.zeros(3, dtype=np.float32)

    calc_center_of_geom(_coords, cog)

    return cog


# #############################################################################
def center_of_mass_purepython(coords, mass):

    """
    Compute the center of mass using pure Python. Periodic boundary
    conditions are not applied.

    Args:
        coords (list or ndarray): Atomic coordinates (unwrapped).
        mass (list or ndarray): Atomic masses.

    Returns:
        ndarray: Center of mass.
    """

    if isinstance(coords, list):
        _coords = np.asarray(coords)
    elif isinstance(coords, np.ndarray):
        _coords = coords
    else:
        return None

    if isinstance(mass, list):
        _mass = np.asarray(mass)
    elif isinstance(mass, np.ndarray):
        _mass = mass
    else:
        return None

    if _coords.shape[0] != mass.shape[0]:
        return None

    mtotal = np.sum(_mass)
    tmp = np.zeros(3)
    natoms = _coords.shape[0]
    for iatom in range(natoms):
        tmp += _coords[iatom, :] * _mass[iatom]

    com = tmp/mtotal
    return com


# #############################################################################
def center_of_mass(coords, mass):

    """
    Compute the center of mass using a C routine. Periodic boundary
    conditions are not applied.

    Args:
        coords (list or ndarray): Coordinates (unwrapped).
        mass (list or ndarray): Masses.

    Returns:
        ndarray: Center of mass.
    """

    from topology.ext_libc.c_unwrap_openmp import calc_center_of_mass

    if isinstance(coords, list):
        _coords = np.asarray(coords)
    elif isinstance(coords, np.ndarray):
        _coords = np.float32(coords)
    else:
        return None

    if isinstance(mass, list):
        _mass = np.asarray(mass, dtype=np.float32)
    elif isinstance(mass, np.ndarray):
        _mass = np.float32(mass)
    else:
        return None

    if _coords.shape[0] != mass.shape[0]:
        return None

    com = np.zeros(3, dtype=np.float32)

    calc_center_of_mass(_coords, _mass, com)

    return com


# #############################################################################
def unwrap_purepython(coords, topology, box_dimensions, iframe=0, write_xyz=False):

    """
    Unwrap coordinates using pure Python. Periodic boundary conditions
    are resolved following molecular connectivity.

    Args:
        coords (ndarray): Wrapped coordinates, shape (natoms, 3).
        topology (Topology): Topology object providing connectivity.
        box_dimensions (ndarray): Box lengths [Lx, Ly, Lz].
        iframe (int, optional): Frame index for XYZ output.
        write_xyz (bool, optional): Write XYZ output if True.

    Returns:
        ndarray: Unwrapped coordinates.
    """

    natoms = coords.shape[0]
    nmols = topology.get_nmols()
    cu = np.zeros((natoms, 3))
    # Store if a node has been visited, and therefore its coordinates
    # are already unwrapped.
    isvisited = np.zeros(natoms, dtype=int)

    # Each molecule can be unwrapped in parallel.
    for imol in nmols:
        # First atom of the molecule
        iatom = imol[0]
        cu[iatom, :] = coords[iatom, :]
        isvisited[iatom] = 1

        # loop over the atoms in the current molecule
        for ind in range(1, len(imol)):
            iatom = imol[ind]
            prev = -1
            for jatom in topology.get_neighbours(imol[ind]):
                if isvisited[jatom] == 1:
                    prev = jatom
            isvisited[iatom] = 1
            d = coords[iatom, :] - cu[prev, :]
            d -= box_dimensions[0:3]*np.round(d/box_dimensions[0:3])
            cu[iatom, :] = cu[prev, :] + d

    if write_xyz:
        filename = "movie_unwrap_{0:06d}.xyz".format(iframe)

        with open(filename, 'w') as f:
            f.writelines("{0:d}\n".format(natoms))
            f.writelines("Frame wrapped: {}\n".format(iframe))
            for iatom in range(natoms):
                line = "{0:s} {1:.4f} {2:.4f} {3:.4f}\n".\
                    format('C', coords[iatom, 0], coords[iatom, 1], coords[iatom, 2])
                f.writelines(line)

            f.writelines("{0:d}\n".format(natoms))
            f.writelines("Frame unwrapped: {}\n".format(iframe))
            for iatom in range(natoms):
                line = "{0:s} {1:.4f} {2:.4f} {3:.4f}\n".\
                    format('C', cu[iatom, 0], cu[iatom, 1], cu[iatom, 2])
                f.writelines(line)

    return cu


# #############################################################################
def unwrap(coords, nmols_array, l_neigh_array, box_dimensions, iframe=0, write_xyz=False):

    """
    Unwrap coordinates using an external C OpenMP routine.

    Args:
        coords (ndarray): Wrapped coordinates, shape (natoms, 3).
        nmols_array (ndarray): Atom indices for each molecule; -1 marks padding.
        l_neigh_array (ndarray): Bonded neighbors for each atom.
        box_dimensions (ndarray): Box lengths.
        iframe (int): Frame number.
        write_xyz (bool): Write XYZ output if True.

    Returns:
        ndarray: Unwrapped coordinates.
    """

    from topology.ext_libc.c_unwrap_openmp import calc_unwrap_openmp

    natoms = coords.shape[0]
    coords_unwrapped = np.zeros((natoms, 3), dtype=np.float32)

    coords = np.ascontiguousarray(coords)

    calc_unwrap_openmp(nmols_array, l_neigh_array, coords, box_dimensions, coords_unwrapped)

    if write_xyz:
        filename = "movie_unwrap_{0:06d}.xyz".format(iframe)

        with open(filename, 'w') as f:
            f.writelines("{0:d}\n".format(natoms))
            f.writelines("Frame wrapped: {}\n".format(iframe))
            for iatom in range(natoms):
                line = "{0:s} {1:.4f} {2:.4f} {3:.4f}\n".\
                    format('C', coords[iatom, 0], coords[iatom, 1], coords[iatom, 2])
                f.writelines(line)

            f.writelines("{0:d}\n".format(natoms))
            f.writelines("Frame wrapped: {}\n".format(iframe))
            for iatom in range(natoms):
                line = "{0:s} {1:.4f} {2:.4f} {3:.4f}\n".\
                    format('C', coords_unwrapped[iatom, 0], coords_unwrapped[iatom, 1], coords_unwrapped[iatom, 2])
                f.writelines(line)

    return coords_unwrapped


# #############################################################################
def unwrap_nojump(coords, coords_0, n_image, nmols_array, l_neigh_array, box_dimensions, iframe=0, write_xyz=False):

    """
    Unwrap coordinates without inter-frame jumps due to PBCs, using a C
    OpenMP routine.

    Args:
        coords (ndarray): Wrapped coordinates for the current frame.
        coords_0 (ndarray): Previous-frame unwrapped coordinates.
        n_image (ndarray): Image vectors.
        nmols_array (ndarray): Molecule definitions.
        l_neigh_array (ndarray): Neighbor list.
        box_dimensions (ndarray): Box lengths.
        iframe (int): Frame index.
        write_xyz (bool): Write XYZ files if True.

    Returns:
        ndarray: Unwrapped coordinates.
    """


    from topology.ext_libc.c_unwrap_openmp import calc_unwrap_nojump_openmp

    natoms = coords.shape[0]
    coords_unwrapped_nojump = np.zeros((natoms, 3), dtype=np.float32)

    coords = np.ascontiguousarray(coords)
    coords_0 = np.ascontiguousarray(coords_0)
    calc_unwrap_nojump_openmp(nmols_array, l_neigh_array, coords, coords_0,
                              n_image, box_dimensions, coords_unwrapped_nojump)

    if write_xyz:
        filename = "movie_unwrap_nojump_{0:06d}.xyz".format(iframe)

        with open(filename, 'w') as f:
            f.writelines("{0:d}\n".format(natoms))
            f.writelines("Frame wrapped: {}\n".format(iframe))
            for iatom in range(natoms):
                line = "{0:s} {1:.4f} {2:.4f} {3:.4f}\n".\
                    format('C', coords[iatom, 0], coords[iatom, 1], coords[iatom, 2])
                f.writelines(line)

            f.writelines("{0:d}\n".format(natoms))
            f.writelines("Frame wrapped: {}\n".format(iframe))
            for iatom in range(natoms):
                line = "{0:s} {1:.4f} {2:.4f} {3:.4f}\n".format('C', coords_unwrapped_nojump[iatom, 0],
                                                                coords_unwrapped_nojump[iatom, 1],
                                                                coords_unwrapped_nojump[iatom, 2])
                f.writelines(line)

    return coords_unwrapped_nojump

