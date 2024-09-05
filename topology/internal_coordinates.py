"""
Miscellany of mathematical functions
"""
# try:
#     from numba import jit
#     import numba
# except ModuleNotFoundError:
#     pass
import math
import numpy as np
import random
import psutil
from scipy.spatial import distance_matrix


# #############################################################################
# noinspection PyUnboundLocalVariable
def distance_array(ref, conf, openmp=False):

    """
    Calculate the distances among the points in two arrays using C and Cython

    PBCs are not corrected in the function

    Parameters:
        * ``ref``: (type: ndarray): Reference array containing a number of points with shape [npoints1, dim],
        where dim is the dimension in 2-D or 3-D points
        * ``conf``: (type: ndarray): Reference array containing a number of points with shape [npoints2, dim],
        where dim is the dimension in 2-D or 3-D points
        * ``openmp``: (type: boolean) Parallel or not

    Returns:
        * ``dist``: (type: 2D-ndarray). Distances in Angstroms between the ith atom (row) and the jth atoms
        * ``rijx``: (type: 2D-ndarray). Component-x of the vector i,j
        * ``rijy``: (type: 2D-ndarray). Component-y of the vector i,j
        * ``rijz``: (type: 2D-ndarray). Component-z of the vector i,j

    ``Examples``:


    """

    if openmp:
        from ext_libc.c_distances_openmp import calc_distance_array_openmp
        # from c_distances_openmp import calc_distance_array_openmp
    else:
        from ext_libc.c_distances import calc_distance_array
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

    Calculate the distances among the points in two arrays using C and Cython

    ref  = [[atom1],[atom2]]
    conf = [[atom3],[atom4]]
    Distances between atom1-atom3 and atom2-atom4

    dist = [ distance_atom1_atom3, distance_atom2_atom4 ]

    PBCs are not corrected in the function

    ``Parameters``:
        * **ref**: (type: ndarray): Reference array containing a number of points with shape
        [npoints, dim], where dim is the dimension in 2-D or 3-D points
        * **conf**: (type: ndarray): Reference array containing a number of points with shape
        [npoints, dim], where dim is the dimension in 2-D or 3-D points
        * **openmp**: (type: boolean) Parallel or not

    ``Returns``:
        * **dist**: (type: 1D-ndarray). Distances in Angstroms between the ith atom and jth atoms in the same row of the ref and coor arrays
        * **rijx**: (type: 1D-ndarray). Component-x of the vector i,j
        * **rijy**: (type: 1D-ndarray). Component-y of the vector i,j
        * **rijz**: (type: 1D-ndarray). Component-z of the vector i,j

    ``Examples``:


    """

    if openmp:
        from ext_libc.c_distances_openmp import calc_distance_diagonal_openmp
        # from c_distances_openmp import calc_distance_array_openmp
    else:
        from ext_libc.c_distances import calc_distance_diagonal
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

    """Calculate the distances among the points in two arrays using pure pyhton and numpy

    PBCs are not corrected in the function

    ``Parameters``:
        * **ref**: (type: ndarray): Reference array containing a number of points with shape
         [npoints1, dim], where dim is the dimension in 2-D or 3-D points
        * **conf**: (type: ndarray): Reference array containing a number of points with shape
         [npoints2, dim], where dim is the dimension in 2-D or 3-D points
        * **openmp**: (type: boolean) Parallel or not

    ``Returns``:
        * **dist**: (type: 2D-ndarray). Distances in Angstroms between the ith atom (row) and the jth atoms
        * **rijx**: (type: 2D-ndarray). Component-x of the vector i,j
        * **rijy**: (type: 2D-ndarray). Component-y of the vector i,j
        * **rijz**: (type: 2D-ndarray). Component-z of the vector i,j
    ``Examples``:
    >>>
    >>>
    >>>
    >>>

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

    """Calculate the distances among the points in two arrays using pure pyhton

    PBCs are not corrected in the function

    ``Parameters``:
        * **ref**: (type: ndarray): Reference array containing a number of points with shape
        [npoints1, dim], where dim is the dimension in 2-D or 3-D points
        * **conf**: (type: ndarray): Reference array containing a number of points with shape
         [npoints2, dim], where dim is the dimension in 2-D or 3-D points
        * **openmp**: (type: boolean) Parallel or not

    ``Returns``:
        * **dist**: (type: 2D-ndarray). Distances in Angstroms between the ith atom (row) and the jth atoms

    ``Examples``:
    >>>
    >>>
    >>>
    >>>

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
    Bend angle

    Finds angle between three atomic positions.
    Periodic boundary conditions are not taken into account, thus, coordinates must be unwrapped.
    This implementation is made in pure python, so it can be slow

        (1) --- (2) --- (3) (Bend angle 1-2-3)

    Parameters:
        * ``c1`` : Coordinates of the point 1
        * ``c2`` : Coordinates of the point 2
        * ``c3`` : Coordinates of the point 3
        * ``radians``: returns value in radians (degrees) if True (False)

    Returns:
        angle between particles in radians or degrees

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

    Calculate the dihedral or improper angle
    PBCs are not corrected in the function

    .. image:: ../../figures/dih_imp.png

    Parameters:
        * ``ci`` :
        * ``cj`` :
        * ``ck`` :
        * ``cl`` :
        * ``radians`` :

    Returns:
        * ``phi`` :
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

    Calculate the dihedral or improper angle
    PBCs are not corrected in the function

    .. image:: ../../figures/dih_imp.png

    Parameters:
        * ``ci`` :
        * ``cj`` :
        * ``ck`` :
        * ``cl`` :
        * ``radians`` :

    Returns:
        * ``phi`` :
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
    Calculate the cos angle, used for bond perception in the paper
    Zhang et al. Journal of Cheminformatics 2012, 4:26 (Figure 4)
    (https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-4-26)

    PBCs are not corrected in the function

    .. image:: ../../figures/cos_angle.png

    Parameters:
        * ``ci`` :
        * ``cj`` :
        * ``ck`` :
        * ``cl`` :
        * ``radians`` :

    Returns:
        angle between particles in radians or degrees

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

    Generate Euler angles.

    To generate an uniform sampling of the three Euler angles:

    .. image:: euler_angles_random.png

    There are many definitions of the Euler angles
    (see: https://en.wikipedia.org/wiki/Euler_angles)

    ``Parameters``:
        * **iseed**  (type: int): Seed for the random number generator

    ``Return``:
        * **A vector** containg the three Euler angles in radians

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

    Create a rotation matrix for a given set of Euler angles
    There are many definitions of the Euler angles
    (see: https://en.wikipedia.org/wiki/Euler_angles)

    The definition here used is that given in:

    .. code-block::

        MATHEMATICAL METHODS FOR PHYSICISTS
        SEVENTH EDITION
        George B. Arfken, Hans J. Weber, Frank E. Harris
        pag: 140-142

    .. image:: euler_book.png

    ``Parameters``:
        * **None**

    ``Returns``:
        * **None**

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
    Calculate the center of geometry of a set of coordinates.
    Periodic boundary conditions are not taken into account.
    This implementation is made in pure python, so it can be slow

    Parameters:
        * ``coords``: (type: list of ndarray-float32 (3)): Coordinates of the atoms. It must be unwrapped.

    Return:
        * ``cog`` (type: ndarray vector): Coordinates of the geometry center.

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
    Calculate the center of geometry of a set of coordinates.
    Periodic boundary conditions are not taken into account.
    This implementation is made in cython

    Parameters:
        * ``coords``: (type: list of ndarray-float32 (3)): Coordinates of the atoms. It must be unwrapped.

    Return:
        * ``cog`` (type: ndarray vector): Coordinates of the geometry center.

    """

    from ext_libc.c_unwrap_openmp import calc_center_of_geom

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
    Calculate the center of mass of a set of coordinates.
    Periodic boundary conditions are not taken into account.
    This implementation is made in pure python, so it can be slow

    Parameters:
        * ``coords``: (type: list of ndarray-float32 (3)): Coordinates of the atoms. It must be unwrapped.
        * ``mass`` : (type: ndarray-float32 (natoms) : Masses of the atoms

    Return:
        * ``com`` (type: ndarray vector): Coordinates of the geometry center.

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
    Calculate the center of geometry of a set of coordinates.
    Periodic boundary conditions are not taken into account.
    This implementation is made in cython

    Parameters:
        * ``coords``: (type: list of ndarray-float32 (3)): Coordinates of the atoms. It must be unwrapped.
        * ``mass`` : (type: ndarray-float32 (natoms) : Masses of the atoms

    Return:
        * ``cog`` (type: ndarray vector): Coordinates of the geometry center.

    """

    from ext_libc.c_unwrap_openmp import calc_center_of_mass

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

    Unwrapping a set of coordinates. Use pure python

    Args:
        * ``coords`` (np.array, [natoms, 3]): Wrapped coordinates
        * ``topology`` (topology object): Topology object
        * ``box_dimensions`` (np.array, [3]): Cubic box length
        * ``iframe`` (integer): Number of the trajectory frame
        * ``write_xyz`` (boolean): Yes to write xyz file for each frame

    Returns:
        * ``coords_unwrap`` (np.array, [natoms, 3]): Unwrapped coordinates
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

    R"""Unwrapping a set of coordinates. Use a external c library.

    Parameters
    ----------
    coords : np.array, [natoms, 3]
        Wrapped coordinates
    nmols_array : np.array [nmols, maxatch]
        Global numbering of the atoms in each molecule (row). The value -1 indicates a position to be not considered.
        This is the case of molecules with different number of atoms.
    l_neigh_array : np.array[natoms,4]
        l_neigh_array[i] = [j, k, -1, -1] j and k are bonded to i.
    box_dimensions :
        Cubic box length
    iframe : integer
        Number of the current trajectory frame
    write_xyz : boolean
        Yes to write xyz file for each frame


    Returns
    -------
    coords_unwrapped :
        Unwrapped coordinates
    """

    from ext_libc.c_unwrap_openmp import calc_unwrap_openmp

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

    R"""Unwrapping without jumps between frames  due to PBC a set of coordinates.
    Use a external c library.

    Parameters
    ----------
    coords : np.array, [natoms, 3]
        Wrapped coordinates
    coords_0 : np.array, [natoms, 3]
        UnWrapped coordinates of the previous frame in the trajectory
    n_image : np.array, [natoms, 3]
        Image
    nmols_array : np.array [nmols, maxatch]
        Global numbering of the atoms in each molecule (row). The value -1 indicates a position to be not considered.
        This is the case of molecules with different number of atoms.
    l_neigh_array : np.array[natoms,4]
        l_neigh_array[i] = [j, k, -1, -1] j and k are bonded to i.
    box_dimensions :
        Cubic box length
    iframe : integer
        Number of the current trajectory frame
    write_xyz : boolean
        Yes to write xyz file for each frame

    Returns
    -------
    coords_unwrapped :
        Unwrapped coordinates
    """

    from ext_libc.c_unwrap_openmp import calc_unwrap_nojump_openmp

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

