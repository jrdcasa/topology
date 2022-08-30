
# import both numpy and the Cython declarations for numpy
import cython
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from copy import copy


# declare the interface to the C code
cdef extern from "calc_unwrap.c":
    #cdef bint USED_OPENMP
    void c_unwrap (int* mols, int* neigh, float* coords, float* box,
                   int nchains, int maxatomsch, int natoms, int maxneigh,
                   float* coords_unwrap)
    void c_center_of_geom(int natoms, float* coord_current_ich, float cg[3])

cdef extern from "calc_unwrap.c":
    void c_unwrap_nojump(int* mols, int* neigh, float* coords, float* coords_0,
                         int* n_image, float* box, int nchains, int maxatomsch, int natoms,
                         int maxneigh, float* coords_unwrap)
    void c_setup(int natoms)

cdef extern from "calc_unwrap.c":
    void c_center_of_geom(int natoms, float* coords, float cog[3])

cdef extern from "calc_unwrap.c":
    void c_center_of_mass(int natoms, float* coords, float* mass, float cog[3])

#OPENMP_ENABLED = True if USED_OPENMP else False

# ========================================================================================
def calc_unwrap_openmp(np.ndarray[int, ndim=2, mode="c"] mols,
                       np.ndarray[int, ndim=2, mode="c"] neigh,
                       np.ndarray[float, ndim=2, mode="c"] coords,
                       np.ndarray[float, ndim=1, mode="c"] box,
                       np.ndarray[float, ndim=2, mode="c"] coords_unwrap):

    cdef int nchains = mols.shape[0]
    cdef int maxatomsch = mols.shape[1]
    cdef int natoms = neigh.shape[0]
    cdef int maxneigh = neigh.shape[1]

    c_unwrap(&mols[0,0], &neigh[0,0], &coords[0,0], &box[0],
             nchains, maxatomsch, natoms, maxneigh, &coords_unwrap[0,0])

    return None

# ========================================================================================
def  calc_unwrap_nojump_openmp(np.ndarray[int, ndim=2, mode="c"] mols,
                              np.ndarray[int, ndim=2, mode="c"] neigh,
                              np.ndarray[float, ndim=2, mode="c"] coords,
                              np.ndarray[float, ndim=2, mode="c"] coords_unwrap_0,
                              np.ndarray[int, ndim=2, mode="c"] n_image,
                              np.ndarray[float, ndim=1, mode="c"] box,
                              np.ndarray[float, ndim=2, mode="c"] coords_nojump):



    cdef int nchains = mols.shape[0]
    cdef int maxatomsch = mols.shape[1]
    cdef int natoms = neigh.shape[0]
    cdef int maxneigh = neigh.shape[1]

    c_unwrap_nojump(&mols[0,0], &neigh[0,0], &coords[0,0], &coords_unwrap_0[0,0], &n_image[0,0], &box[0],
                    nchains, maxatomsch, natoms, maxneigh, &coords_nojump[0,0])

    return None

# ========================================================================================
def calc_center_of_geom(np.ndarray[float, ndim=2, mode="c"] coords,
                        np.ndarray[float, ndim=1, mode="c"] cog):

    cdef int natoms = coords.shape[0]

    c_center_of_geom(natoms, &coords[0,0], &cog[0])

    return None

# ========================================================================================
def calc_center_of_mass(np.ndarray[float, ndim=2, mode="c"] coords,
                        np.ndarray[float, ndim=1, mode="c"] mass,
                        np.ndarray[float, ndim=1, mode="c"] cog):

    cdef int natoms = coords.shape[0]
    cdef int dim2 = mass.shape[0]
    if natoms != dim2:
        print ("ERROR. The number of atoms in the coords and mass arrays must be equal!!!!")
        return None

    c_center_of_mass(natoms, &coords[0,0], &mass[0], &cog[0])

    return None