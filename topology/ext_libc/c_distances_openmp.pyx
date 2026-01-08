# import both numpy and the Cython declarations for numpy
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

# Required for compatibility with NumPy >= 2.0
cdef extern from "numpy/arrayobject.h":
    pass

import numpy as np
cimport numpy as np

# Force use of the new API (no deprecated API allowed)
np.import_array()

# declare the interface to the C code
cdef extern from "calc_distance.c":
    #cdef bint USED_OPENMP
    void c_distance_array (double* ref, double* conf, double* dist,
                           double* rijx, double* rijy, double* rijz,
                           int n, int m, int dim)
    void c_distance_diagonal_array (double* ref, double* conf, double* dist,
                                    double* rijx, double* rijy, double* rijz,
                                    int n, int dim)

#OPENMP_ENABLED = True if USED_OPENMP else False

def calc_distance_array_openmp(np.ndarray[np.float64_t, ndim=2, mode="c"] ref,
                               np.ndarray[np.float64_t, ndim=2, mode="c"] conf,
                               np.ndarray[np.float64_t, ndim=2, mode="c"] dist,
                               np.ndarray[np.float64_t, ndim=2, mode="c"] rijx,
                               np.ndarray[np.float64_t, ndim=2, mode="c"] rijy,
                               np.ndarray[np.float64_t, ndim=2, mode="c"] rijz):

    cdef int rows = ref.shape[0]
    cdef int cols = conf.shape[0]
    dim1 = ref.shape[1]
    dim2 = conf.shape[1]
    if dim1 != dim2:
        print ("ERROR. The dimension of points to calculate distances must be equal!!!!")
        return dist

    c_distance_array(&ref[0,0], &conf[0,0], &dist[0,0], &rijx[0,0],
                     &rijy[0,0], &rijz[0,0], rows, cols, dim1)
    return None

def calc_distance_diagonal_openmp(np.ndarray[np.float64_t, ndim=2, mode="c"] ref,
                                  np.ndarray[np.float64_t, ndim=2, mode="c"] conf,
                                  np.ndarray[np.float64_t, ndim=1, mode="c"] dist,
                                  np.ndarray[np.float64_t, ndim=1, mode="c"] rijx,
                                  np.ndarray[np.float64_t, ndim=1, mode="c"] rijy,
                                  np.ndarray[np.float64_t, ndim=1, mode="c"] rijz):

    cdef int rows = ref.shape[0]
    cdef int cols = conf.shape[0]
    dim1 = ref.shape[1]
    dim2 = conf.shape[1]
    if dim1 != dim2:
        print ("ERROR. The dimension of points to calculate distances must be equal!!!!")
        return dist

    c_distance_diagonal_array(&ref[0,0], &conf[0,0], &dist[0], &rijx[0],
                     &rijy[0], &rijz[0], rows, dim1)

    return None
