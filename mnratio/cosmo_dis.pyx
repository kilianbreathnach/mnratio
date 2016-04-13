import numpy as np
cimport numpy as np
cimport cython

from scipy.integrate import quad
from astropy.cosmology import scalar_inv_efuncs

cdef extern from "math.h":
    long double sin(double x)
    long double cos(double x)

ctypedef np.float64_t DTYPE_t

cdef double pi = 3.14159265359

# for Planck13, with h = 1
cdef double hub_dist = 2997.92458  # in Mpc
cdef double Om0 = 0.30712
cdef double Ode0 = 0.6928382302273349
cdef double Ok0 = 0.
cdef double Ogamma0 = 2.470990199638946e-05
cdef double NeffPerNu = 1.013333333333333
cdef int nmasslessnu = 3
# cdef np.ndarray[DTYPE_t, ndim=1] nu_y = np.array([357.9121574])
nu_y = [357.9121574]


def inv_efunc(double z):

    return scalar_inv_efuncs.lcdm_inv_efunc(z, Om0, Ode0, Ok0, Ogamma0,
                                            NeffPerNu, nmasslessnu, nu_y)


def deg2rad(double deg):
    cdef double rad
    rad = deg * pi / 180.
    return rad


def cosmo_dist(DTYPE_t z):
    cdef double dis
    dis = hub_dist * quad(inv_efunc, 0, z)[0]
    return dis


@cython.boundscheck(False)
@cython.wraparound(False)
def rdz2cart(np.ndarray[DTYPE_t,ndim=2] rdz):

    cdef np.ndarray[DTYPE_t,ndim=2] xyz = np.empty((rdz.shape[0],3))
    cdef double d
    cdef int i

    for i in xrange(xyz.shape[0]):

        d = cosmo_dist(rdz[i, 2])

        xyz[i, 0] = d * cos(deg2rad(rdz[i, 1])) * cos(deg2rad(rdz[i, 0]))
        xyz[i, 1] = d * cos(deg2rad(rdz[i, 1])) * sin(deg2rad(rdz[i, 0]))
        xyz[i, 2] = d * sin(deg2rad(rdz[i, 1]))

    return xyz
