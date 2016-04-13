import numpy as np
import h5py as h5
from astropy import units as u
from astropy.cosmology import Planck13, WMAP5
from astropy.coordinates import Angle, Distance, ICRS


"""
This module contains the function mk_coords which converts the RA, Dec, and
redshift values of galaxy redshift survey data to cartesian
coordinates determined by a particular cosmology. This function is called
by the run_expt script or can be manually called with arguments for which
files to convert and which cosmology to use. The two options are WMAP5
for the original multidark simulation and Planck13 for the newer run.
"""


def mk_coords(radecfile, outfile, cosmology):

    # Set the cosmology with h free
    if cosmology == "Planck":
        Planck13.__init__(100.0, Planck13.Om0)
        cosmo = Planck13
    elif cosmology == "WMAP":
        WMAP5.__init__(100.0, WMAP5.Om0)
        cosmo = WMAP5

    f_in = h5.File(radecfile)
    radecz = f_in["radecz"]

    f_out = h5.File(outfile)
    cart = f_out.create_dataset("cart_pts", shape=(radecz.shape[0], 3),
                                dtype='float64')

    for i in range(radecz.shape[0]):
        ra = Angle(radecz[i, 0], u.deg)
        dec = Angle(radecz[i, 1], u.deg)

        losd = cosmo.comoving_distance(radecz[i, 2])
        dis = Distance(losd)

        coord = ICRS(ra=ra, dec=dec, distance=dis)

        cart[i, :] = np.array([coord.cartesian.x.value,
                               coord.cartesian.y.value,
                               coord.cartesian.z.value])

    f_in.close()
    f_out.close()


if __name__ == "__main__":

    import sys

    if len(sys.argv) != 4:
        print "usage: python cartesian_cosmo.py <radecz file (hdf5)>\
                <output file (hdf5)> <cosmology: Planck | WMAP>"

    mk_coords(sys.argv[1], sys.argv[2], sys.argv[3])
