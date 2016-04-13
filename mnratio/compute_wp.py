import h5py
import numpy as np
from halotools.empirical_models import model_defaults
from halotools.mock_observables import clustering


rbins = np.logspace(-1, np.log10(50), num=25)
rbin_centers = 10 ** ((np.log10(rbins[1:]) + np.log10(rbins[0:-1]))/2.)

pibins = np.linspace(0, 40, num=40)
pibin_centers = (rbins[1:] + rbins[0:-1])/2.


def run_auto(datfile, randfile, outfile):

    fd = h5py.File(datfile)
    fr = h5py.File(randfile)

    cd = fd["cart_pts"]
    cr = fr["cart_pts"]

    dat = np.empty(cd.shape)
    rands = np.empty(cr.shape)

    dat[:, :] = cd[:, :]
    rands[:, :] = cr[:, :]

    wps = clustering.wp(dat, rbins, pibins, randoms=rands,
                        do_cross=False, estimator='Landy-Szalay',
                        N_threads='max')

    outarr = np.empty((rbin_centers.shape[0], 2))
    outarr[:, 0] = rbin_centers
    outarr[:, 1] = wps

    np.savetxt(outfile, outarr)
