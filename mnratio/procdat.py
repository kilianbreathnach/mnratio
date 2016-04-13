import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)
import h5py as h5
from scipy.optimize import minimize
from astropy.io import fits
# from astropy import units as u
# from astropy.cosmology import Planck13, WMAP5
# from astropy.coordinates import Angle, Distance, ICRS
# from cartesian_cosmo import mk_coords
from cosmo_dis import rdz2cart

# Planck13.__init__(100.0, Planck13.Om0)
# cosmo = Planck13


def fits2rdz(fitsfile, outfile):
    """
    This function opens a BOSS fits file and extracts the RA, Dec, and
    redshift of the galaxies and saves it to a hdf5 file.
    """
    dat = fits.open(fitsfile)[1].data
    ra = dat["RA"]
    dec = dat["DEC"]
    z = dat["Z"]

    arr = np.zeros((ra.shape[0], 3))

    arr[:, 0] = ra
    arr[:, 1] = dec
    arr[:, 2] = z

    f_out = h5.File(outfile)
    rdz = f_out.create_dataset("radecz", shape=arr.shape, dtype='float64')
    rdz[:, :] = arr

    f_out.close()


def minfunc(fac, *args):
    """
    function to minimize the amount to downsample the redshift bins by
    """
    gals = args[0]
    clust = args[1]
    return np.abs(np.min(gals - fac * clust))


def down_to_red():
    """
    Take the HDF5 RDZ file and downsample to have a number density profile in
    redshift that is proportional to that of the clusters in redmapper.
    """
    redat = fits.open("../redmapper/redmapper_dr8_public_v5.10_catalog.fits")[1].data
    redzs = redat["Z_LAMBDA"]
    posids = np.where((redzs > 0) * (redzs < 0.45))[0]
    redat = redat[posids]
    redzs = redat["Z_LAMBDA"]

    ngc_ras = np.where((redat["RA"] > 80) * (redat["RA"] < 295))
    sgc_ras = np.where((redat["RA"] < 80) + (redat["RA"] > 295))

    ngc_zs = redzs[ngc_ras]
    sgc_zs = redzs[sgc_ras]

    fitsl = ["../bossdat/lowz-dr12v4-N-Reid.dat.fits",
             "../bossdat/lowz-dr12v4-S-Reid.dat.fits"]
    outfns = ["../data/LOWZ/NGC/rdz.hdf5",
              "../data/LOWZ/SGC/rdz.hdf5"]
    figs = ["NGC_nowbars.pdf", "SGC_nowbars.pdf"]

    for i, cap in enumerate([ngc_zs, sgc_zs]):

        nbar = np.histogram(cap, bins=25)
        nbins = nbar[1]
        nbar_p = nbar[0] / np.float(np.sum(nbar[0]))

        # gal_fits = fits.open("./bossdat/lowz-dr12v4-N-Reid.dat.fits")
        gal_fits = fits.open(fitsl[i])
        galdat = gal_fits[1].data
        zarr = galdat["Z"]

        galn  = np.histogram(zarr, bins=nbins)
        galhist = galn[0]
        gal_p = galhist / np.float(np.sum(galhist))

#         res = minimize(minfunc, 0.5, args=(gal_p, nbar_p))
#
#         # factor to downsample each galaxy redshift bin to
#         downfac = (res.x[0] * nbar_p) / gal_p

        q = 10.

        for k in range(len(galhist)):

            if q * nbar_p[k] > gal_p[k]:

                q *= gal_p[k] / (q * nbar_p[k])

        downfac = (q * nbar_p) / gal_p

        keepgals = np.rint(downfac * galhist)

        mask = np.array([False] * zarr.shape[0])

        for j in range(len(galhist)):

            print j

            bin_zs = np.where((zarr >= nbins[j]) * (zarr < nbins[j + 1]))[0]

            keep_ids = np.random.choice(bin_zs, size=keepgals[j], replace=False)

            mask[keep_ids] = True

        rdz = np.zeros((np.sum(mask), 3))

        rdz[:, 0] = galdat["RA"][mask]
        rdz[:, 1] = galdat["DEC"][mask]
        rdz[:, 2] = galdat["Z"][mask]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        nc, = ax.plot(nbar[1][:-1] + 0.5 * np.diff(nbar[1]), nbar_p * rdz.shape[0],
                      color='r', ls=':', alpha = 0.6,
                      label="redmapper proportional number density")

        finhist = np.histogram(rdz[:, 2], bins=nbins)
        ng, = ax.plot(finhist[1][:-1] + 0.5 * np.diff(finhist[1]), finhist[0],
                      color='b', alpha = 0.6,
                      label="downsampled galaxy number density")

        ofac = float(np.sum(finhist[0])) / np.sum(galn[0])

        no, = ax.plot(galn[1][:-1] + 0.5 * np.diff(galn[1]), ofac * galn[0],
                      color='g', ls='--', alpha = 0.6,
                      label="original lowz galaxy number density")

        ax.set_xlabel(r"$z$")
        ax.set_ylabel(r"$N$")
        ax.set_title("compare redmapper and lowz number densities")

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)

        picname = figs[i]
        plotaddress = "../plotting/" + picname

        fig.savefig(plotaddress)

        os.system("scp {0} broiler:~/public_html/tinker/".format(plotaddress))
        os.system("ssh broiler chmod 644 /home/kilian/public_html/tinker/{0}"
                      .format(picname))

        f_out = h5.File(outfns[i])
        dset = f_out.create_dataset("radecz", shape=rdz.shape, dtype='float64')
        dset[:, :] = rdz
        f_out.close()
        print "done with ", i


def rdz2xyz(infile, outfile):

    f_in = h5.File(infile)

    radecz = f_in["radecz"]

    f_out = h5.File(outfile)

    cart = f_out.create_dataset("cart_pts",
                                shape=(radecz.shape[0], 3), dtype='float64')

#    for i in range(radecz.shape[0]):
#        ra = Angle(radecz[i, 0], u.deg)
#        dec = Angle(radecz[i, 1], u.deg)
#        losd = cosmo.comoving_distance(radecz[i, 2])
#        dis = Distance(losd)
#        coord = ICRS(ra=ra, dec=dec, distance=dis)
#        cart[i, :] = np.array([coord.cartesian.x.value,
#                               coord.cartesian.y.value,
#                               coord.cartesian.z.value])


    rdzarr = np.zeros((radecz.shape[0], 3))
    rdzarr[:, :] = radecz[:, :]
    xyzarr = rdz2cart(rdzarr)
    cart[:, :] = xyzarr

    f_in.close()
    f_out.close()

if __name__ == "__main__":

    import sys

    if len(sys.argv) != 5:
        print "usage: python procdat.py <data fits file> <radecz file (hdf5)>\
                 <cartesian coords file (hdf5)> <cosmology: Planck | WMAP>"

    fits2rdz(sys.argv[1], sys.argv[2])

    mk_coords(sys.argv[2], sys.argv[3], sys.argv[4])

