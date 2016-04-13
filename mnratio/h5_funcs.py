import numpy as np
import h5py as h5
from astropy.io import fits


def mk_h5(infile, outfile, dset_name, mode='a', **kwargs):
    """
    This function takes an ASCII file with space-separated columns of data
    and makes a HDF5 file of the desired columns of data, skipping an
    optional number of rows and assigning a dataset name.
    """
    rows = kwargs.pop("skiprows", 0)
    cols = kwargs.pop("usecols", None)

    dat = np.loadtxt(infile, skiprows=rows, usecols=cols)

    f = h5.File(outfile, mode)

    dset = f.create_dataset(dset_name,
                            shape=dat.shape,
                            dtype=dat.dtype)

    dset[:] = dat[:]

    f.close()


def arr2h5(arr, outfile, dset_name, mode='a'):
    """
    This function takes a numpy array and saves it to a HDF5 file.
    """
    f = h5.File(outfile, mode)

    dset = f.create_dataset(dset_name,
                            shape=arr.shape,
                            dtype=arr.dtype)

    dset[:] = arr[:]

    f.close()


def h5_arr(h5file, dset):
    """
    This function copies the desired array from a hdf5 file and returns
    it.
    """
    f = h5.File(h5file, 'r')

    arr = f[dset][...]

    f.close()

    return arr


def fits2h5(fitsfile, HDU, col_list, h5file, dset):
    """
    This function takes vectors from a fitsfile HDU and puts them into
    columns of an array in a HDF5 dataset.

    Parameters
    ----------

    fitsfile: string
        address of a FITS file
    HDU: int
        the number of the header data unit to get the data from
    col_list: list
        a list of strings that correspond to the names of the columns
        to be extracted
    h5file: string
        name of HDF5 file to be created
    dset: string
        name of the dataset to be created in the HDF5 file to hold the array
        of FITS columns
    """

    dat = fits.open(fitsfile)[HDU].data

    arr = np.zeros((dat.shape[0], len(col_list)))

    for i, col in enumerate(col_list):
        arr[:, i] = dat[col]

    arr2h5(arr, h5file, dset)
