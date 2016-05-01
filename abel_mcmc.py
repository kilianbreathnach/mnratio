import numpy as np
from scipy.integrate import simps
import emcee


def lnprior(theta):
    pass


def lonpost(theta):

    gen_proj = np.dot(matrix, theta)
    pass


def mk_matrix(r_edges, y_edges):

    # get matrix size and shape from inputs
    mat = np.zeros((len(y_edges) - 1, len(r_edges) - 1))

    for i_y in range(len(y_edges) - 1):

        # define projected radius limits for this bin
        y1 = y_edges[i_y]
        y2max = y_edges[i_y + 1]

        # compute the factor for the bin area
        prefac = 4. / (3 * (y2max ** 2 - y1 ** 2))

        k = 0
        # loop over 3D radial density bins that lie in the projected bin
        for i_r in np.where(r_edges >= y1)[0]:

            r1 = r_edges[i_r - 1]
            r2 = r_edges[i_r]

            # ylimit for integral is r2 if it does not extend beyond ymax
            if r2 < y2max:
                y2 = r2
            else:
                y2 = y2max

            # for the first density bin, we integrate from r1 = y1
            if k == 0:
                postfac = (r2 ** 2 - y1 ** 2) ** 1.5 - \
                          (r2 ** 2 - y2 ** 2) ** 1.5
                mat[i_y, i_r - 1] = prefac * postfac
                k = 1
            else:
                postfac = (r2 ** 2 - y1 ** 2) ** 1.5 + \
                          (r1 ** 2 - y2 ** 2) ** 1.5 - \
                          (r2 ** 2 - y2 ** 2) ** 1.5 - \
                          (r1 ** 2 - y1 ** 2) ** 1.5
                mat[i_y, i_r - 1] = prefac * postfac

                print postfac

    return mat


y_edges = np.array([80])
for i in range(36):
    y_edges = np.append(10 ** (np.log10(y_edges[0]) - 1./12), y_edges)

r_edges = np.append(0., np.logspace(-.5, 1, 5))

print mk_matrix(r_edges, y_edges)





# First get the cross-correlation array

xcorr_dat = np.loadtxt("xi_sp_fullcross.dat", usecols=(0, 1, 2))

rp = np.empty((36,))
wpx = np.empty((36,))

for i in range(36):

    rp[i] = xcorr_dat[80 * i, 1]
    wpx[i] = simps(xcorr_dat[i * 80:80 * (i + 1), 2],
                  xcorr_dat[i * 80:80 * (i + 1), 0])



