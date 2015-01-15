import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps
import scaling_relations as sr
from colours import plot_colours
ocols = plot_colours()
from sin_tests import fit_sine
from synthesise import load_data

def flux_rv(y, y_err, teff, t_err):

    # convert y to relative parts per million
    dlL = ((y / np.median(y)) - 1) * 1e6  # ppm
    dlL_err = (y_err / np.median(y)) * 1e6

    rv = dlL/17.7*teff/5777.
    rv_err = np.sqrt((dlL_err/dlL/17.7)**2+(t_err/teff/5777)**2)
    return rv, rv_err, dlL

def closest_point(times, x):
    xs = []
    for i in range(len(times)):
        dist_2 = (x - times[i])**2
        xs.append(np.argmin(dist_2))
    return np.array(xs)

def hd185_rv(x, rv, rv_err, nn, s, start):

    # select an nn+1 day sample at random
    r = np.random.uniform(0, len(x)-int((nn+1)/(x[1]-x[0])))
    l = (x[r] < x) * (x < x[r]+11)
    x, rv, rv_err = x[l], rv[l], rv_err[l]

    # select one point per day at 1/2 a day in
    diff = x[1] - x[0]
    target_times = np.arange(start+s*diff, nn+start+s*diff, 1.)
    l = closest_point((x[0]+target_times), x)
    ts, rvs, rv_errs = x[l], rv[l], rv_err[l]

    return ts, rvs, rv_errs

if __name__ == "__main__":

    DIR = "/Users/angusr/Python/Subgiants"

    # load kids, masses and temperatures
    kid, teff, t_err, m, m_err = \
            np.genfromtxt("%s/data/AMP_subgiants.txt" % DIR, skip_header=1).T

    # load fits file
    fname = str(int(kid[0]))
    x, y, yerr = load_data(fname)

    # convert to rv
    rv, rv_err, dlL = flux_rv(y, y_err, teff[0], t_err[0])

    plt.clf()
    plt.subplot(2, 1, 1)
    plt.errorbar(x, y, yerr=yerr, alpha=.3, **reb)
    plt.subplot(2, 1, 2)
    plt.errorbar(x, rv, yerr=rv_err, alpha=.3, **reb)
    plt.savefig("test")
