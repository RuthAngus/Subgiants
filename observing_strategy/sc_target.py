import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps
import scaling_relations as sr
from colours import plot_colours
ocols = plot_colours()
from sin_tests import fit_sine

def flux_rv(y, y_err, teff, t_err):
    dlL = (y - np.median(y)) * 1e6  # ppm
    dlL_err = (y_err/y) * 1e6
    return dlL/17.7*teff/5777., np.sqrt((dlL_err/dlL)**2+(t_err/teff)**2), dlL

def closest_point(times, x):
    xs = []
    for i in range(len(times)):
        dist_2 = (x - times[i])**2
        xs.append(np.argmin(dist_2))
    return np.array(xs)

def hd185_rv(nn, s):

    DIR = '/Users/angusr/Python/Subgiants'
    x, y = np.genfromtxt('%s/data/hd185351.q16sc.ts' % DIR).T
    y_err = np.ones_like(y)*6.255e-5

#     plt.clf()
#     plt.errorbar(x[:100], y[:100], yerr=y_err[:100])
#     plt.show()
#     raw_input('enter')

    # physical parameters
    teff = 5042
#     t_err = 32
    t_err = 0
    logg = 3.28
    rho = 0.013
    m = 1.99
    r = 5.35
    nm_alt = sr.nu_max_alt(logg, teff)/1e3
    dn_alt = sr.delta_nu_alt(rho)/1e6  # total guess
    nm = sr.nu_max(m, r, teff)/1e3
    dn = sr.delta_nu(m, r)/1e6

    # convert flux to rvs
    rv, rv_err, dlL = flux_rv(y, y_err, teff, t_err)

    # select an nn+1 day sample at random
    np.random.seed(1234)
    r = np.random.uniform(0, len(x)-int((nn+1)/(x[1]-x[0])))
    l = (x[r] < x) * (x < x[r]+11)
    x, rv, rv_err = x[l], rv[l], rv_err[l]

    # select one point per day at 1/2 a day in
    diff = x[1] - x[0]
    target_times = np.arange(.5+s*diff, nn+.5+s*diff, 1.)
    l = closest_point((x[0]+target_times), x)
    ts, rvs, rv_errs = x[l], rv[l], rv_err[l]

#     # plot
#     plt.clf()
#     plt.errorbar(ts, rvs, yerr=rv_errs, fmt='mo', ecolor='m')
#     plt.show()
    return ts, rvs, rv_errs

if __name__ == "__main__":

    hd185_rv(10, 3, 1)
