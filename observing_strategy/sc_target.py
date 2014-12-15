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

def hd185_rv(x, rv, rv_err, nn, s):

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

    return ts, rvs, rv_errs
