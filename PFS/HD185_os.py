import numpy as np
import matplotlib.pyplot as plt
import pyfits
import sc_target
from rc_params import plot_params
reb, fbt = plot_params()
from regularised_sine_fitting import fit_sine_reg
from scipy.signal import butter, lfilter, lombscargle
from sin_tests import show_sine

DIR = "/Users/angusr/Python/Subgiants"

# load kids, masses and temperatures
name = 185351
teff = 5016
m = 1.99
radius = 5.35
import scaling_relations as sr
dnu = sr.delta_nu(m, r)**1e-6
n = sr.nu_max(m, r, t)**1e-3

# calculate exposure times
PFS = (8., 88.)  # 88 secs on a 8 mag star for S/N = 188
Vg, expg = PFS
exptime = exp_time(vmag, Vg, expg)

nmins = 2  # interval between observations in minutes
ndays = 10  # number of nights observed. if this is greater than the no.
# of nights in the simulated data, the rms will increase
ntests = 50  # number of changes in interval
nsamples = [2, 3, 5]  # number of observations per night
nsim = 1  # number of simulations
exptime = 100

os(nmins, ndays, ntests, sample, nsim, exptime, name)
