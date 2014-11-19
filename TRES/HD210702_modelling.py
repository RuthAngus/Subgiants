import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
reb = plot_params()
from astero_modelling import MCMC, lnprior, lnlike, gen_freqs, Gaussian_priors

t0 = 2456900

# multi order
BJD, rv, rv_err = np.genfromtxt('HD210702.ccfSum.txt', skip_header=1,
                                  usecols=(0, 1, 2)).T
rv -= np.median(rv)
BJD -= t0

nfreqs = 12
# par_init = [1., 2., 5500.]  # M, R, teff (complete guess)

# M, R, teff (complete guess)
par_init = [1.64491, 4.47, 5015.45, 0.0604717, 0.165514, 44.]
args = BJD, rv, rv_err, nfreqs, lnlike, Gaussian_priors
burnin = 1000
runs = 20000
fname = 'HD210702'
fig_labels = ['$M$', '$R$', '$T_{eff}$']
# fig_labels = ['$M$', '$T_{eff}$']

mres = MCMC(par_init, args, burnin, runs, fname, fig_labels)
print mres
