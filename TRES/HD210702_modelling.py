import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
reb = plot_params()
from astero_modelling import MCMC, lnlike, Gaussian_priors, lnprior_alt
import scaling_relations as sr

t0 = 2456900

# multi order
BJD, rv, rv_err = np.genfromtxt('HD210702.ccfSum.txt', skip_header=1,
                                  usecols=(0, 1, 2)).T
rv -= np.median(rv)
BJD -= t0

# M, R, teff, logg, rho, teff
m = 1.64491
r = 4.47
teff = 5015.45
nfreqs = 12
logg = 3.36112
rho = m / r**3  # solar units

# with m, r, teff
# sigs = [0.0604717, 0.165514, 44.]
# par_init = [m, r, teff]
# args = BJD, rv, rv_err, nfreqs, lnlike, lnprior, sigs
# burnin = 1000
# runs = 20000
# fname = 'HD210702'
# fig_labels = ['$M$', '$R$', '$T_{eff}$']

# with logg, rho, teff
sigs = [0.0600000, 0.1, 44.]
par_init = [logg, rho, teff]
args = BJD, rv, rv_err, nfreqs, lnlike, Gaussian_priors, sigs
# args = BJD, rv, rv_err, nfreqs, lnlike, lnprior_alt
burnin = 1000
runs = 20000
fname = 'HD210702'
fig_labels = ['$\log~g$', '$\rho$', '$T_{eff}$']
mres = MCMC(par_init, args, burnin, runs, fname, fig_labels)
print mres
