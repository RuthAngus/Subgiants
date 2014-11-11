import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
reb = plot_params()
from astero_modelling import MCMC, Gaussian_priors, lnlike, gen_freqs
from HD82074 import HD
HD = HD()

# SME estimates:
m =.93
m_err = .06
r = 3.95
r_err = 0.13
teff = 4996
t_err = 44
L = 8.64
l_err = 0.62

# asteroseismic/interferometric results
# m = 1.2
# m_err = .11
# r = 3.96
# r_err = 0.12

# multi order
t, rv, rv_err = HD.t-HD.t[0], HD.rv, HD.rv_err
rv -= np.median(rv)

nfreqs = 12
par_init = [m, r, teff, m_err, r_err, t_err]
args = t, rv, rv_err, nfreqs, lnlike, Gaussian_priors
burnin = 1000
runs = 20000
fname = 'HD872074'
fig_labels = ['$M$', '$R$', '$T_{eff}$']

mres, mcmc_result = MCMC(par_init, args, burnin, runs, fname, fig_labels)
print mres
