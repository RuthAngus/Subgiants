import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
reb = plot_params()
from astero_modelling import MCMC, lnlike, Gaussian_priors, lnprior_alt
from astero_modelling import model, gen_freqs_alt
import scaling_relations as sr
import scipy.signal as sps

def plot_sine(par_init, x, A, nfreqs):
    logg, rho, teff = par_init
    w = gen_freqs_alt(logg, rho, teff, nfreqs)
    print w, A
    ys = np.zeros_like(x)
    for i in range(len(w)):
        ys += A[2*i]*np.sin(w[i]*x) + A[2*i+1]*np.cos(w[i]*x)
    ys += A[-1]
    return ys

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
nfreqs = 1
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

# plt.clf()
# plt.errorbar(BJD, rv, yerr=rv_err, **reb[0])
# ys, A = model(par_init, BJD, rv, rv_err, nfreqs)
# plt.plot(BJD, ys)
# xs = np.linspace(min(BJD), max(BJD), 1000)
# ys = plot_sine(par_init, BJD, A, nfreqs)
# # plt.plot(BJD, ys)
# plt.show()
# raw_input('enter')

# periodogram
fs = np.linspace(3, 50, 10000) # c/d
ws = 2*np.pi*fs  # lombscargle uses angular frequencies
pgram = sps.lombscargle(BJD, rv, ws)
plt.clf()
plt.plot((1./fs)*24, pgram)
plt.show()
raw_input('enter')

# args = BJD, rv, rv_err, nfreqs, lnlike, Gaussian_priors, sigs
args = BJD, rv, rv_err, nfreqs, lnlike, lnprior_alt
burnin = 1000
runs = 20000
fname = 'HD210702'
fig_labels = ['$\log~g$', '$\rho$', '$T_{eff}$']
mres = MCMC(par_init, args, burnin, runs, fname, fig_labels)
print mres
