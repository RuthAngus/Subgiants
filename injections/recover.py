import numpy as np
import matplotlib.pyplot as plt
from model import model
import emcee
import triangle
import h5py
from scipy.optimize import minimize
from scipy.signal import lombscargle

# This script recovers the planet parameters

# The likelihood function:
def lnlike(theta, t, rv_obs, rverr, M1, ecc):
    model_rv = model(theta, t, rverr, M1, ecc, wn=False)
    chisq = np.sum(((rv_obs - model_rv) / rverr)**2)
    return -chisq / 2.

# The -likelihood function:
def neglnlike(theta, t, rv_obs, rverr, M1, ecc):
    model_rv = model(theta, t, rverr, M1, ecc, wn=False)
    chisq = np.sum(((rv_obs - model_rv) / rverr)**2)
    return chisq / 2.

# # emcee prior and posterior:
# def lnprior(theta):  # FIXME
#     return 0.

def lnprior(m):
    # log(P), m1, T0, V0, omega
    if 0 < m[0] < 100 and 0 < m[1] < 100 and -10 < m[2] < 1000 and \
            -1000 < m[3] < 1000 and -360 < m[4] < 360:
        return 0.0
    return -np.inf

def lnprob(theta, t, rv_obs, rverr, M1, ecc):
    return lnprior(theta) + lnlike(theta, t, rv_obs, rverr, M1, ecc)

def MCMC(theta, x, y, yerr, M1, ecc, fname, n, sub):

    nwalkers, ndim = 32, len(theta)
    p0 = [theta+1e-4*np.random.rand(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
                                    args=(x, y, yerr, M1, ecc))
    p0, lp, state = sampler.run_mcmc(p0, 500)
    sampler.reset()
    p0, lp, state = sampler.run_mcmc(p0, 5000)

#     fig_labels = ['P', 'M2', 'T0', 'V0', 'omega']
#     flatchain = sampler.chain[:, 50:, :].reshape((-1, ndim))
#     fig = triangle.corner(flatchain, truths=theta, labels=fig_labels)
#     plt.savefig("%s_triangle" % fname)

    print "saving samples"
    f = h5py.File("%s/results/%s_%s_%s_samples" % (DIR, n, fname, sub), "w")
    data = f.create_dataset("samples", np.shape(sampler.chain))
    data[:, :] = np.array(sampler.chain)
    f.close()

    flat = sampler.chain[:, 50:, :].reshape((-1, ndim))
    mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                      zip(*np.percentile(flat, [16, 50, 84], axis=0)))
    np.savetxt("%s/results/%s_%s_%s_results.txt" % (DIR, n, fname, sub),
               mcmc_result)

    return mcmc_result

def periodogram(x, y, fname, n, sub):
    fs = np.linspace(2, 25, 1000)
    pgram = lombscargle(x, y, fs*2*np.pi)
    theta_true = np.genfromtxt("%s/params/%s_%s_params.txt" %
                               (DIR, n, fname)).T
    truep = theta_true[0]
    print truep
    plt.clf()
    plt.plot(fs, pgram)
    plt.title("%s" % str(truep))
    plt.xlabel("%s" % str(truep))
    plt.savefig("%s/results/%s_%s_%s_pgram" % (DIR, n, fname, sub))
    l = pgram==max(pgram)
    period = fs[l]
    np.savetxt("%s/results/%s_%s_%s_pgram.txt" % (DIR, n, fname, sub), period)
    return period

def opt(theta, t, rv_obs, rverr, M1, ecc, fname, n):
    print theta
    results = minimize(neglnlike, theta, args=(t, rv_obs, rverr, M1, ecc))
    print results
    raw_input('enter')

if __name__ == "__main__":

    fname = "HD185"

    # load data and params
    DIR = "/Users/angusr/Python/Subgiants/injections"

    N = 100
    sub = 20
    for n in range(N):
        print n
#         x, y, yerr = np.genfromtxt("%s/rv_curves/%s_%s_rvs.txt"
#                                    % (DIR, n, fname)).T
        x, y, yerr = np.genfromtxt("%s/rv_curves/%s_%s_rvs_%s.txt"
                                   % (DIR, n, fname, sub)).T

        # P, M2, T0, V0, omega
        theta_true = np.genfromtxt("%s/params/%s_%s_params.txt" %
                                   (DIR, n, fname))
        theta_init = [10, 1, 0, 0, 0]
        M1, M1_err = np.genfromtxt("%s/params/%s_mass.txt" % (DIR, fname))
        ecc = 0.

        # sample in log period
        theta_init[0] = np.log(theta_init[0])
#         results = MCMC(theta_init, x, y, yerr, M1, ecc, fname, n, sub)
        periodogram(x, y, fname, n, sub)
