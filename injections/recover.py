import numpy as np
import matplotlib.pyplot as plt
from model import model
import emcee
import triangle
import h5py

# This script recovers the planet parameters

# The likelihood function:
def lnlike(theta, t, rv_obs, rverr, M1, ecc):
    model_rv = model(theta, t, rverr, M1, ecc, wn=False)
    chisq = np.sum(((rv_obs - model_rv) / rverr)**2)
    return -chisq / 2.

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

def MCMC(theta, x, y, yerr, M1, ecc, fname, n):

    nwalkers, ndim = 32, len(theta)
    p0 = [theta+1e-4*np.random.rand(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
                                    args=(x, y, yerr, M1, ecc))
    p0, lp, state = sampler.run_mcmc(p0, 300)
    sampler.reset()
    p0, lp, state = sampler.run_mcmc(p0, 500)

    fig_labels = ['P', 'M2', 'T0', 'V0', 'omega']
    flatchain = sampler.chain[:, 50:, :].reshape((-1, ndim))
    fig = triangle.corner(flatchain, truths=theta, labels=fig_labels)
    plt.savefig("%s_triangle" % fname)

    print "saving samples"
    f = h5py.File("%s/results/%s_%s_samples" % (DIR, n, fname), "w")
    data = f.create_dataset("samples", np.shape(sampler.chain))
    data[:, :] = np.array(sampler.chain)
    f.close()

    flat = sampler.chain[:, 50:, :].reshape((-1, ndim))
    mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                      zip(*np.percentile(flat, [16, 50, 84], axis=0)))
    np.savetxt("%s/results/%s_%s_results.txt" % (DIR, n, fname), mcmc_result)

    return mcmc_result

if __name__ == "__main__":

    fname = "HD185"

    # load data and params
    DIR = "/Users/angusr/Python/Subgiants/injections"

    N = 10
    for n in range(N):
        x, y, yerr = np.genfromtxt("%s/rv_curves/%s_%s_rvs.txt"
                                   % (DIR, n, fname)).T

        # P, M2, T0, V0, omega
        theta_init = np.genfromtxt("%s/params/%s_%s_params.txt" %
                                   (DIR, n, fname))
        M1, M1_err = np.genfromtxt("%s/params/%s_mass.txt" % (DIR, fname))
        ecc = 0.

        # sample in log period
        theta_init[0] = np.log(theta_init[0])
        results = MCMC(theta_init, x, y, yerr, M1, ecc, fname, n)

        print "true P = ", np.exp(theta_init[0]), "recovered P = ", \
                np.exp(results[0])
        print "true m1 = ", theta_init[1], "recovered m1 = ", results[1]
