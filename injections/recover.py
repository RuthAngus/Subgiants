import numpy as np
import matplotlib.pyplot as plt
from model import model
import emcee
import triangle
import h5py

# The likelihood function:
def lnlike(theta, t, rv_obs, rverr, M1, ecc):
    model_rv = model(theta, t, rverr, M1, ecc, wn=False)
    chisq = np.sum(((rv_obs - model_rv) / rverr)**2)
    return -chisq / 2.

# emcee prior and posterior:
def lnprior(theta):
    return 0.

def lnprob(theta, t, rv_obs, rverr, M1, ecc):
    return lnprior(theta) + lnlike(theta, t, rv_obs, rverr, M1, ecc)

def MCMC(theta, x, y, yerr, M1, ecc, fname):

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
    f = h5py.File("%s_samples" % fname, "w")
    data = f.create_dataset("samples", np.shape(sampler.chain))
    data[:, :] = np.array(sampler.chain)
    f.close()

    flat = sampler.chain[:, 50:, :].reshape((-1, ndim))
    mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                      zip(*np.percentile(flat, [16, 50, 84], axis=0)))
    mres = np.array(mcmc_result)[:, 0]
    print 'mcmc_result = ', mres

if __name__ == "__main__":

    fname = "HD185"

    # load data and params
    DIR = "/Users/angusr/Python/Subgiants/injections"
    x, y, yerr = np.genfromtxt("%s/rv_curves/0_%s_rvs.txt" % (DIR, fname)).T
    P, M2 = np.genfromtxt("%s/params/0_%s_params.txt" % (DIR, fname))
    M1, M1_err = np.genfromtxt("%s/params/%s_mass.txt" % (DIR, fname))

    # P, M2, T0, V0, omega
    theta_init = [P, M2, 0., 0., 0.]
    ecc = 0.

    MCMC(theta_init, x, y, yerr, M1, ecc, fname)
