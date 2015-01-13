import numpy as np
import matplotlib.pyplot as plt
from model import model
import emcee

# The likelihood function:
def lnlike(theta, t, rv_obs, rverr, M1):
    model_rv = model(theta, t, rverr, M1, wn=False)
    chisq = np.sum(((rv_obs - model_rv) / rverr)**2)
    return -chisq / 2.

# emcee prior and posterior:
def lnprior(theta):
    return 0.

def lnprob(theta, t, rv_obs, rverr, M1):
    return lnprior(theta) + lnlike(theta, t, rv_obs, rverr, M1)

def MCMC():

    nwalkers, ndim = 32, len(theta)
    p0 = [theta+1e-4*np.random.rand(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
    p0, lp, state = sampler.run_mcmc(p0, 300)

    fig_labels = ['P', 'M2', 'T0', 'V0', 'Ecc', 'omega']
    flatchain = sampler.chain[:, 50:, :].reshape((-1, ndim))
    fig = triangle.corner(flatchain, truths=theta, labels=fig_labels)
