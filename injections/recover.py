import numpy as np
import matplotlib.pyplot as plt
from model import model
import emcee

# The likelihood function:
def lnlike(theta, t, rv_obs, rverr, M1, ecc):
    model_rv = model(theta, t, rverr, M1, wn=False)
    chisq = np.sum(((rv_obs - model_rv) / rverr)**2)
    return -chisq / 2.

# emcee prior and posterior:
def lnprior(theta):
    return 0.

def lnprob(theta, t, rv_obs, rverr, M1, ecc):
    return lnprior(theta) + lnlike(theta, t, rv_obs, rverr, M1, ecc)

def MCMC(theta, x, y, yerr, M1, ecc):

    nwalkers, ndim = 32, len(theta)
    p0 = [theta+1e-4*np.random.rand(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
                                    args=(x, y, yerr, M1, ecc))
    p0, lp, state = sampler.run_mcmc(p0, 300)

    fig_labels = ['P', 'M2', 'T0', 'V0', 'omega']
    flatchain = sampler.chain[:, 50:, :].reshape((-1, ndim))
    fig = triangle.corner(flatchain, truths=theta, labels=fig_labels)
    plt.savefig("triangle")

if __name__ == "__main__":

    # load data and params
    DIR = "/Users/angusr/Python/Subgiants/injections"
    x, y, yerr = np.genfromtxt("%s/rv_curves/0_HD185_rvs.txt" % DIR).T
    P, M2 = np.genfromtxt("%s/params/0_HD185_params.txt" % DIR)
    M1 = np.genfromtxt("%s/params/HD185_mass.txt" % DIR)

    raw_input('enter')

    # P, M2, T0, V0, omega
    theta_init = [P, M2, 0., 0., 0.]
    ecc = 0.

    MCMC(theta_init, x, y, yerr, M1, ecc)
