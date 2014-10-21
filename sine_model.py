import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin
import emcee
import triangle

# take arrays of freqencies, amplitudes and phis
# to fit a series of sine waves to your data
def model3(par, x):
    # example with 3 sine waves
    return par[3]*np.sin(2*np.pi*par[0]*x + par[6]) + \
            par[4]*np.sin(2*np.pi*par[1]*x + par[7]) + \
            par[5]*np.sin(2*np.pi*par[2]*x + par[8])

def model1(par, x):
    # example with 1 sine wave
    return par[1]*np.sin(2*np.pi*par[0]*x + par[2])

def model_freq(par, x, freq):
    # example with 1 sine wave with known freq
    return par[0]*np.sin(2*np.pi*freq*x + par[1])

def neglnlike(par, x, y, yerr):
    return -np.logaddexp.reduce(-0.5*(y - model3(par, x))**2/yerr**2, axis=0)

def neglnlike_freq(par, x, y, yerr):
    return -np.logaddexp.reduce(-0.5*(y - model_freq(par, x))**2/yerr**2, axis=0)

def lnlike(par, x, y, yerr):
    return np.sum(-0.5*(y - model3(par, x))**2/yerr**2)

def lnprior1(par):
    if -10 < par[0] < 10 and -10 < par[1] < 10 \
            and -10 < par[2] < 10:
        return 0.
    else: return -np.inf

def lnprior3(par):
    if -10 < par[0] < 10 and -10 < par[1] < 10 \
            and -10 < par[2] < 10 and -10 < par[3] < 10 \
            and -10 < par[4] < 10 and -10 < par[5] < 10 \
            and -10 < par[6] < 10 and -10 < par[7] < 10:
        return 0.
    else: return -np.inf

def lnprob(par, x, y, yerr):
    return lnlike(par, x, y, yerr) + lnprior3(par)

def MCMC(par_init):
    ndim, nwalkers = len(par_init), 32
    p0 = [par_init + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
    p0, lp, state = sampler.run_mcmc(p0, 1000)
    sampler.reset()
    p0, lp, state = sampler.run_mcmc(p0, 2000)

    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    fig_labels = ["f", "a", "phi", "f2", "a2", "phi2", "f3", "a3", "phi3"]
    fig = triangle.corner(samples, labels=fig_labels,
                          truths=par_true)
    fig.savefig("triangle")

if __name__ == "__main__":

    # generate some data
    par_true3 = [1., 2., 3., 1., 2., 3., 0., 0., 0.]
    par_true1 = [1., 1., 0.]
    par_true = par_true3

    x = np.arange(0, 10, 0.01)
    y = model3(par_true, x)
    yerr = np.ones_like(x)*0.01
    y += np.random.rand(len(y)) * yerr

    par_init = par_true + .01*np.random.rand(len(par_true))

    args = (x, y, yerr)
    results = fmin(neglnlike, par_init, args=args)
    print par_true
    print par_init

    plt.clf()
    plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
    plt.plot(x, model3(par_init, x), 'g')
    plt.show()

    MCMC(par_init)
