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

def model3_freq(par, x, freq):
    # example with 3 sine waves
    return par[0]*np.sin(2*np.pi*freq[0]*x + par[3]) + \
            par[1]*np.sin(2*np.pi*freq[1]*x + par[4]) + \
            par[2]*np.sin(2*np.pi*freq[2]*x + par[5])

def model_freq(par, x, freq):
    sin_sum, n = 0, len(freq)
    for i in range(n):
        sin_sum += par[i]*np.sin(2*np.pi*freq[i]*x + par[i+n])
    return sin_sum

def lnlike_freq(par, x, y, yerr, freq):
    return np.sum(-0.5*(y - model_freq(par, x, freq))**2/yerr**2)

def lnlike3_freq(par, x, y, yerr, freq):
    return np.sum(-0.5*(y - model3_freq(par, x, freq))**2/yerr**2)

def lnlike(par, x, y, yerr):
    return np.sum(-0.5*(y - model3(par, x))**2/yerr**2)

def lnprior3(par):
    if -10 < par[0] < 10 and -10 < par[1] < 10 \
            and -10 < par[2] < 10 and -10 < par[3] < 10 \
            and -10 < par[4] < 10 and -10 < par[5] < 10 \
            and -10 < par[6] < 10 and -10 < par[7] < 10:
        return 0.
    else: return -np.inf

def lnprior3_freq(par):
    if -10 < par[0] < 10 and -10 < par[1] < 10 \
            and -10 < par[2] < 10 and -10 < par[3] < 10 \
            and -10 < par[4] < 10:
        return 0.
    else: return -np.inf

def lnprior4_freq(par):
    if -10 < par[0] < 10 and -10 < par[1] < 10 \
            and -10 < par[2] < 10 and -10 < par[3] < 10 \
            and -10 < par[4] < 10 and -10 < par[5] < 10 \
            and -10 < par[6] < 10:
        return 0.
    else: return -np.inf

def lnprior5_freq(par):
    if -10 < par[0] < 10 and -10 < par[1] < 10 \
            and -10 < par[2] < 10 and -10 < par[3] < 10 \
            and -10 < par[4] < 10 and -10 < par[5] < 10 \
            and -10 < par[6] < 10 and -10 < par[7] < 10 \
            and -10 < par[8] < 10:
        return 0.
    else: return -np.inf

def lnprob(par, x, y, yerr):
    return lnlike(par, x, y, yerr) + lnprior3(par)

def lnprob_freq(par, x, y, yerr, freq):
    return lnlike_freq(par, x, y, yerr, freq) + lnprior_freq(par)

def lnprob3_freq(par, x, y, yerr, freq):
    return lnlike3_freq(par, x, y, yerr, freq) + lnprior3_freq(par)

def lnprob4_freq(par, x, y, yerr, freq):
    return lnlike_freq(par, x, y, yerr, freq) + lnprior4_freq(par)

def lnprob5_freq(par, x, y, yerr, freq):
    return lnlike_freq(par, x, y, yerr, freq) + lnprior5_freq(par)

def MCMC(par_init, args, lnlike, lnprob, lnprior):

    ndim, nwalkers = len(par_init), 100
    p0 = [par_init + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=args)
    print 'burning in...'
    p0, lp, state = sampler.run_mcmc(p0, 1000)
    sampler.reset()
    print 'production run...'
    p0, lp, state = sampler.run_mcmc(p0, 5000)

    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
#     fig_labels = ["a1", "a2", "a3", "phi1", "phi2", "phi3"]
#     fig_labels = ["a1", "a2", "a3", "a4", "a5", "phi1",
#                   "phi2", "phi3", "phi4", "phi5", "1", "2", "3"]
    fig_labels = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
    fig = triangle.corner(samples, labels=fig_labels,
                          truths=par_init)
    fig.savefig("triangle")

    mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                      zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    mres = np.array(mcmc_result)[:, 0]
    return mres

if __name__ == "__main__":

    # generate some data
    par_true3 = [1., 2., 3., 1., 2., 3., 0., 0., 0.]
    par_true3_freq = [1., 2., 3., 0., 0., 0.]
    freqs = [1., 2., 3.]
    par_true1 = [1., 1., 0.]
    par_true = par_true3_freq

    x = np.arange(0, 10, 0.01)
#     y = model3(par_true, x)
    y = model3_freq(par_true3_freq, x, freqs)
    yerr = np.ones_like(x)*0.01
    y += np.random.rand(len(y)) * yerr

    par_init = par_true + .01*np.random.rand(len(par_true))

#     plt.clf()
#     plt.errorbar(x, y, yerr=yerr, capsize=0, ecolor='.8', fmt='k.')
#     plt.plot(x, model3_freq(par_true3_freq, x, freqs))
#     plt.show()

#     args = (x, y, yerr)
    args = (x, y, yerr, freqs)

#     MCMC(par_init, args, lnlike, lnprob, lnprior1)
    mres = MCMC(par_init, args, lnlike3_freq, lnprob3_freq, lnprior3_freq)
    print mres
