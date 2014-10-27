import numpy as np
import matplotlib.pyplot as plt
from scaling_relations import nu_max, delta_nu
from sin_tests import fit_sine
import emcee
import triangle

def gen_freqs(m, r, t, nfreqs):
    nm = nu_max(m, r, t)
    dn = delta_nu(m, r)
    return np.arange(nm-nfreqs*dn, nm+nfreqs, dn)

def model(pars, x, y, yerr, nfreqs):
    m, r, t = pars
    freqs = gen_freqs(m, r, t, nfreqs)
    return fit_sine(x, y, 2*np.pi*freqs)

def lnlike(pars, x, y, yerr, nfreqs):
    m, r, teff = pars
    return np.sum(-0.5*(y - model(pars, x, y, yerr, nfreqs))**2/yerr**2)

def lnprior(pars):
    m, r, t = pars
    if 0. < m < 10. and 0. < r < 10. and 500. < t < 20000:
        return 0.
    else:
        return -np.inf

def lnprob(pars, x, y, yerr, nfreqs):
    return lnlike(pars, x, y, yerr, nfreqs) + lnprior(pars)

def MCMC(par_init, args, burnin, runs, fname, fig_labels):

    x, y, yerr, nfreqs = args
    print 'initial likelihood = ', lnlike(par_init, x, y, yerr, nfreqs)
    ndim, nwalkers = len(par_init), 32
    p0 = [par_init + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=args)
    print "burning in..."
    p0, lp, state = sampler.run_mcmc(p0, burnin)
    sampler.reset()
    print "running..."
    p0, lp, state = sampler.run_mcmc(p0, runs)

    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    fig = triangle.corner(samples, labels=fig_labels,
                          truths=par_init)
    fig.savefig("triangle%s" % fname)

    mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                      zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    mres = np.array(mcmc_result)[:, 0]
    return mres
