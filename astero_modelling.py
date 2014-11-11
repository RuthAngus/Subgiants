import numpy as np
import matplotlib.pyplot as plt
from scaling_relations import nu_max, delta_nu
from sin_tests import fit_sine
import emcee
import triangle
import h5py

def gen_freqs(m, r, t, nfreqs):
    nm = nu_max(m, r, t)
    dn = delta_nu(m, r)
    return np.arange(nm-nfreqs*dn, nm+nfreqs, dn)

def model(pars, x, y, yerr, nfreqs):
    freqs = gen_freqs(pars[0], pars[1], pars[2], nfreqs)
    return fit_sine(x, y, 2*np.pi*freqs)

def lnlike(pars, x, y, yerr, nfreqs):
    return np.sum(-0.5*(y - model(pars, x, y, yerr, nfreqs))**2/yerr**2)

def lnprior(pars):
    m, r, t = pars
    if 0. < m < 10. and 0. < r < 10. and 500. < t < 20000:
        return 0.
    else:
        return -np.inf

def Gaussian_priors(pars):
    m, r, t, m_sig, r_sig, t_sig = pars
    return - m**2/(2*m_sig**2) - r**2/(2*r_sig) - t**2/(2*t_sig)

def lnprob(pars, x, y, yerr, nfreqs, like, prior):
    return like(pars, x, y, yerr, nfreqs) + prior(pars)

def MCMC(par_init, args, burnin, runs, fname, fig_labels):

    x, y, yerr, nfreqs, like, prior = args
    print 'initial likelihood = ', like(par_init, x, y, yerr, nfreqs)
    ndim, nwalkers = len(par_init), 100
    p0 = [par_init + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=args)
    print "burning in..."
    p0, lp, state = sampler.run_mcmc(p0, burnin)
    sampler.reset()
    print "running..."
    p0, lp, state = sampler.run_mcmc(p0, runs)

    samples = sampler.chain[:, 50:, :]
    flatchain = samples.reshape((-1, ndim))
    mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                      zip(*np.percentile(flatchain, [16, 50, 84], axis=0)))
    mres = np.array(mcmc_result)[:, 0]

    fig = triangle.corner(flatchain, labels=fig_labels, truths=mres)
    fig.savefig("triangle%s" % fname)

    print("Plotting traces")
    plt.figure()
    for i in range(ndim):
        plt.clf()
        plt.plot(samples[:, :, i].T, 'k-', alpha=0.3)
        plt.savefig("%s%s.png" %(i, fname))

    print "saving samples"
    f = h5py.File("samples_%s" % fname, "w")
    data = f.create_dataset("samples", np.shape(sampler.chain))
    data[:,:] = np.array(sampler.chain)
    f.close()

    return mres, mcmc_result
