import numpy as np
import matplotlib.pyplot as plt
from scaling_relations import nu_max, delta_nu
from sin_tests import fit_sine_err, fit_sine
import emcee
import triangle
import h5py
# import isochrone_calcs

def gen_freqs(m, r, t, nfreqs):
    nm = nu_max(m, r, t)
    dn = delta_nu(m, r)
    return np.arange(nm-nfreqs*dn, nm+nfreqs, dn)

def gen_freqs_alt(logg, rho, t, nfreqs):
    nm = nu_max_alt(logg, rho, t)
    dn = delta_nu_alt(logg, rho)
    return np.arange(nm-nfreqs*dn, nm+nfreqs, dn)

# def iso_gen_freqs(m, t, nfreqs):
#     r = isochrone_calcs.rad(m, t)
#     nm = nu_max(m, r, t)
#     dn = delta_nu(m, r)
#     return np.arange(nm-nfreqs*dn, nm+nfreqs, dn)

def model(pars, x, y, yerr, nfreqs):
    freqs = gen_freqs(pars[0], pars[1], pars[2], nfreqs)
    ys, A = fit_sine_err(x, y, yerr, 2*np.pi*freqs)
# #     ys, A = fit_sine(x, y, 2*np.pi*freqs)
    return ys
#     return np.ones_like(y)

def lnlike(pars, x, y, yerr, nfreqs):
    return np.sum(-0.5*(y - model(pars, x, y, yerr, nfreqs))**2/yerr**2)

def lnprior(pars):
    m, r, t = pars
    if 0. < m < 10. and 0. < r < 10. and 500. < t < 20000:
        return 0.
    else:
        return -np.inf

def HD21072_priors(pars, sigs):
    return - pars[0]**2/(2*sigs[0]**2) - pars[1]**2/(2*sigs[1]) - \
            pars[2]**2/(2*sigs[2])

def lnprob(pars, x, y, yerr, nfreqs, like, prior, sigs):
    return like(pars, x, y, yerr, nfreqs) + prior(pars, sigs)

def MCMC(par_init, args, burnin, runs, fname, fig_labels):

    x, y, yerr, nfreqs, like, prior, sigs = args
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
