import numpy as np
import matplotlib.pyplot as plt
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel, WhiteKernel
from colours import plot_colours
ocols = plot_colours()
from rc_params import plot_params
reb, fbk = plot_params()
from BGdata import BetaGem
BG = BetaGem()
from scaling_relations import nu_max, delta_nu
import scipy.signal as sps
from sampling import dumb_sampling

plotpar = {'legend.fontsize': 10}
plt.rcParams.update(plotpar)

mins = 5

def sample_prior(theta, xs, yerr, nsamp):
    # Compute GP prior sample
    k = theta[0] * ExpSquaredKernel(theta[1]) * \
            ExpSine2Kernel(theta[2], theta[3])
    gp = george.GP(k)
    gp.compute(xs, yerr)
    np.random.seed(1234)
    samples = gp.sample(xs, nsamp)
    return samples

def smart_sampling(fname):

    nmins = 50.  # interval between observations in minutes
    ndays, ntests, nsamples = 12, 10, 3
    xs = obs_times(nmins, ndays, ntests, nsamples)
    yerr = np.ones_like(xs)*.01

    # Compute GP prior sample
    theta = [8.6969, 1.725e-3, 1.654, P]
    nsamp = 1
    samples = sample_prior(theta, xs[0], yerr[0], nsamp)
    plt.clf()
    plt.plot(xs[0], samples, '.')
#     xss = np.linspace(min(xs[0]), max(xs[0]), 1000)
#     samples = sample_prior(theta, xss, np.ones_like(xss)*.1, nsamp)
#     plt.plot(xss, samples, 'b')
    plt.plot(xs[1], samples, '.')
    plt.plot(xs[2], samples, '.')
    plt.plot(xs[3], samples, '.')
    plt.plot(xs[4], samples, '.')
    plt.plot(xs[5], samples, '.')
    plt.plot(xs[6], samples, '.')
    plt.plot(xs[7], samples, '.')
    plt.show()
    raw_input('enter')

    return rms2, mean_rms, lab, best_time


def sampling_method(P, fname):
    rms2, mean_rms, lab, best_time = smart_sampling(fname)
    return np.array(best_time)

if __name__ == "__main__":

    print 'Beta Gem = ', 1./BG.nm/3600.

#     ps = [1, 2, 3]
    ps = [2]
    for p in ps:
        times = []
        print p
        P = p/24.
        times.append(sampling_method(P, 'test'))

        plt.clf()
        plt.hist(times, color='.3', edgecolor='w')
        plt.savefig('hist%s' % p)
