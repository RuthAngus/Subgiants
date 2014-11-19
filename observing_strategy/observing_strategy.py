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
plotpar = {'legend.fontsize': 10}
plt.rcParams.update(plotpar)
from scipy import interpolate

def sample_prior(theta, xs, yerr, nsamp):
    # Compute GP prior sample
    k = theta[0] * ExpSquaredKernel(theta[1]) * \
            ExpSine2Kernel(theta[2], theta[3])
    gp = george.GP(k)
    gp.compute(xs, yerr)
    np.random.seed(1234)
    samples = gp.sample(xs, nsamp)
    return samples

def obs_times(nmins, ndays, ntests, nsamples):

    t = nmins/60./24  # 10 minutes in days

    times = np.zeros((nsamples*(ndays-2), ntests))
    for i in range(ntests):
        t1 = np.arange(1, ndays-1)
        t2 = np.arange(1-(t*i), ndays-1-(t*i), 1)
        t3 = np.arange(1+(t*i), ndays-1+(t*i), 1)
        times[:, i] = np.sort(np.concatenate((t1, t2, t3)))
    return times.T

def interp(x, y, times):
    tck = interpolate.splrep(x, y, s=0)
    ynew = interpolate.splev(times, tck, der=0)
    return ynew

def smart_sampling(P, fname):

    nmins = 1.  # interval between observations in minutes
    ndays = 12  # number of nights observed
    ntests = 10  # number of changes in interval
    nsamples = 3  # number of observations per night
    xs = obs_times(nmins, ndays, ntests, nsamples)
    yerr = np.ones_like(xs)*.01

    # Compute GP prior sample
    theta = [8.6969, 1.725e-3, 1.654, P]
    nsamp = 1
    xgrid = np.linspace(0, ndays, 1000)
    samples = sample_prior(theta, xgrid, np.ones_like(xgrid)*.01, nsamp)

    # calculate y values at observation position
    ys = np.zeros_like(xs)
    rms_ys = np.zeros((ntests, ndays-2))
    rms_xs = np.zeros((ntests, ndays-2))
    for i in range(len(xs)):
        ys[i, :] = interp(xgrid, samples, xs[i])

    mean_ys = ys[0][::nsamples]  # because this is the same for all tests
    rms_per_night = np.sqrt(np.mean(mean_ys**2))
    print rms_per_night

def sampling_method(P, fname):
    smart_sampling(P, fname)

if __name__ == "__main__":

    print 'Beta Gem = ', 1./BG.nm/3600.

    P = 1
    sampling_method(P, 'test')
# #     ps = [1, 2, 3]
#     for p in ps:
#         times = []
#         print p
#         P = p/24.
#         times.append(sampling_method(P, 'test'))
#
#         plt.clf()
#         plt.hist(times, color='.3', edgecolor='w')
#         plt.savefig('hist%s' % p)
