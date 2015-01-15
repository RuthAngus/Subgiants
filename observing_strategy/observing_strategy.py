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
from sampling import dumb_sampling, sample_prior, sc_sampling
from sine_wave_gen import kepler_sine_synth, HDsine_synth

# an array of ntests of 3 observing times, separated by nmins, over ndays
def obs_times(nmins, ndays, ntests, nsamples):
    t = nmins/60./24  # time interval (days)
    times = np.zeros((nsamples*(ndays-2), ntests))  # construct empty array
    for i in range(ntests):
        t1 = np.arange(1, ndays-1)  # one obs per day
        t2 = np.arange(1-(t*i), ndays-1-(t*i), 1)  # 2 per day, separated by t
        t3 = np.arange(1+(t*i), ndays-1+(t*i), 1)  # 3 per day, separated by t
        times[:, i] = np.sort(np.concatenate((t1, t2, t3)))  # sort the times
    return times.T

# interpolation function
def interp(x, y, times):
    tck = interpolate.splrep(x, y, s=0)
    ynew = interpolate.splev(times, tck, der=0)
    return ynew

# simulate nsamp rv curves with the method defined by stype
def simulate(stype, nsim):
    if stype == "GP":
        # Compute GP prior sample
        theta = [8.6969, 1.725e-3, 1.654, P]
        # points at which to produce the simulated data
        xgrid = np.linspace(0, ndays, 1000)
        return sample_prior(theta, xgrid, np.ones_like(xgrid)*.01, nsim), xgrid
    elif stype == "sine":
        x, y = np.genfromtxt("%s/injections/HD185_rvs.txt" % DIR).T
        yerr = np.ones_like(y)*2.  # make up uncertainties
        samples = np.zeros((len(x), nsim))
        for i in range(nsim):
            samples[:, i] = y
        return samples, x

# sample simulated data at arbitrary observation times
def smart_sampling(nmins, ndays, ntests, nsamples, nsim, fname):

    print "generate an array of the observing times"
    # ts is an array of ntests, observing times
    ts = obs_times(nmins, ndays, ntests, nsamples)
    yerr = np.ones_like(ts)*.01

    # samples is a 2d array containing simulated data generated at xgrid times
    samples, xgrid = simulate("sine", nsim)

    print "calculate y values at observation positions"
    # number of observations, number of tests, number of simulations
    ys = np.ndarray((len(ts.T), ntests, nsim))
    for i in range(nsim):
        for j in range(ntests):
            # xs you have, ys, xs you want
            ys[:, j, i] = interp(xgrid, samples[:, i], ts[j, :])

#     plt.clf()
#     plt.plot(xgrid, samples[:, 0])
#     plt.plot(ts[j, :], ys[:, j, i], "r.")
#     plt.savefig('test')

    # store just the single observation values
    # don't loop because this is the same for all tests

#     central_ys = ys[:, 0, 0]
#     central_ts = ts[0, :]
#     print "rms = ", np.sqrt(np.mean(central_ys[::3]**2))

    # total number of obs, ntests, nsims
        for i in range(ntests):
            nightly_obs = np.reshape(ys[:, i, 0], (len(ys[:, i, 0])/nsamples,
                                     nsamples))
            nightly_av = np.mean(nightly_obs, axis=1)
            rms = np.sqrt(np.mean(nightly_av**2))
            print "rms = ", rms
            raw_input('enter')

    # calculate the rms for just the nightly observations
    rms_per_night = np.sqrt(np.mean(mean_ys**2))

    print rms_per_night

if __name__ == "__main__":

    DIR = '/Users/angusr/Python/Subgiants'

    nmins = 5  # interval between observations in minutes
    ndays = 10  # number of nights observed
    ntests = 10  # number of changes in interval
    nsamples = 3  # number of observations per night
    nsim = 2  # number of simulations
    fname = "test"

    smart_sampling(nmins, ndays, ntests, nsamples, nsim, fname)

#     P, nsamp, fname = 1, 1, "HD185"
#     nmins, ndays = 5, 10  # number of minutes, ndays
#     sample_type = 'sine'
#     interval = (60./nmins)*24  # samples per day
#     xs = np.linspace(0, ndays, interval*ndays) # one point every nmins

    # generate synthetic data
#     HDsine_synth(xs, ndays, train=True)
#     HDsine_synth(xs, ndays, train=False)

#     kid = "9002278"
#     kepler_sine_synth(kid, xs, ndays, train=True)

#     dumb_sampling(P, nsamp, nmins, ndays, sample_type, fname)

#     x, y = np.genfromtxt('%s/data/hd185351.q16sc.ts' % DIR).T
#     y_err = np.ones_like(y)*6.255e-5
#     fname = 'hd'
#     sc_sampling(fname)
