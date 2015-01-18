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
# for a given starting time. start = integer
def obs_times(nmins, ndays, ntests, nsamples, start):

    t = nmins/60./24  # time interval (days)
    times = np.zeros((nsamples*(ndays-2), ntests))  # construct empty array

    # starting point is 'start' intervals before the first day
    st = 1-start*nmins/24./60.

    # calculate observing times
    separations = []
    for i in range(ntests):
        t1 = np.arange(st, ndays-st, 1)  # one obs per day
        t2 = np.arange(st-(t*i), ndays-st-(t*i), 1)  # 2/day, separated by t
        t3 = np.arange(st+(t*i), ndays-st+(t*i), 1)  # 3/day, separated by t

        # make sure there are only ndays-2 times so you don't get edge effects
        if len(t1) > ndays-2 or len(t2) > ndays-2 or len(t3) > ndays-2:
            t1 = t1[:ndays-2]
            t2 = t2[:ndays-2]
            t3 = t3[:ndays-2]

        times[:, i] = np.sort(np.concatenate((t1, t2, t3)))  # sort the times
        separations.append(t*i)
    return times.T, separations

# interpolation function
def interp(x, y, times):
    tck = interpolate.splrep(x, y, s=0)
    ynew = interpolate.splev(times, tck, der=0)
    return ynew

# simulate nsamp rv curves with the method defined by stype
def simulate(stype, nsim, fname):
    if stype == "GP":
        # Compute GP prior sample
        theta = [8.6969, 1.725e-3, 1.654, P]
        # points at which to produce the simulated data
        xgrid = np.linspace(0, ndays, 1000)
        return sample_prior(theta, xgrid, np.ones_like(xgrid)*.01, nsim), xgrid
    elif stype == "sine":
        x, y = np.genfromtxt("%s/injections/%s_rvs.txt" % (DIR, fname)).T
        yerr = np.ones_like(y)*2.  # make up uncertainties
        samples = np.zeros((len(x), nsim))
        for i in range(nsim):
            samples[:, i] = y
        return samples, x

# sample simulated data at arbitrary observation times
# input:
# nmins = interval between observations in minutes
# ndays = number of nights observed
# ntests = number of changes in interval
# nsamples = number of observations per night
# nsim = number of simulations
# start = starting position, int. should be range(0, 24*60/nmins)
# output:
# an array of rms values for ntests of nsims
def smart_sampling(nmins, ndays, ntests, nsamples, nsim, start, fname):

    # generate an array of the observing times
    # ts is an array of ntests, observing times
    ts, separations = obs_times(nmins, ndays, ntests, nsamples, start)
    yerr = np.ones_like(ts)*.01

    # samples is a 2d array containing simulated data generated at xgrid times
    samples, xgrid = simulate("sine", nsim, fname)

    # calculate y values at observation positions
    # number of observations, number of tests, number of simulations
    ys = np.ndarray((len(ts.T), ntests, nsim))
    for i in range(nsim):
        for j in range(ntests):
            # xs you have, ys, xs you want
            ys[:, j, i] = interp(xgrid, samples[:, i], ts[j, :])

    # calculate rms
    # ys = total number of obs, ntests, nsims
    rms = np.zeros((ntests, nsim))
    for j in range(nsim):
        for i in range(ntests):
            nightly_obs = np.reshape(ys[:, i, j], (len(ys[:, i, j])/nsamples,
                                     nsamples))
            nightly_av = np.mean(nightly_obs, axis=1)
            rms[i, j] = np.sqrt(np.mean(nightly_av**2))

    return rms, separations, samples, xgrid

# calculate rms over all possible starting times
def all_start_times(start_times, nmins, ndays, ntests, nsamples, nsim, fname):

    all_rms = np.zeros((len(start_times), ntests))
    for i, st in enumerate(start_times):
        print i, "of ", len(start_times)

        # calculate rms over all observing separations ntests, nsim
        rms, separations, samples, xgrid = smart_sampling(nmins, ndays, ntests,
                                                          nsamples, nsim, st,
                                                          fname)
        s = np.array(separations) *24*60
        all_rms[i, :] = rms[:, 0]

    # calculate mean and minimum rms
    mean_rms = np.mean(all_rms, axis=0)
    l = mean_rms == min(mean_rms)

    return all_rms, s, mean_rms, l, samples, xgrid

if __name__ == "__main__":

    DIR = '/Users/angusr/Python/Subgiants'

    nmins = 1  # interval between observations in minutes
    ndays = 4  # number of nights observed. if this is greater than the number
    # of nights in the simulated data, the rms will increase
    ntests = 100  # number of changes in interval
    nsamples = 3  # number of observations per night
    nsim = 1  # number of simulations

#     fname = "HD185"
    fname = "3424541"

    # range of start times. Go from 0 to the number of start times in a day
    start_times = np.arange(0, 24*60/nmins/40, 1)

    # calculate rms (smart sampling)
    all_rms, s, mean_rms, l, samples, xgrid = \
            all_start_times(start_times, nmins, ndays, ntests, nsamples, nsim,
                            fname)

    # plot
    plt.clf()
    plt.subplot(2, 1, 1)
    plt.plot(xgrid, samples[:, 0], color=ocols.blue)

    plt.subplot(2, 1, 2)
    for i in range(len(start_times)):
        plt.plot(s, all_rms[i, :], color=ocols.orange, alpha=.05)
    plt.plot(mean_rms, color=ocols.orange, linewidth=2,
             label="Minimum=%s" % s[l])
    plt.legend(loc="best")
    plt.savefig("%s_os" % fname)

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
