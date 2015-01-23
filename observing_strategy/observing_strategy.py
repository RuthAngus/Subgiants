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
from __future__ import print_function

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

# an array of ntests of 2 observing times, separated by nmins, over ndays
# for a given starting time. start = integer
def obs_times2(nmins, ndays, ntests, nsamples, start):
    t = nmins/60./24  # time interval (days)
    times = np.zeros((nsamples*(ndays-2), ntests))  # construct empty array
    # starting point is 'start' intervals before the first day
    st = 1-start*nmins/24./60.
    # calculate observing times
    separations = []
    for i in range(ntests):
        t1 = np.arange(st, ndays-st, 1)  # one obs per day
        t2 = np.arange(st+(t*i), ndays-st+(t*i), 1)  # 3/day, separated by t
        # make sure there are only ndays-2 times so you don't get edge effects
        if len(t1) > ndays-2 or len(t2) > ndays-2:
            t1 = t1[:ndays-2]
            t2 = t2[:ndays-2]
        times[:, i] = np.sort(np.concatenate((t1, t2)))  # sort the times
        separations.append(t*i)
    return times.T, separations

# an array of ntests of 5 observing times, separated by nmins, over ndays
# for a given starting time. start = integer
def obs_times5(nmins, ndays, ntests, nsamples, start):
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
        t4 = np.arange(st-(t*i*2), ndays-st-(t*i), 1)
        t5 = np.arange(st+(t*i*2), ndays-st+(t*i*2), 1)
        # make sure there are only ndays-2 times so you don't get edge effects
        if len(t1) > ndays-2 or len(t2) > ndays-2 or len(t3) > ndays-2 or \
                len(t4) > ndays-2 or len(t5) > ndays-2:
            t1 = t1[:ndays-2]
            t2 = t2[:ndays-2]
            t3 = t3[:ndays-2]
            t4 = t4[:ndays-2]
            t5 = t3[:ndays-2]
        times[:, i] = np.sort(np.concatenate((t1, t2, t3, t4, t5)))
        separations.append(t*i)
    return times.T, separations

# interpolation function
def interp(x, y, times, exptime):
    tck = interpolate.splrep(x, y, s=0)
    yold = interpolate.splev(times, tck, der=0)

    # convert exptime to days
    exptime_days = exptime/24./3600.

    # calculate rvs at every second of exposure
    yexp = np.zeros((exptime, len(times)))
    for i in range(exptime):
        yexp[i, :] = interpolate.splev(times+(i*exptime_days), tck, der=0)

    ynew = np.mean(yexp, axis=0)
    return ynew

# simulate nsamp rv curves with the method defined by stype
def simulate(stype, nsim, fname):
    DIR = '/Users/angusr/Python/Subgiants'
    if stype == "GP":
        # Compute GP prior sample
        theta = [8.6969, 1.725e-3, 1.654, P]
        # points at which to produce the simulated data
        xgrid = np.linspace(0, ndays, 1000)
        return sample_prior(theta, xgrid, np.ones_like(xgrid)*.01, nsim), xgrid
    elif stype == "sine":
        x, y = np.genfromtxt("%s/injections/%s_rvs.txt" % (DIR, fname)).T
        e = 2.
        yerr = np.ones_like(y)*e  # make up uncertainties
#         y += np.random.randn(len(y))*e  # add white noise
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
# start = float between 0 and 1
# output:
# an array of rms values for ntests of nsims
def smart_sampling(nmins, ndays, ntests, nsamples, nsim, start, fname,
                   exptime):

    # generate an array of the observing times
    # ts is an array of ntests, observing times
    if nsamples == 3:
        ts, separations = obs_times(nmins, ndays, ntests, nsamples, start)
    elif nsamples == 2:
        ts, separations = obs_times2(nmins, ndays, ntests, nsamples, start)
    elif nsamples == 5:
        ts, separations = obs_times5(nmins, ndays, ntests, nsamples, start)
    yerr = np.ones_like(ts)*.01

    # samples is a 2d array containing simulated data generated at xgrid times
    samples, xgrid = simulate("sine", nsim, fname)

    # calculate y values at observation positions
    # number of observations, number of tests, number of simulations
    ys = np.ndarray((len(ts.T), ntests, nsim))
    for i in range(nsim):
        for j in range(ntests):
            # xs you have, ys, xs you want
            ys[:, j, i] = interp(xgrid, samples[:, i], ts[j, :], exptime)

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
def all_start_times(start_times, nmins, ndays, ntests, nsamples, nsim, fname,
                    exptime):

    all_rms = np.zeros((len(start_times), ntests))
    for i, st in enumerate(start_times):
        print(i, "of ", len(start_times))

        # calculate rms over all observing separations ntests, nsim
        rms, separations, samples, xgrid = smart_sampling(nmins, ndays, ntests,
                                                          nsamples, nsim, st,
                                                          fname, exptime)
        s = np.array(separations) *24*60
        all_rms[i, :] = rms[:, 0]

    # calculate mean and minimum rms
    mean_rms = np.mean(all_rms, axis=0)
    l = mean_rms == min(mean_rms)

    return all_rms, s, mean_rms, s[l][0], samples, xgrid

def os(nmins, ndays, ntests, nsamples, nsim, exptime, fnames):

    times, min_rms = [], []
    for i, fname in enumerate(fnames):
        print(i, fname)

        # range of start times. Go from 0 to the number of start times in a day
        start_times = np.arange(0, 24*60/nmins/10, 1)

        # calculate rms (smart sampling)
        all_rms, s, mean_rms, best_time, samples, xgrid = \
                all_start_times(start_times, nmins, ndays, ntests, nsamples,
                                nsim, str(fname), exptime)

        times.append(best_time)
        min_rms.append(min(mean_rms))

        # plot
        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(xgrid, samples[:, 0], color=ocols.blue)
        plt.ylabel("$\mathrm{RV}~(ms^{-1})$")
        plt.xlabel("$\mathrm{Time~(days)}$")

        plt.subplot(2, 1, 2)
        for i in range(len(start_times)):
            plt.plot(s, all_rms[i, :], color=ocols.orange, alpha=.2)
        plt.plot(s, mean_rms, color=ocols.orange, linewidth=2,
                 label="$\mathrm{Minimum}=%s$" % best_time)
        plt.ylabel("$\mathrm{RMS}~(ms^{-1})$")
        plt.xlabel("$\mathrm{Interval,}~\Delta t~\mathrm{(mins)}$")
        plt.legend(loc="best")
#         plt.subplots_adjust("hspace=1")
        plt.savefig("%s_%s_os.pdf" % (fname, nsamples))

    fnames = np.array([int(filter(str.isdigit, fname)) for fname in fnames])
    np.savetxt("named_best_times_%s.txt" % nsamples,
               np.transpose((fnames, times, min_rms)))
#     np.savetxt("AMP_best_times.txt", np.transpose((fnames, times)))

if __name__ == "__main__":

    DIR = '/Users/angusr/Python/Subgiants'

    nmins = 2  # interval between observations in minutes
    ndays = 10  # number of nights observed. if this is greater than the number
    # of nights in the simulated data, the rms will increase
    ntests = 100  # number of changes in interval
    nsamples = 2  # number of observations per night
    nsim = 1  # number of simulations
    exptime = 100

#     fnames = [1430163, 1435467, 1725815, 2010607, 2309595, 3424541, 5955122,
#               7747078, 7976303, 8026226, 8524425, 10018963, 11026764]
#     fnames = [3424541, 5955122, 7747078, 7976303, 8026226, 8524425, 10018963,
#               11026764]

    # luan's subgiants
    data = np.genfromtxt("%s/proposal/sample_luan.out" % DIR,
                         skip_header=1, dtype=str).T
    fnames = data[0]

    os(nmins, ndays, ntests, nsamples, nsim, exptime, fnames)
