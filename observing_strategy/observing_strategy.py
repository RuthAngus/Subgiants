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

mins = 5

def train_BG():
    # Load Beta Gem data
    t, rv, rv_err = BG.rvHJD-BG.rvHJD[0], BG.rv/np.median(BG.rv), \
            BG.rv_err/np.median(BG.rv)
    nm_days = BG.nm*3600*24
    P = 1./nm_days
    l = 38
    # l = len(t)
    t, rv, rv_err = t[:l], rv[:l], rv_err[:l]

    # train GP on Beta Gem (without white noise)
    theta = [15., 5., 1., P]
    k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[3])
    gp = george.GP(k)
    gp.compute(t, rv_err)
    # results = gp.optimize(t, rv, rv_err, dims=[0, 1, 2])[0]
    results = gp.optimize(t, rv, rv_err)[0]
    print np.exp(results)

    # plot BG data and result
    xs = np.linspace(min(t), max(t), 1000)
    mu, cov = gp.predict(rv, xs)
    plt.clf()
    plt.errorbar(t, rv, yerr=rv_err, **reb)
    plt.plot(xs, mu, color=ocols.blue)
    plt.xlabel('$\mathrm{Time~(days)}$')
    plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
    std = np.sqrt(np.diag(cov))
    plt.fill_between(xs, mu+std, mu-std, **fbk)
    plt.savefig('BG_trained')

def sample_prior(theta, xs, yerr, nsamp):
    # Compute GP prior sample
    k = theta[0] * ExpSquaredKernel(theta[1]) * \
            ExpSine2Kernel(theta[2], theta[3])
    gp = george.GP(k)
    gp.compute(xs, yerr)
    np.random.seed(1234)
    samples = gp.sample(xs, nsansamp)

    return samples

def dumb_sampling(fname):

    # Generate time series with a GP
    interval = 6*24  # 10 mins
    interval = 12*24  # 5 mins
    ndays = 10
    xs = np.linspace(0, ndays, interval*ndays) # one point every 10 minutes
    yerr = np.ones_like(xs)*.01

    # Compute GP prior sample
    theta = [8.6969, 1.725e-3, 1.654, P]
    samples = sample_prior(theta, xs, yerr)

    l = interval  # sample every day
    best_time = []
    for i, s in enumerate(samples):

        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(xs, s, color=ocols.blue)
        plt.xlabel("$\mathrm{Time~(days)}$")
        plt.ylabel("$\mathrm{RV~(ms}^{-1}\mathrm{)}$")
        plt.subplot(2, 1, 2)

        times = 48
        rms2 = np.zeros((times, times))
        ms = np.array(range(times)) # from 0 to 120 minute intervals
        for j in ms:
            stds, rms = [], []
            for m in ms:
                xs1, xs2, xs3 = xs[j::l], xs[m+j::l], xs[2*m+j::l]
                ss1, ss2, ss3 = s[j::l], s[m+j::l], s[2*m+j::l]
                yerrs1, yerrs2, yerrs3 = yerr[j::l], yerr[m+j::l], \
                        yerr[2*m+j::l]

                xss = np.vstack((xs1, xs2, xs3))
                ss = np.vstack((ss1, ss2, ss3))
                yerrs = np.vstack((yerrs1, yerrs2, yerrs3))
                ymean = np.mean(ss, axis=0)

                stds.append(np.std(ymean))
                rms.append(np.sqrt(np.mean(ymean**2)))

    #             plt.clf()
    #             plt.plot(xs, s, color=ocols.blue)
    #             plt.plot(xs1, ss1, 'k.')
    #             plt.plot(xs2, ss2, 'k.')
    #             plt.plot(xs3, ss3, 'k.')
    #             plt.plot(xs2, ss2, 'mo')
    #             plt.plot(xs2, np.mean(ss, axis=0), '.', color=ocols.pink)
    #             plt.xlabel('$\mathrm{Time~(days)}$')
    #             plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
    #             plt.savefig('observe')

            stds, rms = np.array(stds), np.array(rms)

            plt.plot(ms*mins, rms, color=ocols.orange, alpha=.3)
            rms2[:][j] = rms

        mean_rms = np.mean(rms2, axis=0)
        ll = mean_rms==min(mean_rms)
        lab = ms[ll][0]*mins
        plt.plot(ms*mins, mean_rms, color=ocols.orange, linewidth=2,
                label='$\mathrm{min~RMS}=%s$' % lab)
        plt.subplots_adjust(hspace=.3)
        plt.xlabel("$\mathrm{Time~between~samples~(minutes)}$")
        plt.ylabel("$\mathrm{RMS}$")
        plt.xlim(0, times*mins)
        plt.legend()

#         fms = sps.find_peaks_cwt(mean_rms, np.arange(1, 100))
#         for fm in fms:
#             print mean_rms[fm]*mins
#             plt.axvline(mean_rms[fm]*mins, color=ocols.blue)

        plt.savefig('GP_%s_%s' % (i, fname))

        best_time.append(lab)

    return rms2, mean_rms, lab, best_time

def obs_times(nmins, ndays, ntests, nsamples):
    t = nmins/60./24  # 10 minutes in days
    times = np.zeros((nsamples*ndays, ntests))
    for i in range(ntests):
        times[:, i] = np.sort(np.concatenate((np.arange(1, ndays-1),
                        np.arange(1-(t*i), ndays-1-(t*i), 1),
                        np.arange(1+(t*i), ndays-1+(t*i), 1))))
    return times

def smart_sampling(fname):

    nmins = 10.  # interval between observations in minutes
    ndays, ntests, nsamples = 12, 10, 3
    xs = obs_times(nmins, ndays, ntests, nsamples)
    yerr = np.ones_like(xs)*.01

    # Compute GP prior sample
    theta = [8.6969, 1.725e-3, 1.654, P]
    nsamp = 1
    samples = sample_prior(theta, xs[0], yerr, nsamp)
    plt.clf()
    plt.plot(xs[0], samples, colour=ocols.blue)
    plt.show()
    raw_input('enter')

    for i, s in enumerate(samples):

        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(xs, s, color=ocols.blue)
        plt.xlabel("$\mathrm{Time~(days)}$")
        plt.ylabel("$\mathrm{RV~(ms}^{-1}\mathrm{)}$")
        plt.subplot(2, 1, 2)

        times = 48
        rms2 = np.zeros((times, times))
        ms = np.array(range(times)) # from 0 to 120 minute intervals
        for j in ms:
            stds, rms = [], []
            for m in ms:
                xs1, xs2, xs3 = xs[j::l], xs[m+j::l], xs[2*m+j::l]
                ss1, ss2, ss3 = s[j::l], s[m+j::l], s[2*m+j::l]
                yerrs1, yerrs2, yerrs3 = yerr[j::l], yerr[m+j::l], \
                        yerr[2*m+j::l]

                xss = np.vstack((xs1, xs2, xs3))
                ss = np.vstack((ss1, ss2, ss3))
                yerrs = np.vstack((yerrs1, yerrs2, yerrs3))
                ymean = np.mean(ss, axis=0)

                stds.append(np.std(ymean))
                rms.append(np.sqrt(np.mean(ymean**2)))

    #             plt.clf()
    #             plt.plot(xs, s, color=ocols.blue)
    #             plt.plot(xs1, ss1, 'k.')
    #             plt.plot(xs2, ss2, 'k.')
    #             plt.plot(xs3, ss3, 'k.')
    #             plt.plot(xs2, ss2, 'mo')
    #             plt.plot(xs2, np.mean(ss, axis=0), '.', color=ocols.pink)
    #             plt.xlabel('$\mathrm{Time~(days)}$')
    #             plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
    #             plt.savefig('observe')

            stds, rms = np.array(stds), np.array(rms)

            plt.plot(ms*mins, rms, color=ocols.orange, alpha=.3)
            rms2[:][j] = rms

        mean_rms = np.mean(rms2, axis=0)
        ll = mean_rms==min(mean_rms)
        lab = ms[ll][0]*mins
        plt.plot(ms*mins, mean_rms, color=ocols.orange, linewidth=2,
                label='$\mathrm{min~RMS}=%s$' % lab)
        plt.subplots_adjust(hspace=.3)
        plt.xlabel("$\mathrm{Time~between~samples~(minutes)}$")
        plt.ylabel("$\mathrm{RMS}$")
        plt.xlim(0, times*mins)
        plt.legend()

#         fms = sps.find_peaks_cwt(mean_rms, np.arange(1, 100))
#         for fm in fms:
#             print mean_rms[fm]*mins
#             plt.axvline(mean_rms[fm]*mins, color=ocols.blue)

        plt.savefig('GP_%s_%s' % (i, fname))

        best_time.append(lab)

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
