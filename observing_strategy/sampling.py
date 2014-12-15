import numpy as np
import matplotlib.pyplot as plt
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel, WhiteKernel
from colours import plot_colours
ocols = plot_colours()
from rc_params import plot_params
reb, fbk = plot_params()
from sine_wave_gen import sine_synth
from fit_gaussians import GP_mix

def sample_prior(theta, xs, yerr, nsamp):
    # Compute GP prior sample
    k = theta[0] * ExpSquaredKernel(theta[1]) * \
            ExpSine2Kernel(theta[2], theta[3])
    gp = george.GP(k)
    gp.compute(xs, yerr)
    np.random.seed(1234)
    samples = gp.sample(xs, nsamp)
    return samples

def sc_sampling(fname, times=100):

    from sc_target import flux_rv, hd185_rv
    DIR = '/Users/angusr/Python/Subgiants'
    x, y = np.genfromtxt('%s/data/hd185351.q16sc.ts' % DIR).T
    y_err = np.ones_like(y)*6.255e-5
    teff = 5042
    t_err = 0
    # convert flux to rvs
    rv, rv_err, dlL = flux_rv(y, y_err, teff, t_err)
    mins = 5

    best_time = []
    for i in range(10):
        print i

        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(x, rv, color=ocols.blue)
        plt.xlabel("$\mathrm{Time~(days)}$")
        plt.ylabel("$\mathrm{RV~(ms}^{-1}\mathrm{)}$")
        plt.subplot(2, 1, 2)

        rms2 = np.zeros((times, times))
        ms = np.array(range(times))
        js = np.linspace(0.01, 1, times)
        for j in ms:  # once for every test
            stds, rms = [], []
            for m in ms:  # once for every starting pos
                np.random.seed(ms)
                xs1, ss1, yerrs1 = hd185_rv(x, rv, rv_err, 10, -mins*m, j)
                np.random.seed(ms)
                xs2, ss2, yerrs2 = hd185_rv(x, rv, rv_err, 10, 0, j)
                np.random.seed(ms)
                xs3, ss3, yerrs3 = hd185_rv(x, rv, rv_err, 10, +mins*m, j)

                xss = np.vstack((xs1, xs2, xs3))
                ss = np.vstack((ss1, ss2, ss3))
                yerrs = np.vstack((yerrs1, yerrs2, yerrs3))
                ymean = np.mean(ss, axis=0)  # take the mean of the 3

                stds.append(np.std(ymean))
                rms.append(np.sqrt(np.mean(ymean**2)))

#             raw_input('neter')
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
        plt.savefig('GP_%s_%s' % (i, fname))
        best_time.append(lab)
    return rms2, mean_rms, lab, best_time

def dumb_sampling(P, nsamp, mins, ndays, sample_type, fname):

    # Generate time series with a GP
    interval = (60./mins)*24  #

    xs = np.linspace(0, ndays, interval*ndays) # one point every few minutes
    yerr = np.ones_like(xs)*.01

    if sample_type == 'GP':
        # Compute GP prior sample
        theta = [8.6969, 1.725e-3, 1.654, P]
        samples = sample_prior(theta, xs, yerr, nsamp)
    elif sample_type == "sine":
        samples = np.zeros((nsamp, len(xs)))
        for i in range(nsamp):
            samples[i, :] = sine_synth(xs)
    elif sample_type == "GPmix":
        samples = np.zeros((nsamp, len(xs)))
        for i in range(nsamp):
            samples[i, :] = GP_mix(xs)

    l = interval  # sample every day
    best_time = []
    for i, s in enumerate(samples):

        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(xs, s, color=ocols.blue)
        plt.xlabel("$\mathrm{Time~(days)}$")
        plt.ylabel("$\mathrm{RV~(ms}^{-1}\mathrm{)}$")
        plt.subplot(2, 1, 2)

#         times = 48
        times = 96  # number of times tested
        rms2 = np.zeros((times, times))
        ms = np.array(range(times))
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
                ymean = np.mean(ss, axis=0)  # take the mean of the 3

                stds.append(np.std(ymean))
                rms.append(np.sqrt(np.mean(ymean**2)))

                  # plot
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

