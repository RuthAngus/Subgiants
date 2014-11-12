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

def sampling_method(P, fname):

    # Generate time series with a GP
    interval = 6*24  # 10 mins
    interval = 12*24  # 5 mins
    ndays = 3
    xs = np.linspace(0, ndays, interval*ndays) # one point every 10 minutes
    yerr = np.ones_like(xs)*.01

    # Compute GP prior sample
    theta = [8.6969e+2, 1.725e-3, 1.654, P]
#     theta = [8.6969e+2, 1.725e+5, 1.654, P]
    # theta = results
    # theta = np.zeros(len(results)+1)
    # theta[:-1] = results
    # theta[-1] = P
    k = theta[0] * ExpSquaredKernel(theta[1]) * \
            ExpSine2Kernel(theta[2], theta[3])
    gp = george.GP(k)
    gp.compute(xs, yerr)
    np.random.seed(1234)
    samples = gp.sample(xs, 10)

    l = interval  # sample every day
    best_time = []
    for i, s in enumerate(samples):

        ms = np.array(range(18)) # from 0 to 120 minute intervals
        stds = []
        for m in ms:
            xs1, xs2, xs3 = xs[::l], xs[m::l], xs[2*m::l]
            ss1, ss2, ss3 = s[::l], s[m::l], s[2*m::l]
            yerrs1, yerrs2, yerrs3 = yerr[::l], yerr[m::l], yerr[2*m::l]

            xss = np.vstack((xs1, xs2, xs3))
            ss = np.vstack((ss1, ss2, ss3))
            yerrs = np.vstack((yerrs1, yerrs2, yerrs3))
            xmean = np.mean(xss, axis=0)
            ymean = np.mean(ss, axis=0)

            stds.append(np.std(ymean))

#             plt.clf()
#             plt.plot(xs, s, color=ocols.blue)
#             plt.plot(xs1, ss1, 'k.')
#             plt.plot(xs2, ss2, 'k.')
#             plt.plot(xs3, ss3, 'k.')
#             plt.plot(xs2, np.mean(ss, axis=0), '.', color=ocols.pink)
#             plt.show()
#             raw_input('enter')

        stds = np.array(stds)
        ll = stds==min(stds)
        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(xs, s, color=ocols.blue)
        plt.xlabel("$\mathrm{Time~(days)}$")
        plt.ylabel("$\mathrm{RV~(ms}^{-1}\mathrm{)}$")
        plt.subplot(2, 1, 2)
        mins = 5
        lab = ms[ll][0]*mins
        plt.plot(ms*mins, stds, color=ocols.orange,
                  label="$%s$" % lab)
        plt.xlabel("$\mathrm{Time~between~samples~(minutes)}$")
        plt.ylabel("$\mathrm{RMS}$")
        plt.subplots_adjust(hspace=.3)
        plt.legend()
        plt.savefig('GP_%s_%s' % (i, fname))

        best_time.append(ms[ll][0])

if __name__ == "__main__":

#     masses = np.linspace(1., 3., 10)
#     teffs =
#     Ps = np.linspace(1e-6, .1, 5)
#     for i, P in enumerate(Ps):
#         sampling_method(P, i)
    P = 2./24.
    sampling_method(P, 'test')
