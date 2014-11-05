import numpy as np
import matplotlib.pyplot as plt
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel, WhiteKernel
from colors import plot_colors
ocols = plot_colors()
from rc_params import plot_params
reb = plot_params()
from BGdata import BetaGem
BG = BetaGem()
from HD82074 import HD
HD = HD()

def sampling(t, xs, ys, rv_err):
    # try calculating different rms
    ppd = 1./(max(t) - min(t)) * 1000.  # points per day
    ppt = np.array([24, 24*2, 24*3, 24*4, 24*5, 24*10,
                   24*20, 24*30, 24*60])  # pph, pphh, pptm, ppm
    rms = []
    for pts in ppt[::-1]:
        l = int(ppd/pts)
        x, y, yerr = xs[::l], ys[::l], np.ones_like(ys[::l])*np.mean(rv_err)
        print np.std(y)
        rms.append(np.sqrt(sum((y)**2)/len(y)))
    return ppt, rms

def prior_sample(t, rv, rv_err, theta, i):

    # sample from the prior
    k = theta[0] * ExpSquaredKernel(theta[1]) * \
            ExpSine2Kernel(theta[2], theta[3])
    gp = george.GP(k)
    xs = np.linspace(min(t), max(t), 1000)
    yerr_inds = np.random.randint(0, len(rv_err), len(xs))
    yerrs = rv_err[yerr_inds]
    gp.compute(xs, yerrs)
    np.random.seed(12)
    ys = gp.sample()

    ppt, rms = sampling(t, xs, ys, rv_err)

    plt.clf()
    plt.subplot(3, 1, 1)
    plt.errorbar(t, rv, yerr=rv_err, **reb)
    plt.subplot(3, 1, 2)
    plt.plot(xs, ys, ocols.maroon)
    plt.errorbar(x, y, yerr=yerr, **reb)
    plt.subplot(3, 1, 3)
    plt.plot(ppt/24, rms[::-1], color=ocols.orange)
    plt.xlabel("$\mathrm{Samples~per~hour}$")
    plt.ylabel("$\mathrm{RMS}~(ms^{-1})$")
    plt.savefig("sampling_strategy%s" % i)

if __name__ == "__main__":

    DIR = "/Users/angusr/Python/Subgiants"

    # load trained GP parameters
    BGtheta = np.exp(np.genfromtxt('%s/data/BGGP_results.txt' % DIR).T)
    i = 8
    HDtheta = np.exp(np.genfromtxt('%s/data/HDGP_results%s.txt' % (DIR, i)).T)

#     prior_sample(BG.rvHJD, BG.rv, BG.rv_err, BGtheta)
#     prior_sample(HD.t0, HD.rv0, HD.rv_err0, HDtheta, i)
#     prior_sample(HD.t1, HD.rv1, HD.rv_err1, HDtheta, i)
#     prior_sample(HD.t2, HD.rv2, HD.rv_err2, HDtheta, i)
#     prior_sample(HD.t3, HD.rv3, HD.rv_err3, HDtheta, i)
#     prior_sample(HD.t4, HD.rv4, HD.rv_err4, HDtheta, i)
#     prior_sample(HD.t5, HD.rv5, HD.rv_err5, HDtheta, i)
#     prior_sample(HD.t6, HD.rv6, HD.rv_err6, HDtheta, i)
#     prior_sample(HD.t7, HD.rv7, HD.rv_err7, HDtheta, i)
#     prior_sample(HD.t8, HD.rv8, HD.rv_err8, HDtheta, i)
