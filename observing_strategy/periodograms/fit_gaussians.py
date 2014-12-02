# fit gaussians to the periodogram in fourier space
import numpy as np
import matplotlib.pyplot as plt
from lombscargle import gen_pgram, float64ize
from scaling_relations import nu_max, delta_nu
from BGdata import BetaGem
BG = BetaGem()
from rc_params import plot_params
reb, rfb = plot_params()
from colours import plot_colours
ocols = plot_colours()
import scipy.optimize as so
import george
from george.kernels import ExpSquaredKernel, CosineKernel

def resids(pars, x, y):
    return np.sqrt(sum((y - model(pars, x))**2))

def model(pars, x):
    mu1, sig1, A1, mu2, sig2, A2 = pars
    G1 = A1 * np.exp(-.5*(x-mu1)**2/sig1**2)
    G2 = A2 * np.exp(-.5*(x-mu2)**2/sig2**2)
    return G1 + G2

def plot(pars, x, y, nm, dn):
    # calculate periodogram
    fs, pgram = gen_pgram(x, y, nm, dn, 'BG')
    # plot pgram
    plt.clf()
    plt.plot(fs, pgram, color=ocols.blue)
    plt.xlabel("$\mathrm{Frequency~(Hz)}$")
    plt.ylabel("$\mathrm{Power}$")
    # plot modes
    ndn = 4
    plt.axvline(nm, color=ocols.orange)
    for i in range(1, ndn):
        plt.axvline(nm+(i*dn), color=ocols.orange, alpha=1-(i/5.))
        plt.axvline(nm-(i*dn), color=ocols.orange, alpha=1-(i/5.))
    xs = np.linspace(min(fs), max(fs), 1000)
    ys = model(pars, xs)
    plt.plot(xs, ys, color=ocols.green)
    plt.xlim(0, 2e-4)
#     plt.show()
    plt.savefig('fit_gauss')

def fitting(pars, args):
    result = so.minimize(resids, pars, arg)
    result = so.fmin(resids, pars, args=(x, y))
    return result

if __name__ == "__main__":
    DIR = "/Users/angusr/Python/Subgiants"

    # Beta Gem
    BGx, BGy, BGyerr = BG.rvHJD-BG.rvHJD[0], BG.rv, BG.rv_err
    x = BGx*24*3600
    x = float64ize(x)
    y = float64ize(BGy)
    y -= np.median(BGy)

    nm = nu_max(BG.m, BG.r, BG.teff)/1e3  # calculate nu_max
    dn = delta_nu(BG.m, BG.r)/1e6  # calculate delta_nu

    # plot initial guess Gaussians
    # mu1, sig1, A1, mu2, sig2, A2
    theta = (.9e-4, 2e-5, 1e3, 1.5e-5, 2e-5, 1e3)
    plot(theta, x, y, nm, dn)

    k = theta[2] * ExpSquaredKernel(1./(2*np.pi*theta[1]**2)) \
            * CosineKernel(theta[0]) \
            + theta[5] * ExpSquaredKernel(1./(2*np.pi*theta[4]**2)) \
            * CosineKernel(theta[3])
    gp = george.GP(k)
    l = 38
    gp.compute(x[:l], BGyerr[:l])

    xs = np.linspace(min(x), max(x[:l]), 1000)
    mu, cov = gp.predict(y[:l], xs)

    plt.clf()
    plt.errorbar(x[:l], y[:l], yerr=BGyerr[:l], **reb)
    plt.plot(xs, mu, color=ocols.blue)
#     plt.xlim(0, .01e7)
    plt.show()
