# fit gaussians to the periodogram in fourier space
# to do next: pick a different gaussian and try optimising
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
    G1 = pars[2] * np.exp(-.5*(x-pars[0])**2/pars[1]**2)
#     G2 = pars[5] * np.exp(-.5*(x-pars[3])**2/pars[4]**2)
#     return G1 + G2
    return G1

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

def GP_mix(xs):

    xs *= 24*3600
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
    thetad = (.9e-4, 2e-5, 1e3, 1.5e-5, 2e-5, 1e3)
    thetas = (1./.9e-4, 2e-5, 1e3)
    theta_mix = (1./.9e-4, 2e-5, 1e3, 1e2, 1e6)
#     plot(thetas, x, y, nm, dn)

    k_double = thetad[2] * ExpSquaredKernel(1./(2*np.pi*thetad[1]**2)) \
            * CosineKernel(thetad[0]) \
            + thetad[5] * ExpSquaredKernel(1./(2*np.pi*thetad[4]**2)) \
            * CosineKernel(thetad[3])
    k_single = thetas[2] * ExpSquaredKernel(1./(2*np.pi*thetas[1]**2)) \
            * CosineKernel(thetas[0])
    k_mix = theta_mix[2] * ExpSquaredKernel(1./(2*np.pi*theta_mix[1]**2)) \
            * CosineKernel(theta_mix[0]) + theta_mix[3] \
            * ExpSquaredKernel(theta_mix[4])

    k, theta = k_single, thetas

    gp = george.GP(k)
    l = 38
    gp.compute(x[:l], BGyerr[:l])
#     gp.optimize(x[:l], y[:l], BGyerr[:l])

#     xs = np.linspace(min(x), max(x[:l]), 1000)
    mu, cov = gp.predict(y[:l], xs)
    np.random.seed(3)
    s = gp.sample(xs, size=3)

    plt.clf()
    plt.subplot(2, 1, 1)
    plt.errorbar(x[:l], y[:l], yerr=BGyerr[:l], **reb)
    plt.plot(xs, mu, color=ocols.blue)
    plt.subplot(2, 1, 2)
    plt.plot(xs, s[0], color=ocols.blue)
    plt.plot(xs, s[1], color=ocols.pink)
    plt.plot(xs, s[2], color=ocols.green)
    plt.savefig('GP_mixture')
#     plt.show()
    return s[0]
