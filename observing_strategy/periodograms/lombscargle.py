import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps
from BGdata import BetaGem
BG = BetaGem()
from read_idl_file import read_file
from rc_params import plot_params
reb, rfb = plot_params()
from colours import plot_colours
ocols = plot_colours()
from scaling_relations import delta_nu, nu_max
import scipy.optimize as so

# compute figure 1 periodogram
def lombscar(x, y, fs, fname):
    ws = 2*np.pi*fs  # lombscargle uses angular freqs
    pgram = sps.lombscargle(x, y, ws)
    return pgram

# compute figure 1 periodogram
def lombscar_fig1(x, y, fname):
    fs = Xfreqs1()
    ws = 2*np.pi*fs  # lombscargle uses angular freqs
    pgram = sps.lombscargle(x, y, ws)
    return pgram, fs

# compute figure 2 periodogram
def lombscar_fig2(x, y, fname):
    fs = Xfreqs2()
    ws = 2*np.pi*fs  # lombscargle uses angular freqs
    pgram = sps.lombscargle(x, y, ws)
    return pgram, fs

# compute same frequencies as in Dumusque2014 Fig.1
def Xfreqs1():
    # freqs(in Xavier's paper they range from 10e-6 to 10e-2 Hz
    fmin, fmax, fstep = .0001, 0.014, .00001
    freqs = np.arange(fmin, fmax, fstep) # Hz
    return freqs

# compute same frequencies as in Dumusque2014 fig.2
def Xfreqs2():
    # freqs(in Xavier's paper they range from 10e-6 to 10e-2 Hz
    log_fmin, log_fmax, fstep = -6, -2, .01
    freqs = 10**np.arange(log_fmin, log_fmax, fstep) # Hz
    return freqs

# convert ys to float64
def float64ize(y):
    y2 = np.empty((len(y)))
    for i in range(len(y)):
        y2[i] = y[i].astype('float64')
    return y2

def VPSD(pars, f):

    A0, A1, A2, B0, B1, B2, C0, C1, C2, Al, Gamma, f_0, c = np.exp(pars)

    # granulation
    P0 = A0 / (1 + (B0*f)**C0)  # granulation
    P1 = A1 / (1 + (B1*f)**C1)  # mesogranulation
    P2 = A2 / (1 + (B2*f)**C2)  # mesogranulation

    # lorentzian
    Pl = Al * Gamma**2/((f - f_0)**2 + Gamma**2)

    return P0, P1, P2, Pl, P0+P1+P2+Pl+c

def resid(pars, x, y):
    return sum((y - VPSD(pars, x)[4])**2)

def spectral_analysis(t_days, y, yerr, m, r, t, fname):

    A0, A1, A2 = 1, 1, 0.5
    B0, B1, B2 = 30*60, 17*3600, 14*3600
    C0, C1, C2 = 10., 10., 10.
    Al = 3.
    Gamma = 1e-3
    f_0 = 1e-3
    c = 0.1
    pars = [A0, A1, A2, B0, B1, B2, C0, C1, C2, Al, Gamma, f_0, c]

    x = t_days*24*3600
    x = float64ize(x)
    y = float64ize(y)

    nm = nu_max(m, r, t)/1e3
    dn = delta_nu(m, r)/1e6
    print nm, dn

    pgram1, fs1 = lombscar_fig1(x, y, fname)
    plt.clf()
    plt.plot(fs1, pgram1, color=ocols.blue)
    plt.xlabel("$\mathrm{Frequency~(Hz)}$")
    plt.ylabel("$\mathrm{Power}$")
    plt.savefig('%spgram_fig1' % fname)

    pgram2, fs2 = lombscar_fig2(x, y, fname)
    # fit model using minimize
    results = so.minimize(resid, np.log(pars),
                          args=(np.log10(fs2), np.log10(pgram2)),
                          method='L-BFGS-B')
    print results.x

    plt.clf()
    plt.plot(np.log10(fs2), np.log10(pgram2), color=ocols.blue)
    plt.xlabel("$\log_{10}\mathrm{Frequency~(Hz)}$")
    plt.ylabel("$\mathrm{Power}$")
    plt.axvline(np.log10(nm), color=ocols.orange)
    p = VPSD(pars, fs2)
    p1, p2, p3, pl, p = VPSD(pars, fs2)
#     plt.plot(np.log10(fs2), p1, color='.5')
#     plt.plot(np.log10(fs2), p2, color='.5')
#     plt.plot(np.log10(fs2), p3, color='.5')
#     plt.plot(np.log10(fs2), pl, color='.5', linestyle='--')
    plt.plot(np.log10(fs2), VPSD(pars, fs2)[4], color='.5')
    plt.plot(np.log10(fs2), VPSD(results.x, fs2)[4], color=ocols.pink)
    plt.savefig('%spgram_fig2' % fname)

    fs = np.linspace(1e-6, 600e-6, 10000)
    pgram = lombscar(x, y, fs, fname)
    plt.clf()
    plt.plot(fs, pgram, color=ocols.blue)
    plt.xlabel("$\mathrm{Frequency~(Hz)}$")
    plt.ylabel("$\mathrm{Power}$")
    ndn = 4
    plt.axvline(nm, color=ocols.orange)
    for i in range(1, ndn):
        plt.axvline(nm+(i*dn), color=ocols.orange, alpha=1-(i/5.))
        plt.axvline(nm-(i*dn), color=ocols.orange, alpha=1-(i/5.))
    plt.savefig('%spgram' % fname)

if __name__ == "__main__":
    DIR = "/Users/angusr/Python/Subgiants"

    # Beta Gem
    BGx, BGy, BGyerr = BG.rvHJD-BG.rvHJD[0], BG.rv, BG.rv_err
    spectral_analysis(BGx, BGy, BGyerr, BG.m, BG.r, BG.teff, 'BG')

#     # Ashley's star
#     hdt, hdrv, hdrv_err = \
#             read_file('%s/data/vsthd82074.dat' % DIR)
#     l = 146  # truncate time series
# #     l = len(hdt)
#     data = np.genfromtxt('%s/data/HD82074_params.txt' % DIR, skip_header=1).T
#     m, r, t = data[0], data[2], data[4]
#     spectral_analysis(hdt[:l], hdrv[:l], hdrv_err[:l], m, r, t, 'HD')
#
#    # TRES star
#     BJD, rv, rv_err = np.genfromtxt('%s/TRES/HD210702.ccfSum.txt' % DIR,
#                                     skip_header=1, usecols=(0, 1, 2)).T
#     data = np.genfromtxt('%s/data/HD21072_params.txt' % DIR, skip_header=1,
#                          usecols=(2)).T
#     m, r, t = data[17], data[19], data[0]
#     print m, r, t
#     spectral_analysis(BJD, rv, rv_err, m, r, t, 'TRES')
