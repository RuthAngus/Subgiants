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

# compute periodogram
def lombscar(x, y, fs, fname):
    ws = 2*np.pi*fs  # lombscargle uses angular freqs
    pgram = sps.lombscargle(x, y, ws)
    plt.clf()
    plt.plot(fs, pgram, color=ocols.blue)
#     plt.plot(np.log10(fs), pgram, color=ocols.blue)
#     plt.xlabel("$\log_{10}\mathrm{Frequency~(Hz)}$")
    plt.xlabel("$\mathrm{Frequency~(Hz)}$")
    plt.ylabel("$\mathrm{Power}$")
    plt.savefig('%spgram' % fname)
    return pgram

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

if __name__ == "__main__":
    DIR = "/Users/angusr/Python/Subgiants"

    # reproducing figure 1
    freqs = Xfreqs1()

    # Beta Gem
    BGx, BGy, BGyerr = BG.rvHJD-BG.rvHJD[0], BG.rv, BG.rv_err
    BGx *= 24*3600  # convert to seconds
    pgram = lombscar(BGx, BGy, freqs, 'BG')

#     # Ashley's star
#     hdt, hdrv, hdrv_err = \
#             read_file('%s/data/vsthd82074.dat' % DIR)
#     # convert ys to float64
#     y = float64ize(hdrv)
#     pgram = lombscar(hdt, y, freqs, 'HD')

   # TRES star
    BJD, rv, rv_err = np.genfromtxt('%s/TRES/HD210702.ccfSum.txt' % DIR,
                                    skip_header=1, usecols=(0, 1, 2)).T
    BJD -= BJD[0]
    BJD *= 24*3600  # convert to seconds
    pgram = lombscar(BJD, rv, freqs, 'TRES')
