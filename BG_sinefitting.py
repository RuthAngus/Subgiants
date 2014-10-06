import numpy as np
from BGdata import BetaGem
import matplotlib.pyplot as plt
from scipy.signal import periodogram, lombscargle
from scipy.misc import derivative
from rc_params import plot_params
from colors import plot_colors
ocol = plot_colors()

def sinefit(freqs, amps, xs):
    sine = np.zeros_like(xs)
    for i in range(len(freqs)):
        sine += amps[i] * np.sin(2*np.pi*freqs[i]*xs)
    return sine

dfreq, dfreq_err, freq, amp, amp_err = \
        np.genfromtxt("/Users/angusr/Python/Subgiants/BG_freqs.txt",
                      skip_header=1).T

# load flux data
BG = BetaGem()
x = BG.fHJD
y, yerr = BG.flux-np.median(BG.flux), BG.flux_err

# load rv data
rvx, rv, rverr = BG.rvHJD-2450000, BG.rv, BG.rv_err
xs = np.linspace(min(rvx), max(rvx), 100)

plt.clf()
plt.errorbar(rvx, rv, yerr=rverr, fmt='k.', capsize=0, ecolor='.8')
plt.plot(xs, sinefit(dfreq, amp, xs))
plt.xlim(77.48, 77.7)
plt.show()
