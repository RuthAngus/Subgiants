import numpy as np
from BGdata import BetaGem
import matplotlib.pyplot as plt
from scipy.signal import periodogram, lombscargle
from scipy.misc import derivative
from rc_params import plot_params
from colors import plot_colors
ocol = plot_colors()

def sinefit(xs, freq, amp, phase):
    signal = np.zeros(len(xs))
    for i in range(len(freq)):
        print amp[i], freq[i], phase[i]
        signal += amp[i] * np.sin(2*np.pi*freq[i]*xs - phase[i])
        print signal
    return signal

# dfreq, dfreq_err, freq, amp, amp_err = \
#         np.genfromtxt("/Users/angusr/Python/Subgiants/BG_freqs.txt",
#                       skip_header=1).T

# load period04 fitted frequencies
freq, amp, phase = np.genfromtxt('period04_results.txt',
                                 skip_header=1).T

# load flux data
BG = BetaGem()
x = BG.fHJD
y, yerr = BG.flux-np.median(BG.flux), BG.flux_err

plt.clf()
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
plt.show()

# load rv data
rvx, rv, rverr = BG.rvHJD-2450000, BG.rv, BG.rv_err
xs = np.linspace(min(rvx), max(rvx), 100000)

plt.clf()
plt.errorbar(rvx, rv, yerr=rverr, fmt='k.', capsize=0, ecolor='.8')
plt.plot(xs, sinefit(xs, freq, amp, phase))
plt.xlim(77.48, 77.7)
plt.show()
