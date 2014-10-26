import numpy as np
import matplotlib.pyplot as plt
import george
from sin_tests import fit_sine
from rotation import before_and_after
from BGdata import BetaGem
import scipy.signal as sps
from scipy.misc import derivative
from scaling_relations import nu_max, delta_nu
from rc_params import plot_params
from colors import plot_colors
ocol = plot_colors()
BG = BetaGem()

# load period04 fitted frequencies
# freq, amp, phase = \
#         np.genfromtxt('/Users/angusr/Python/Subgiants/data/period04_results.txt',
#                       skip_header=1).T

# load frequencies and amplitudes
freq_day, f_err, freq_mHz, amp, err = \
        np.genfromtxt('/Users/angusr/Python/Subgiants/data/BG_freqs.txt',
                      skip_header=1).T

BGm = 1.91
BGm_err = 0.09
BGr = 8.8
BGr_err = 0.1
BGteff = 4750
BGteff_err = 150
nm = nu_max(BGm, BGr, BGteff)/1e3
dnu = delta_nu(BGm, BGr)/1e6

# load rv data
l = 38
rvx, rv, rverr = BG.rvHJD[:l]-2450000, BG.rv[:l], BG.rv_err[:l]
ys = fit_sine(rvx, rv, freq_day)

plt.clf()
plt.errorbar(rvx, rv, yerr=rverr, fmt='k.', capsize=0, ecolor='.8')
# plt.plot(rvx, ys-7681250)
plt.plot(rvx, ys)
plt.xlabel('Time (s)')
plt.ylabel('RV')
plt.show()

# load flux data
# BG = BetaGem()
# x = (BG.fHJD - BG.fHJD[0]) *24.*60.*60. # time zeroed and in seconds
# y, yerr = BG.flux-np.median(BG.flux), BG.flux_err

# # compute lomb scargle
# fs = np.arange(100e-6, 2000e-6, 1e-7) # Hz
# ws = 2*np.pi*fs  # lombscargle uses angular frequencies
# pgram = sps.lombscargle(x, y, ws)
# plt.clf()
# plt.plot(ws/2./np.pi, pgram)
# plt.axvline(2*np.pi*nm, color='r')
# print ws[pgram==max(pgram)]/2./np.pi, 'Hz'
# peaks = sps.find_peaks_cwt(pgram, np.arange(0.00001, 0.0002, 10))
# plt.show()

# knowing the freqency of nu_max, fit one sine wave to the data
# print nu_max, dnu
# freqs = np.array([nu_max, nu_max+dnu, nu_max-dnu])
# ys = fit_sine(x, y, freqs)

# plt.clf()
# plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
# plt.plot(x, ys, color=ocol.blue)
# plt.xlabel('Time (s)')
# plt.ylabel('Flux')
# plt.savefig('results')
