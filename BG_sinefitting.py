import numpy as np
from BGdata import BetaGem
import matplotlib.pyplot as plt
from scipy.signal import periodogram, lombscargle
from scipy.misc import derivative
from scaling_relations import nu_max, delta_nu
from rc_params import plot_params
from colors import plot_colors
from sine_model import model_freq, lnlike_freq, lnprior_freq
from sine_model import lnprob_freq, MCMC
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
freq, amp, phase = \
        np.genfromtxt('/Users/angusr/Python/Subgiants/data/period04_results.txt',
                      skip_header=1).T

BGm = 1.91
BGm_err = 0.09
BGr = 8.8
BGr_err = 0.1
BGteff = 4750
BGteff_err = 150
nu_max = nu_max(BGm, BGr, BGteff)/1e3
dnu = delta_nu(BGm, BGr)

# load flux data
BG = BetaGem()
x = (BG.fHJD - BG.fHJD[0]) *24.*60.*60. # time zeroed and in seconds
y, yerr = BG.flux-np.median(BG.flux), BG.flux_err

# knowing the freqency of nu_max, fit one sine wave to the data
print nu_max
xs = np.linspace(min(x), max(x), 1000)

par_init = np.array([1., 0.])
args = (x, y, yerr, nu_max)
results = MCMC(par_init, args, lnlike_freq, lnprob_freq, lnprior_freq)
print results

plt.clf()
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
plt.plot(xs, model_freq(par_init, xs, nu_max), 'b')
plt.plot(xs, model_freq(results, xs, nu_max), 'r')
plt.show()

# # load rv data
# rvx, rv, rverr = BG.rvHJD-2450000, BG.rv, BG.rv_err
# xs = np.linspace(min(rvx), max(rvx), 100000)

# plt.clf()
# plt.errorbar(rvx, rv, yerr=rverr, fmt='k.', capsize=0, ecolor='.8')
# plt.plot(xs, sinefit(xs, freq, amp, phase))
# plt.xlim(77.48, 77.7)
# plt.show()

# plt.clf()
# plt.subplot(3, 1, 1)
# plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
# plt.subplot(3, 1, 2)
# t = x
# sp = np.fft.fft(y)
# # freq = np.fft.fftfreq(t.shape[-1])
# # freq = np.linspace(1e-6, 1e-5, t.shape[-1])
# freq = np.linspace(1e-6, 1e-5, t.shape[-1])
# plt.plot(freq, sp.real)
# # plt.xlim(0, .6)
# plt.xlim(0, 1e-5)
# plt.subplot(3, 1, 3)
# freq = np.linspace(0, 1, len(x))
# pgram = lombscargle(x, y, freq)
# plt.plot(freq, pgram)
# plt.show()
#
