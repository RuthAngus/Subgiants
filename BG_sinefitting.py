import numpy as np
import matplotlib.pyplot as plt
import george
from rotation import before_and_after
from BGdata import BetaGem
from scipy.signal import periodogram, lombscargle
from scipy.misc import derivative
from scaling_relations import nu_max, delta_nu
from rc_params import plot_params
from colors import plot_colors
from sine_model import model_freq, lnlike_freq
from sine_model import lnprob5_freq, lnprior5_freq
from sine_model import lnprob4_freq, lnprior4_freq
from sine_model import model3_freq, lnlike3_freq, lnprior3_freq
from sine_model import lnprob3_freq, MCMC
ocol = plot_colors()

# dfreq, dfreq_err, freq, amp, amp_err = \
#         np.genfromtxt("/Users/angusr/Python/Subgiants/BG_freqs.txt",
#                       skip_header=1).T

# load period04 fitted frequencies
# freq, amp, phase = \
#         np.genfromtxt('/Users/angusr/Python/Subgiants/data/period04_results.txt',
#                       skip_header=1).T

BGm = 1.91
BGm_err = 0.09
BGr = 8.8
BGr_err = 0.1
BGteff = 4750
BGteff_err = 150
nu_max = nu_max(BGm, BGr, BGteff)/1e3
dnu = delta_nu(BGm, BGr)/1e6

# load flux data
BG = BetaGem()
x = (BG.fHJD - BG.fHJD[0]) *24.*60.*60. # time zeroed and in seconds
y, yerr = BG.flux-np.median(BG.flux), BG.flux_err

# knowing the freqency of nu_max, fit one sine wave to the data
print nu_max, dnu

freqs = np.array([nu_max, nu_max+dnu, nu_max-dnu])
print results

plt.clf()
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
plt.plot(xs, model_freq(results, xs, freqs), color=ocol.blue)
plt.xlabel('Time (s)')
plt.ylabel('Flux')
plt.savefig('results')

# # load rv data
# rvx, rv, rverr = BG.rvHJD-2450000, BG.rv, BG.rv_err
# xs = np.linspace(min(rvx), max(rvx), 100000)
