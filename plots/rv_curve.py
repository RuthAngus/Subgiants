import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
reb, fbt = plot_params()
from colours import plot_colours
cols = plot_colours()
import scipy.signal as sps

DIR = "/Users/angusr/Python/Subgiants"

# load rv curve
t, rv = np.genfromtxt("%s/injections/hd1100_rvs.txt" % DIR).T
e = 2.
# rv += e*np.random.randn(len(rv))

# load freqs
freqs = np.genfromtxt("%s/synthesise/freqs/hd1100_freqs.txt" % DIR).T

plt.clf()
plt.subplot(2, 1, 2)
fs = np.linspace(min(freqs)-1, max(freqs)+1, 1000)
pgram = sps.lombscargle(t, rv, fs*2*np.pi)
plt.plot(fs, pgram, color=cols.blue)
plt.xlabel("$\mathrm{Frequency~(days}^{-1}\mathrm{)}$")
plt.ylabel("$\mathrm{Power}$")

plt.subplot(2, 1, 1)
plt.plot(t, rv, color=cols.blue, zorder=1)
rv_err = np.ones_like(rv)*e
plt.fill_between(t, rv+rv_err, rv-rv_err, color=cols.blue, alpha=.2,
                 edgecolor=cols.blue)
plt.xlim(0, 3)
plt.xlabel("$\mathrm{Time~(days)}$")
plt.ylabel("$\mathrm{RV~(ms}^{-1}\mathrm{)}$")
plt.subplots_adjust(hspace=.3)
plt.savefig("%s/paper/hd1100_rvs.pdf" % DIR)

# load rv curve
t, rv = np.genfromtxt("%s/injections/3424541_rvs.txt" % DIR).T
e = 2.

# load freqs
freqs = np.genfromtxt("%s/data/3424541_freqs.txt" % DIR).T[1]

plt.clf()
plt.subplot(2, 1, 2)
fs = np.linspace(min(freqs)-1, max(freqs)+1, 1000)
pgram = sps.lombscargle(t, rv, fs*2*np.pi)
plt.plot(fs, pgram, color=cols.blue)
plt.xlabel("$\mathrm{Frequency~(days}^{-1}\mathrm{)}$")
plt.ylabel("$\mathrm{Power}$")

plt.subplot(2, 1, 1)
plt.plot(t, rv, color=cols.blue, zorder=1)
rv_err = np.ones_like(rv)*e
plt.fill_between(t, rv+rv_err, rv-rv_err, color=cols.blue, alpha=.2,
                 edgecolor=cols.blue)
plt.xlim(0, 3)
plt.xlabel("$\mathrm{Time~(days)}$")
plt.ylabel("$\mathrm{RV~(ms}^{-1}\mathrm{)}$")
plt.subplots_adjust(hspace=.3)
plt.savefig("%s/paper/3424541_rvs.pdf" % DIR)
