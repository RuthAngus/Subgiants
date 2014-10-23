import numpy as np
import matplotlib.pyplot as plt
import idlsave
from rc_params import plot_params
from colors import plot_colors
params = plot_params()
ocols = plot_colors()

DIR = '/Users/angusr/Python/Subgiants'
data = idlsave.read("%s/data/vsthd82074.dat" % DIR)

t0 = 16700
t = data.cf3['JD'] - t0
rv = data.cf3['MNVEL']
rv_err = data.cf3['ERRVEL']

plt.clf()
plt.errorbar(t*24, rv, yerr=rv_err, ecolor='.8', capsize=0, fmt='k.')
plt.xlabel('$\mathrm{JD-%s}$' % t0)
plt.ylabel('$\mathrm{RV}~(ms^{-1})$')
# plt.xlim(1.6*24, 1.85*24)
# plt.xlim(85, 95)
# plt.xlim(110, 120)
# plt.xlim(135, 145)
# plt.xlim(155, 170)
plt.xlim(205, 215)
# plt.xlim(170, 300)
plt.ylim(-30, 30)
plt.savefig("%s/figures/HD82074" % DIR)
plt.show()
