import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
reb, fbt = plot_params()
from colours import plot_colours
cols = plot_colours()

# load times
DIR = "/Users/angusr/Python/Subgiants"
fnames, times = np.genfromtxt("%s/observing_strategy/best_times.txt" % DIR).T

# load masses
KID, T, T_err, m, m_err, r, r_err, dnu, nm = \
        np.genfromtxt("%s/data/AMP_subgiants.txt" % DIR, skip_header=1).T

print 1/(nm/1e6*60)
print times

plt.clf()
plt.errorbar(times, m, yerr=m_err, fmt="ko")
plt.ylabel("$\mathrm{Mass}~(M_{Earth})$")
plt.xlabel("$\mathrm{Time~(mins)}$")
plt.xlim(0, 15)
plt.savefig("mass_vs_time")

plt.clf()
plt.errorbar(times, T, yerr=T_err, fmt="ko")
plt.ylabel("$\mathrm{T}_{eff}~(K)$")
plt.xlabel("$\mathrm{Time~(mins)}$")
plt.xlim(0, 15)
plt.savefig("T_vs_time")

plt.clf()
plt.errorbar(times, r, yerr=r_err, fmt="ko")
plt.errorbar(25, 5.35, yerr=.2, fmt="ro")
p = np.polyfit(times, r, 1)
xs = np.linspace(0, 30, 100)
ys = np.polyval(p, xs)
plt.plot(xs, ys, color=cols.blue, label="$y=%.2f~x+%.2f$" % (p[0], p[1]))
plt.legend(loc="best")
plt.ylabel("$\mathrm{Radius~(R}_{Earth}\mathrm{)}$")
plt.xlabel("$\mathrm{Optimum~time~interval~(mins)}$")
plt.xlim(0, 30)
plt.savefig("r_vs_time")

plt.clf()
plt.plot(times, dnu, "ko")
plt.ylabel("$\mathrm{dnu}$")
plt.xlabel("$\mathrm{Time~(mins)}$")
plt.xlim(0, 15)
plt.savefig("dnu_vs_time")

plt.clf()
plt.plot(times, nm, "ko")
plt.ylabel("$\mathrm{nm}$")
plt.xlabel("$\mathrm{Time~(mins)}$")
plt.xlim(0, 15)
plt.savefig("nm_vs_time")
