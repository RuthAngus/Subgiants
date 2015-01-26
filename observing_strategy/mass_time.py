import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
reb, fbt = plot_params()
from colours import plot_colours
cols = plot_colours()

DIR = "/Users/angusr/Python/Subgiants"

# load masses
T, T_err, m, m_err, r, r_err = \
        np.genfromtxt("%s/proposal/sample_luan.out" % DIR, skip_header=1,
                      usecols=(1, 2, 9, 10, 13, 14)).T
# # load AMP stars
# fnames, times = np.genfromtxt("%s/observing_strategy/AMP_best_times.txt"
#                               % DIR).T
# KID, T, T_err, m, m_err, r, r_err, dnu, nm = \
#         np.genfromtxt("%s/data/AMP_subgiants.txt" % DIR, skip_header=1).T

N = 5
# load luan's subgiants
fnames, times, min_rms = \
        np.genfromtxt("%s/observing_strategy/named_best_times_%s.txt"
                      % (DIR, N)).T
l = len(fnames)
T, T_err, m, m_err, r, r_err = T[:l], T_err[:l], m[:l], m_err[:l], r[:l], \
        r_err[:l]

plt.clf()
plt.errorbar(times, m, yerr=m_err, fmt="ko")
plt.ylabel("$\mathrm{Mass}~(M_{\odot})$")
plt.xlabel("$\mathrm{Time~(mins)}$")
plt.savefig("mass_vs_time")

plt.clf()
plt.errorbar(times, T, yerr=T_err, fmt="ko")
plt.ylabel("$\mathrm{T}_{eff}~(K)$")
plt.xlabel("$\mathrm{Time~(mins)}$")
plt.savefig("T_vs_time")

plt.clf()
fig, ax1 = plt.subplots()
ax1.errorbar(times, r, yerr=r_err, fmt="ko")
# plt.errorbar(25, 5.35, yerr=.2, fmt="ro")
l = times < 100
p = np.polyfit(times[l], r[l], 3)
xs = np.linspace(0, max(times), 100)
ys = np.polyval(p, xs)
# ax1.plot(xs, ys, color=cols.blue, label="$y=%.2f~x+%.2f$" % (p[0], p[1]))
ax1.set_ylabel("$\mathrm{Radius~(R}_{Earth}\mathrm{)}$")
ax1.set_xlabel("$\mathrm{Optimum~time~interval~(mins)}$")
# ax2 = ax1.twinx()
# ax2.axhline(np.mean(min_rms), color=cols.orange,
#             label="$\mathrm{Mean~minimum~RMS}=%s$" % np.mean(min_rms))
# xs = range(0, 200)
# ax2.fill_between(xs, np.ones(200)*min(min_rms),
#                  np.ones(200)*max(min_rms), color=cols.orange, alpha=.3)
# ax2.set_ylim(0, .2)
# plt.legend(loc="best")
plt.savefig("/Users/angusr/Python/Subgiants/paper/r_vs_time_%s.pdf" % N)

# plt.clf()
# plt.plot(times, dnu, "ko")
# plt.ylabel("$\mathrm{dnu}$")
# plt.xlabel("$\mathrm{Time~(mins)}$")
# plt.xlim(0, 15)
# plt.savefig("dnu_vs_time")

# plt.clf()
# plt.plot(times, nm, "ko")
# plt.ylabel("$\mathrm{nm}$")
# plt.xlabel("$\mathrm{Time~(mins)}$")
# plt.xlim(0, 15)
# plt.savefig("nm_vs_time")
