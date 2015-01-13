import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
reb, fbt = plot_params()
from colours import plot_colours
ocols = plot_colours()

N = 10
true_p, true_m = [], []
success = []
p = np.zeros((3, N))
m = np.zeros((3, N))
for n in range(N):

    fname = "HD185"
    DIR = "/Users/angusr/Python/Subgiants/injections"

    # load truths
    theta_true = np.genfromtxt("%s/params/%s_%s_params.txt" %
                               (DIR, n, fname)).T

    # load results
    results = np.genfromtxt("%s/results/%s_%s_results.txt" %
                               (DIR, n, fname)).T

    true_p.append(theta_true[0])
    true_m.append(theta_true[1])
    p[:, n] = results[:, 0]
    m[:, n] = results[:, 1]

    # check that the result falls within 1 sigma of the truth
    if p[0, n]-p[1, n] < true_p[n] and true_p[n] < p[0, n]+p[2, n] and \
            m[0, n]-m[1, n] < true_m[n] and true_m[n] < m[0, n]+m[2, n]:

        success.append(true_p[n])

plt.clf()
plt.errorbar(p[0, :], m[0, :], yerr=(p[1, :], p[2, :]),
             xerr=(m[1, :], m[2, :]), **reb)
plt.plot(true_p, true_m, '.', color=ocols.pink)
plt.xlabel("$\mathrm{P_{orb}~(Days)}$")
plt.ylabel("$\mathrm{M_{planet}~(M_earth)}$")
plt.savefig("compare")
