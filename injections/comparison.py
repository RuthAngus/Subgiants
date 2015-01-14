import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
reb, fbt = plot_params()
from colours import plot_colours
ocols = plot_colours()

# This script compares the injected and recovered planets

def compare(N, fname, MCMC=False):

    true_p, true_m = [], []
    success_p, success_m = [], []
    success_psub, success_msub = [], []
    p = np.zeros((3, N))
    m = np.zeros((3, N))
    psub = np.zeros((3, N))
    msub = np.zeros((3, N))
    for n in range(N):

        # load truths
        theta_true = np.genfromtxt("%s/params/%s_%s_params.txt" %
                                   (DIR, n, fname)).T

        if MCMC == True:
            # load results
            results = np.genfromtxt("%s/results/%s_%s_results.txt" %
                                       (DIR, n, fname)).T
            p[:, n] = np.exp(results[:, 0])
            m[:, n] = results[:, 1]
        else:
            # load results
            results = np.genfromtxt("%s/results/%s_%s_1_pgram.txt" %
                                       (DIR, n, fname)).T
            p[0, n] = results

        true_p.append(theta_true[0])
        true_m.append(theta_true[1])
        print true_p[n], p[0, n]

        if MCMC == True:
            # check that the result falls within 1 sigma of the truth
            if p[0, n]-p[1, n] < true_p[n] and true_p[n] < p[0, n]+p[2, n] and \
                    m[0, n]-m[1, n] < true_m[n] and true_m[n] < m[0, n]+m[2, n]:
                success_p.append(true_p[n])
                success_m.append(true_m[n])
        else:
            # check that the period lies within 10% of the truth
            f = .2
            if p[0, n] < true_p[n]+true_p[n]*f and \
                    true_p[n]-true_p[n]*f < p[0, n]:
                success_p.append(true_p[n])
                success_m.append(true_m[n])

    plt.clf()
    plt.plot(true_p, true_m, '.', color=ocols.pink, markersize=20, alpha=.2)
    plt.plot(success_p, success_m, '.', color=ocols.blue, markersize=15, alpha=.3)

    cols = (ocols.green, ocols.red)
    alphas = [.3, .5]
    ms = [10, 5]
    subs = [10, 20]
    for i in range(len(subs)):
        for n in range(N):
            # load sub results
            if MCMC == True:
                results10 = np.genfromtxt("%s/results/%s_%s_%s_results.txt" %
                                          (DIR, n, fname, subs[i])).T
                psub[:, n] = np.exp(results10[:, 0])
                msub[:, n] = results10[:, 1]

            else:
                results10 = np.genfromtxt("%s/results/%s_%s_%s_pgram.txt" %
                                          (DIR, n, fname, subs[i])).T
                psub[0, n] = results10

            if MCMC == True:
                # check that the result falls within 1 sigma of the truth
                if psub[0, n]-psub[1, n] < true_p[n] and \
                        true_p[n] < psub[0, n]+psub[2, n] and \
                        msub[0, n]-msub[1, n] < true_m[n] and \
                        true_m[n] < msub[0, n]+msub[2, n]:
                    success_psub.append(true_p[n])
                    success_msub.append(true_m[n])

            else:
                # check that the period lies within 10% of the truth
                if psub[0, n] < true_p[n]+true_p[n]*f and \
                        true_p[n]-true_p[n]*f < psub[0, n]:
                    success_p.append(true_p[n])
                    success_m.append(true_m[n])

        plt.plot(success_psub, success_msub, '.', color=ocols.blue,
                 alpha=alphas[i], markersize=ms[i])

    plt.xlabel("$\mathrm{P_{orb}~(Days)}$")
    plt.ylabel("$\mathrm{M_{planet}~(M_{earth})}$")
    if MCMC == True:
        plt.savefig("compare_subs_MCMC")
    else:
        plt.savefig("compare_subs_pgram")

if __name__ == "__main__":

    N = 100
    fname = "HD185"
    DIR = "/Users/angusr/Python/Subgiants/injections"
    compare(N, fname, MCMC=False)
