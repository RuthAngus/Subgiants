import numpy as np
import matplotlib.pyplot as plt
from params import plot_params
reb = plot_params()

def mean(x, y, yerr):
    return np.mean(x), np.mean(y), np.sqrt(sum(yerr**2))

def rms(y):
    return np.sqrt(np.mean(y**2))

# given a time series return a list of lists
def nights(t, rv, rverr):
    ts, rvs, rverrs = [], [], []
    diff = np.diff(t)
    x = np.where(diff > .5)[0] + 1
    if x[0] == 1:
        ts.append([t[0]])
        rvs.append([rv[0]])
        rverrs.append([rverr[0]])
    for i in range(len(x)-1):
        ts.append(t[x[i]:x[i+1]])
        rvs.append(rv[x[i]:x[i+1]])
        rverrs.append(rverr[x[i]:x[i+1]])
    return ts, rvs, rverrs

def cm(x):
    s = 3
    n = 2
    ncombinations = s**n
    if len(x) > 6:
        n = 3
        c = range(1, s+1) + range(1, s+1) + range(1, s+1)
        comb = np.zeros((ncombinations, n))
        comb[:, 1] = np.sort(np.array(c))
        c = range(s+1, 2*s+1) + range(s+1, 2*s+1) + range(s+1, 2*s+1)
        comb[:, 2] = np.array(c)
    if len(x) < 7:
        c = range(s) + range(s) + range(s)
        comb = np.zeros((ncombinations, n))
        comb[:, 0] = np.sort(np.array(c))
        c = range(s, 2*s) + range(s, 2*s) + range(s, 2*s)
        comb[:, 1] = np.array(c)
    return comb

if __name__ == "__main__":

    D = "/Users/angusr/Python/Subgiants/paper"

    t1, rv1, rverr1 = np.genfromtxt("HD95089_PFS2.vels").T
    t2, rv2, rverr2 = HD98219 = np.genfromtxt("HD98219_PFS2.vels").T

    t0 = 2457061
    t1 -= t0
    t2 -= t0

    # split into list of nights
    ts1, rvs1, rverrs1 = nights(t1, rv1, rverr1)
    ts2, rvs2, rverrs2 = nights(t2, rv2, rverr2)

    bt1, brv1, brverr1 = [t1[0]], [rv1[0]], [rverr1[0]]
    b = mean(t1[1:4], rv1[1:4], rverr1[1:4])
    bt1.append(b[0])
    brv1.append(b[1])
    brverr1.append(b[2])
    b = mean(t1[5:], rv1[5:], rverr1[5:])
    bt1.append(b[0])
    brv1.append(b[1])
    brverr1.append(b[2])

    bt1, brv1, brverr1 = np.array(bt1), np.array(brv1), np.array(brverr1)

    c1 = cm(t1)
    rmss1 = []
    for c in c1:
        l = np.array([int(i) for i in c])
        rmss1.append(rms(rv1[l]))

    plt.clf()
    plt.errorbar(t1, rv1, yerr=rverr1, **reb)
    plt.errorbar(bt1, brv1, yerr=brverr1, fmt="r.", capsize=0, ecolor=".7")
    plt.axhline(rms(brv1), color="r")
    for r in rmss1:
        plt.axhline(r, color=".7")
    plt.xlabel("$\mathrm{JD}-%s$" % t0)
    plt.ylabel("$\mathrm{RV~(ms}^{-1}\mathrm{)}$")
    plt.savefig("%s/HD95089.pdf" % D)

    bt2, brv2, brverr2 = [], [], []
    b = mean(t2[:3], rv2[:3], rverr2[:3])
    bt2.append(b[0])
    brv2.append(b[1])
    brverr2.append(b[2])
    b = mean(t2[4:], rv2[4:], rverr2[4:])
    bt2.append(b[0])
    brv2.append(b[1])
    brverr2.append(b[2])

    bt2, brv2, brverr2 = np.array(bt2), np.array(brv2), np.array(brverr2)

    c2 = cm(t2)
    rmss2 = []
    for c in c2:
        l = np.array([int(i) for i in c])
        rmss2.append(rms(rv2[l]))

    plt.clf()
    plt.errorbar(t2, rv2, yerr=rverr2, **reb)
    plt.errorbar(bt2, brv2, yerr=brverr2, fmt="r.", capsize=0, ecolor=".7")
    plt.axhline(rms(brv2), color="r")
    for r in rmss2:
        plt.axhline(r, color=".7")
    plt.xlabel("$\mathrm{JD}-%s$" % t0)
    plt.ylabel("$\mathrm{RV~(ms}^{-1}\mathrm{)}$")
    plt.savefig("%s/HD98219.pdf" % D)

    print rv1
    print rv2
