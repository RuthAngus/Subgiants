import numpy as np
import matplotlib.pyplot as plt
from params import plot_params
reb = plot_params()
import itertools
from gatspy.periodic import LombScargle

def mean(x, y, yerr):
    x, y, yerr = np.array(x), np.array(y), np.array(yerr)
    return np.mean(x), np.mean(y), np.sqrt(sum(yerr**2))/len(x)

def rms(y):
    y = np.array(y)
    return np.sqrt(np.mean(y**2))

# given a list of lists, compare each item to each other item
def comb(x, y, yerr):
    arrays = y
    return list(itertools.product(*arrays))

# given a time series return a list of lists
def nights(t, rv, rverr, name):
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
    if name == "HD98219":
        ts.append(t[:3])
        rvs.append(rv[:3])
        rverrs.append(rverr[:3])
    ts.append(t[-3:])
    rvs.append(rv[-3:])
    rverrs.append(rverr[-3:])
    return ts, rvs, rverrs

def make_plot(t, rv, rverr, D, fname):

    t0 = 2457061
    t -= t0

    # split into list of nights
    ts, rvs, rverrs = nights(t, rv, rverr, fname)

    c = comb(ts, rvs, rverrs)

    # bin nightly observations
    bt = [mean(ts[i], rvs[i], rverrs[i])[0] for i in range(len(ts))]
    brv = [mean(ts[i], rvs[i], rverrs[i])[1] for i in range(len(ts))]
    brverr = [mean(ts[i], rvs[i], rverrs[i])[2] for i in range(len(ts))]

    # make plot
    plt.clf()
    plt.errorbar(t, rv, yerr=rverr, **reb)
    plt.errorbar(bt, brv, yerr=brverr, fmt="r.", capsize=0, ecolor=".7")
    rms_combinations = []
    for i in c:
#         plt.axhline(rms(i), color="k", alpha=.01)
#         plt.axhline(rms(i)/np.sqrt(3), color="b", alpha=.01)
        rms_combinations.append(rms(i))
#     plt.axhline(rms(brv), color="r")

    plt.xlabel("$\mathrm{JD}-%s$" % t0)
    plt.ylabel("$\mathrm{RV~(ms}^{-1}\mathrm{)}$")
    print "%s/%s.png" % (D, fname)
    plt.savefig("%s/%s" % (D, fname))
    return rms(brv), rms_combinations

def make_histogram(true_rms, rms_combinations, fname):
    plt.clf()
    plt.hist(rms_combinations, histtype="stepfilled", color="w")
    plt.hist(rms_combinations/np.sqrt(3), histtype="stepfilled", color="w",
             linestyle="dashed", fc=(0, 0, 0, 0))
    plt.axvline(true_rms, color="r")
    plt.xlabel("$\mathrm{RMS~(ms}^{-1}\mathrm{)}$")
    plt.ylabel("$\mathrm{Number~of~combinations}$")
    plt.savefig("results_hist_%s" % fname)
    print true_rms, np.median(rms_combinations), np.mean(rms_combinations)
    print np.median(rms_combinations)/true_rms

def pgram(t, y, dy, fname):
#     t *= 24*3600  # convert to seconds
#     fs = np.arange(50, 1000, .1)*1e-6  # Hz
#     period = 1. / fs
#     plt.plot(fs, power)
    model = LombScargle().fit(t, y, dy)
    period = np.linspace(.25/24, 5./24, 1000)
    fs = 1./period
    power = model.periodogram(period)
    plt.clf()
    plt.plot(period*24., power)
    print "%s_pgram" % fname
    plt.savefig("%s_pgram" % fname)
    print period[power==max(power)]*24

if __name__ == "__main__":

    D = "/Users/angusr/Python/Subgiants/paper"
    fnames = ["HD95089", "HD98219"]
    fnames = ["HD98219"]
    for fname in fnames:
        t, rv, rverr = np.genfromtxt("%s_PFS3.vels" % fname).T
        # fit a straight line to the data and subtract off
        Y = (np.matrix(rv)).T
        A = np.matrix((np.ones(len(t)), t)).T
        C = rverr**2 * np.identity(len(t))
        line = np.array((np.linalg.inv(A.T * np.linalg.inv(C) * A))
                * (A.T * np.linalg.inv(C) * Y))
        rv -= (line[0] + line[1]*t)
        true_rms, rms_combinations = make_plot(t, rv, rverr, D, fname)
        make_histogram(true_rms, rms_combinations, fname)
