import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# an array of ntests of 3 observing times, separated by nmins, over ndays
# for a given starting time. start = integer
def obs_times(nmins, ndays, nsamples, start):
    t = nmins/60./24  # time interval (days)
    times = np.zeros(nsamples*(ndays-2))  # construct empty array
    # starting point is 'start' intervals before the first day
    st = start

    # calculate observing times
    t1 = np.arange(st, ndays-st, 1)  # one obs per day
    t2 = np.arange(st-(t), ndays-st-(t), 1)  # 2/day, separated by t
    t3 = np.arange(st+(t), ndays-st+(t), 1)  # 3/day, separated by t
    # make sure there are only ndays-2 times so you don't get edge effects
    if len(t1) > ndays-2 or len(t2) > ndays-2 or len(t3) > ndays-2:
        t1 = t1[:ndays-2]
        t2 = t2[:ndays-2]
        t3 = t3[:ndays-2]
    times = np.sort(np.concatenate((t1, t2, t3)))  # sort the times
    return times.T

# interpolation function
def interp(x, y, times):
    tck = interpolate.splrep(x, y, s=0)
    ynew = interpolate.splev(times, tck, der=0)
    return ynew

# simulate nsamp rv curves with the method defined by stype
def simulate(fname):
    DIR = '/Users/angusr/Python/Subgiants'
    x, y = np.genfromtxt("%s/injections/%s_rvs.txt" % (DIR, fname)).T
    e = 2.
    yerr = np.ones_like(y)*e  # make up uncertainties
    # add noise
#     y += e*np.random.randn(len(y))
    return x, y, yerr

def model(pars, x, y, yerr, args):
    nmins, start = pars
    ndays, nsamples, fname = args

    # generate an array of the observing times
    # ts is an array of observing times
    ts = obs_times(nmins, ndays, nsamples, start)
    print ts

    # calculate y values at observation positions
    # number of observations, number of tests
    ys = np.ndarray(len(ts.T))

    # xs you have, ys, xs you want
    ys = interp(x, y, ts)

    # calculate rms
    e = 2.
    nightly_av = np.mean(y, axis=1)

        # add observation errors in quadrature
        dx = sum(np.sqrt(2*(e/nightly_av)**2))/ndays
        x = np.mean(nightly_av**2)
        rms[i] = np.sqrt(np.mean(nightly_av**2))
        rms_err[i] = .5*(dx/x)*rms[i]
    l = rms == min(rms)
    return rms[l][0], rms_err[l][0]

# The likelihood function:
def lnlike(pars, x, y, yerr, args):
    min_rms, min_rms_err = model(pars, x, y, yerr, args)
    chisq = np.sum(((min_rms) / min_rms_err)**2)
    return -chisq / 2.

# uniform priors
def lnprior(pars):
    nmins, start = pars
    if 0 < nmins < 100 and 0 < start < 1:
        return 0.0
    else: return -np.inf

if __name__ == "__main__":

    nmins = 1  # in minutes
    ndays = 10
    ntests = 100
    nsamples = 3
    nsim = 1
    np.random.seed(123)
    start = np.random.uniform(0, 1)
    fname = 3424541

    # load data
    x, y, yerr = simulate(fname)

    pars = [nmins, start]
    args = [ndays, ntests, nsamples, fname]
    min_rms, min_rms_err = model(pars, x, y, yerr, args)
