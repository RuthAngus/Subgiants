import numpy as np
import matplotlib.pyplot as plt
import orbit
from model import model

# This script injects planets into RV curves

# create N planets with random periods and masses
# return theta
def random_planet(N, kid, p, m, dist="uniform"):

    np.random.seed(1234)
    pmin, pmax = p
    if dist == "log":
        # draw period from log uniform between pmin and pmax days
        P = np.exp(np.random.uniform(np.log(pmin), np.log(pmax), N))
    elif dist == "uniform":
        P = np.random.uniform(pmin, pmax, N)

    np.random.seed(5678)
    mmin, mmax = m
    if dist == "log":
        # draw mass from log uniform between mmin and mmax Earth masses
        m2 = np.exp(np.random.uniform(np.log(mmin), np.log(mmax), N))
    elif dist == "uniform":
        m2 = np.random.uniform(mmin, mmax, N)

    theta = np.zeros((5, N))
    theta[0, :] = P
    theta[1, :] = m2
    return theta

# Add white noise
def white_noise(y, yerr):
    return y + yerr*np.random.rand(len(y))

# take the parameters and inject the planet
def inject(theta, kid):
    # load data
    x, y = np.genfromtxt("%s/%s_rvs.txt" % (DIR, kid)).T
    if kid == "HD185":
        x, y = np.genfromtxt("%s/%s_rvs_100.txt" % (DIR, kid)).T
    yerr = np.ones_like(y)*2.  # make up uncertainties
    M1, M1_err = np.genfromtxt("%s/params/%s_mass.txt"
                               % (DIR, kid)).T
    ecc = 0.
    rv = model(theta, x, yerr, M1, ecc)  # compute planet rvs
    noisy_rv = white_noise(y+rv, yerr[0])  # add planet to lc with noise
    return x, noisy_rv, yerr

# generate a bunch of planets and inject them all
def rv_gen(N, kid, p, m, sub=1):

    theta = random_planet(N, kid, p, m)
    P, m2, T0, V0, omega = theta

    for i in range(N):
        print "Period = ", P[i], "M2 = ", m2[i]
        x, rv, rv_err = inject(theta[:, i], kid)

        # subsample
        x, rv, rv_err = x[::sub], rv[::sub], rv_err[::sub]

        plt.clf()
        plt.plot(x, rv, "k.")
        plt.savefig("%s/rv_curves/%s_%s_rvs" % (DIR, i, kid))

        np.savetxt("%s/params/%s_%s_params.txt" % (DIR, i, kid),
                theta[:, i])
#         if sub == 1:
#             np.savetxt("%s/rv_curves/%s_%s_rvs.txt" % (DIR, i, kid),
#                        np.transpose((x, rv, rv_err)))
#         else:
        np.savetxt("%s/rv_curves/%s_%s_rvs_%s.txt" % (DIR, i, kid, sub),
                   np.transpose((x, rv, rv_err)))

if __name__ == "__main__":

    DIR = "/Users/angusr/Python/Subgiants/injections"

    kid = "HD185"
    N = 1000
    sub = 1
    p = (3, 100)
    m = (1, 100)
    rv_gen(N, kid, p, m, sub=sub)
