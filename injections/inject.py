import numpy as np
import matplotlib.pyplot as plt
import orbit
from model import model

# This script injects planets into RV curves

# create N planets with random periods and masses
# return theta
def random_planet(N, kid):

    np.random.seed(1234)
    # draw period from log uniform between 0.5 and 2 days
    P = np.random.uniform(np.log(0.1), np.log(20), N)

    np.random.seed(1234)
    # draw mass from log uniform between 0.5 and 4 Earth masses
    m2 = np.random.uniform(np.log(0.5), np.log(4), N)

    theta = np.zeros((5, N))
    theta[0, :] = np.exp(P)
    theta[1, :] = np.exp(m2)
    return theta

# Add white noise
def white_noise(y, yerr):
    return y + yerr*np.random.rand(len(y))

# take the parameters and inject the planet
def inject(theta, kid):
    # load data
    x, y = np.genfromtxt("%s/%s_rvs.txt" % (DIR, kid)).T
    yerr = np.ones_like(y)*2.  # make up uncertainties
    M1, M1_err = np.genfromtxt("%s/params/%s_mass.txt"
                               % (DIR, kid)).T
    ecc = 0.
    rv = model(theta, x, yerr, M1, ecc)  # compute planet rvs
    noisy_rv = white_noise(y+rv, yerr[0])  # add planet to lc with noise
    return x, noisy_rv, yerr

# generate a bunch of planets and inject them all
def rv_gen(N, kid):

    theta = random_planet(N, kid)
    P, m2, T0, V0, omega = theta

    for i in range(N):
        print "Period = ", P[i], "M2 = ", m2[i]
        x, rv, rv_err = inject(theta[:, i], kid)

        plt.clf()
        plt.plot(x, rv, "k.")
        plt.savefig("%s/rv_curves/%s_%s_rvs" % (DIR, i, kid))
        np.savetxt("%s/rv_curves/%s_%s_rvs.txt" % (DIR, i, kid),
                   np.transpose((x, rv, rv_err)))
        np.savetxt("%s/params/%s_%s_params.txt" % (DIR, i, kid),
                theta[:, i])

if __name__ == "__main__":

    DIR = "/Users/angusr/Python/Subgiants/injections"

    kid = "HD185"
    N = 100
    rv_gen(N, kid)
