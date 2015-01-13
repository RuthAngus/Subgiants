import numpy as np
import matplotlib.pyplot as plt
import orbit
from model import model
from recovery import MCMC

# create N planets with random periods and masses
# return theta
def random_planet(N, kid):
    # draw period from log uniform between 0.5 and 100 days
    P = np.random.uniform(np.log(0.5), np.log(100), N)

    # draw mass from log uniform between 0.5 and 15 Earth masses
    m2 = np.random.uniform(np.log(0.5), np.log(15), N)

    theta = np.zeros((6, N))
    theta[0, :] = np.exp(P)
    theta[1, :] = np.exp(m2)
    return theta

# Add white noise
def white_noise(y, yerr):
    return y + yerr*np.random.rand(len(y))

# take the parameters and inject the planet
def inject(theta, kid):
    # load data
    DIR = "/Users/angusr/Python/Subgiants"
    x, y = np.genfromtxt("%s/observing_strategy/%s_rvs.txt" % (DIR, kid)).T
    yerr = np.ones_like(y)*2.  # make up uncertainties
    M1, M1_err = np.genfromtxt("%s/injections/params/%s_mass.txt"
                               % (DIR, kid)).T
    rv = model(theta, x, yerr, M1)  # compute planet rvs
    noisy_rv = white_noise(y+rv, yerr[0])  # add planet to lc with noise
    return x, noisy_rv, yerr

# generate a bunch of planets and inject them all
def rv_gen(N, kid):

    theta = random_planet(N, kid)
    P, m2, T0, V0, ecc, omega = theta

    for i in range(N):
        print "Period = ", P[i], "M2 = ", m2[i]
        x, rv, rv_err = inject(theta[:, i], kid)

        plt.clf()
        plt.plot(x, rv, "k.")
        plt.savefig("%s/rv_curves/%s_%s_rvs" % (DIR, i, kid))

        np.savetxt("%s/rv_curves/%s_%s_rvs.txt" % (DIR, i, kid),
                   np.transpose((x, rv, rv_err)))
        np.savetxt("%s/params/%s_%s_params.txt" % (DIR, i, kid),
                   np.transpose((P, m2)))

if __name__ == "__main__":

    DIR = "/Users/angusr/Python/Subgiants/injections"

    kid = "HD185"
    N = 2
    rv_gen(N)
