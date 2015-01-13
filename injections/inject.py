import numpy as np
import matplotlib.pyplot as plt
from orbit import radvel

# Assume circular. m1 is star's mass in solar masses, m2 is planet's
# mass in earth masses
# def model(x, P, m1, m2):
#     G = 6.67e-11
#     M_sun = 1.9891e30
#     M_earth = 5.972e24
#     m1, m2 = m1*M_sun, m2*M_earth
#     K = (2*np.pi*G/P)**(1./3.) * m2 * np.sin(90)/m1**(2./3.)
#     return radvel(x, P, K)

# Simple one planet model:
def model(theta, t, rverr, wn=False):
    P, M1, M2, T0, V0, Ecc, omega = theta
    G = 6.67e-11
    M_sun = 1.9891e30
    M_earth = 5.972e24
    m1, m2 = m1*M_sun, m2*M_earth
    K = (2*np.pi*G/P)**(1./3.) * m2 * np.sin(90)/m1**(2./3.)
    rv = orbit.radvel(t, P1, K1, T01, V01, Ecc1, omega1)
    if wn:
        rv += np.random.randn(len(t)) * rverr
    return rv

# Add white noise
def white_noise(y, yerr):
    return y + yerr*np.random.rand(len(y))

def inject(theta, kid):
    # load data
    DIR = "/Users/angusr/Python/Subgiants"
    x, y = np.genfromtxt("%s/observing_strategy/%s_rvs.txt" % (DIR, kid)).T
    yerr = np.ones_like(y)*2.  # make up uncertainties

    rv = model(theta, x, yerr)  # compute planet rvs
    noisy_rv = white_noise(y+rv, yerr[0])  # add planet to lc with noise
    return x, noisy_rv, yerr

# create N planets with random periods and masses
def random_planet(N, kid):
    m1, m1_err = np.genfromtxt("%s/params/%s_mass.txt" % (DIR, kid)).T

    # draw period from log uniform between 0.5 and 100 days
    P = np.random.uniform(np.log(0.5), np.log(100), N)

    # draw mass from log uniform between 0.5 and 15 Earth masses
    m2 = np.random.uniform(np.log(0.5), np.log(15), N)
    return np.exp(P), np.ones(N)*m1, np.exp(m2), np.zeros(N), np.zeros(N), \
            np.zeros(N)

def rv_gen(N):

    DIR = "/Users/angusr/Python/Subgiants/injections"
    kid = "HD185"

    theta = random_planet(N, kid)

    for i in range(N):
        print "Period = ", P[i], "M1 = ", m1, "M2 = ", m2[i]

        x, rv, rv_err = inject(kid, P[i], m1, m2[i])

        plt.clf()
        plt.plot(x, rv, "k.")
        plt.savefig("%s/rv_curves/%s_%s_rvs" % (DIR, i, kid))

        np.savetxt("%s/rv_curves/%s_%s_rvs.txt" % (DIR, i, kid),
                   np.transpose((x, rv, rv_err)))
        np.savetxt("%s/params/%s_%s_params.txt" % (DIR, i, kid),
                   np.transpose((P, m2)))

if __name__ == "__main__":

    N = 2
    rv_gen(N)
