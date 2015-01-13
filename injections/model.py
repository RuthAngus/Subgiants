import numpy as np
import orbit

# Assume circular. m1 is star's mass in solar masses, m2 is planet's
# mass in earth masses
# Simple one planet model
def model(theta, t, rverr, M1, ecc, wn=False):
    P, M2, T0, V0, omega = theta
    P = np.exp(P)
    G = 6.67e-11
    M_sun = 1.9891e30
    M_earth = 5.972e24
    M1, M2 = M1*M_sun, M2*M_earth
    K = (2*np.pi*G/P)**(1./3.) * M2 * np.sin(90)/M1**(2./3.)
    rv = orbit.radvel(t, P, K, T0, V0, ecc, omega)
    if wn:
        rv += np.random.randn(len(t)) * rverr
    return rv

