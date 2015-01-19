import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, a, m, s):
    return a*np.exp(-(x - m)**2/(2*s**2))

def amps_and_freqs(dn, nm):

    # load amplitudes
    amps = np.genfromtxt("%s/Subgiants/observing_strategy/HD185_amps.txt"
                         % D).T

    for k in range(len(kid)):
        # generate individual frequencies
        nf = 12
        freqs = np.arange(nm[k]-nf*dn[k], nm[k]+nf*dn[k], dn[k])
        # generate amps
        amps = gaussian(freqs, max(amps), nm[k], 3*dn[k])
        np.savetxt("%s_freqs_amps.txt", np.transpose((freqs, amps)))

if __name__ == "__main__":

    D = "/Users/angusr/Python/"
    data = np.genfromtxt("%s/Gyro/data/ApJtable_zeros.txt" % D,
                         skip_header=30, skip_footer=1892-549).T
    kid, nm, nm_err, dn, dn_err = data[:5, :]
