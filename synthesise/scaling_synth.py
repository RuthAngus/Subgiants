import numpy as np
import matplotlib.pyplot as plt
from sin_tests import show_sine
import scaling_relations as sc

def gaussian(x, a, m, s):
    return a*np.exp(-(x - m)**2/(2*s**2))

# dn, nm in uHz
def amps_and_freqs(kid, dn, nm):

    # generate individual frequencies
    nf = 6
    freqs = 24*3600*1e-6*np.arange(nm-nf*dn, nm+nf*dn, dn)

    # make sure you don't get nf+1 freqs
    freqs = freqs[:nf*2]

    # generate amps
    a = gaussian(freqs, .5, nm, 3*dn)  # FIXME amplitude is made up
    amps = np.zeros(2*nf*2+1)
    for i in range(nf*2):
        amps[i*2] = a[i]
        amps[i*2+1] = a[i]
    amps[-1] = 1.

    # amplitudes and frequencies
    np.savetxt("%s/Subgiants/synthesise/freqs/%s_freqs.txt" % (D, kid),
               np.transpose(freqs))
    np.savetxt("%s/Subgiants/synthesise/freqs/%s_amps.txt" % (D, kid),
               np.transpose(amps))

def fake_rvs(x, kid):

    # load freqs and amps
    freqs = np.genfromtxt("%s/Subgiants/synthesise/freqs/%s_freqs.txt"
                          % (D, kid)).T
    amps = np.genfromtxt("%s/Subgiants/synthesise/freqs/%s_amps.txt"
                          % (D, kid)).T

    print freqs
    raw_input('meter')
    # generate rv curve
    ys = show_sine(x, freqs*2*np.pi, amps)
    np.savetxt("%s/Subgiants/injections/%s_rvs.txt" % (D, kid),
               np.transpose((x, ys)))

if __name__ == "__main__":

    D = "/Users/angusr/Python/"

    # load examples x values
    x, y = np.genfromtxt("3424541_rvs.txt").T

#     # load delta nu and nu_max
#     data = np.genfromtxt("%s/Gyro/data/ApJtable_zeros.txt" % D,
#                          skip_header=30, skip_footer=1892-549).T
#     kids, nm, nm_err, dn, dn_err = data[:5, :]

    # load Luan's subgiants
    data = np.genfromtxt("%s/Subgiants/proposal/sample_luan.out" % D,
                         skip_header=1, dtype=str).T
    kids = data[0]
    T, T_err = data[1], data[2]
    M, M_err = data[9], data[10]
    R, R_err = data[13], data[14]

    for i, kid in enumerate(kids[:5]):
        dn = sc.delta_nu(float(M[i]), float(R[i]))
        nm = sc.nu_max(float(M[i]), float(R[i]), float(T[i]))
        amps_and_freqs(kid, dn, nm*1e3)
        fake_rvs(x, kid)
