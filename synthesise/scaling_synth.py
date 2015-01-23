import numpy as np
import matplotlib.pyplot as plt
from sin_tests import show_sine
import scaling_relations as sc

def gaussian(x, a, m, s):
    return a*np.exp(-(x-m)**2/(2*s**2))

# dn, nm in uHz
def amps_and_freqs(kid, dn, nm):

    # try using existing amps
    # A = np.genfromtxt("%s/Subgiants/synthesise/5955122_amps.txt" % D).T
    # nf = (len(A)-1)/4

    # generate individual frequencies in days-1
    nf = 6
    freqs = 24*3600*1e-6*np.arange(nm-(nf*dn), nm+(nf*dn), dn)

    # generate amps
    a = gaussian(freqs, 1, 24*3600*1e-6*nm, .5*dn)  # FIXME amp is made up

    amps = np.zeros(2*len(freqs)+1)
    for i in range(len(freqs)):
        amps[i*2] = a[i]
        amps[i*2+1] = a[i]
    amps[-1] = 0.

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

    # generate rv curve
    ys = show_sine(x, freqs*2*np.pi, amps)
    np.savetxt("%s/Subgiants/injections/%s_rvs.txt" % (D, kid),
               np.transpose((x, ys)))

if __name__ == "__main__":

    D = "/Users/angusr/Python/"

    # load examples x values
    x, y = np.genfromtxt("3424541_rvs.txt").T
    xs = np.linspace(min(x), max(x), len(x))

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

    for i, kid in enumerate(kids):
        dn = sc.delta_nu(float(M[i]), float(R[i]))
        nm = sc.nu_max(float(M[i]), float(R[i]), float(T[i]))
        amps_and_freqs(kid, dn, nm*1e3)
        fake_rvs(xs, kid)
