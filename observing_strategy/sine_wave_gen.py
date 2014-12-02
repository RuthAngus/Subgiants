import numpy as np
import matplotlib.pyplot as plt
from astero_modelling import gen_freqs_nm_dn
from BGdata import BetaGem
BG = BetaGem()
from asteroseismology import fit_sines
from sin_tests import show_sine
from rc_params import plot_params
reb, fbt = plot_params()
from colours import plot_colours
ocols = plot_colours()

def sine_synth():

    # load BG data
    x, y, yerr, nm, dn = BG.rvHJD*24*3600, BG.rv, BG.rv_err, BG.nm, BG.dnu

    # train on BG data
    fs, A = fit_sines(x, y, yerr, nm, dn, nfreqs 'BG')

    # generate time series
    xs = np.linspace(min(x), max(x), 1000)
    ys = show_sine(xs, fs*2*np.pi, A)

    plt.clf()
    plt.plot(xs, ys, color=ocols.blue)
    plt.xlabel('$\mathrm{Time~(s)}$')
    plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
    plt.savefig('sine_wave_gen')

    return xs, ys

if __name__ == "__main__":
    sine_synth()
