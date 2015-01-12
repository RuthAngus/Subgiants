import numpy as np
import matplotlib.pyplot as plt
from astero_modelling import gen_freqs_nm_dn
from BGdata import BetaGem
BG = BetaGem()
from asteroseismology import fit_sines
from sin_tests import show_sine, fit_sine_err
from rc_params import plot_params
reb, fbt = plot_params()
from colours import plot_colours
ocols = plot_colours()

DIR = "/Users/angusr/Python/Subgiants"

# this code uses the top frequencies of HD185 to generate sines
# it fits to the first 1000 data points in the light curve
def HDsine_fit():

    # load HD data
    from sc_target import flux_rv, hd185_rv
    DIR = '/Users/angusr/Python/Subgiants'
    x, y = np.genfromtxt('%s/data/hd185351.q16sc.ts' % DIR).T
    y_err = np.ones_like(y)*6.255e-5

    # load frequencies and fit sines
    fs, f_err = np.genfromtxt("%s/data/top_HD185_freqs.txt" % DIR).T
    fs *= 1e-6
    f_err *= 1e-6

#     l = 5000
#     ys, a = fit_sine_err(x[:l], y[:l], y_err[:l], fs*2*np.pi)
#     print "trained on = ", x[l]-x[0], "days"

    l = x < x[0]+10
    ys, a = fit_sine_err(x[l], y[l], y_err[l], fs*2*np.pi)
    print "trained on = ", max(x[l])-min(x[l]), "days"

    np.savetxt("HD185_amps.txt", A)

    return ys, A

# this function uses the amplitudes found above to generate sine time series
def HDsine_synth(xs, train=False):

    xs2 = xs*24*3600  # convert to seconds

    if train == True:
        HDsine_fit()

    #  load amplitudes and freqs
    fs, f_err = np.genfromtxt("%s/data/top_HD185_freqs.txt" % DIR).T
    fs *= 1e-6
    f_err *= 1e-6
    A = np.genfromtxt("HD185_amps.txt").T

    # generate sine waves
    ys = show_sine(xs2, fs*2*np.pi, A)

#     plt.clf()
#     plt.plot(xs, ys, color=ocols.blue)
#     plt.xlabel('$\mathrm{Time~(s)}$')
#     plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
#     plt.savefig('HDsine_wave_gen')

    return ys

# This code uses BG data to produce sine waves
def BGsine_synth(xs):
    xs *= 24*3600

    # load BG data
    x, y, yerr, nm, dn = BG.rvHJD*24*3600, BG.rv, BG.rv_err, BG.nm, BG.dnu

    # train on BG data
    nfreqs = 5
    fs, A = fit_sines(x, y, yerr, nm, dn, nfreqs, 'BG')

    # generate time series
#     xs = np.linspace(min(x), max(x), 1000)
    ys = show_sine(xs, fs*2*np.pi, A)

    plt.clf()
    plt.plot(xs, ys, color=ocols.blue)
    plt.xlabel('$\mathrm{Time~(s)}$')
    plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
    plt.savefig('sine_wave_gen')

    return ys

if __name__ == "__main__":
    BGsine_synth()
