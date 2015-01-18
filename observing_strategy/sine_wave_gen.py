import numpy as np
import pyfits
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
from sc_target import flux_rv, hd185_rv
from synthesise import load_data

# a series of tools for generating and loading synthetic rv curves

DIR = "/Users/angusr/Python/Subgiants"

def HD185_data(rv=True):

    # load HD data and convert fluxes to rvs
    DIR = '/Users/angusr/Python/Subgiants'
    x, flux = np.genfromtxt('%s/data/hd185351.q16sc.ts' % DIR).T
    flux_err = np.ones_like(flux)*6.255e-5
    teff = 5042
    t_err = 0
    med = np.median(flux)
    flux /= med
    y, y_err, dlL = flux_rv(flux, flux_err, teff, t_err)
    y_err = np.ones_like(y)*2.

    # load 6 top frequencies and fit sines FIXME: use more than 6?
    fs, f_err = np.genfromtxt("%s/data/top_HD185_freqs.txt" % DIR).T
    fs *= 1e-6
    f_err *= 1e-6

    plt.clf()
    plt.errorbar(x, y, yerr=y_err, alpha=.3, **reb)
    plt.savefig("HD185")

    if rv==True:  # return either rvs or flux
        return x, y, y_err, fs, f_err
    return x, flux, flux_err, fs, f_err

# this code uses the top frequencies of HD185 to generate sines
# it fits to the first ndays in the light curve
def HDsine_fit(ndays):

    x, y, y_err, fs, f_err = HD185_data()
    print fs
    raw_input('enter')

    l = x < x[0]+ndays
    ys, A = fit_sine_err(x[l], y[l], y_err[l], fs*2*np.pi)
    print "trained on = ", max(x[l])-min(x[l]), "days"

    np.savetxt("HD185_amps.txt", A)

    return ys, A

# this function uses the amplitudes found above to generate sine time series
def HDsine_synth(xs, ndays, train=False):

    xs2 = xs*24*3600  # convert to seconds

    # train on the HD185 data and save the amplitudes
    # (do this the first time)
    if train == True:
        ys, A = HDsine_fit(ndays)

    # just use amplitudes already calculated to generate RV curve
    elif train == False:
        #  load amplitudes and freqs
        fs, f_err = np.genfromtxt("%s/data/top_HD185_freqs.txt" % DIR).T
        fs *= 1e-6
        f_err *= 1e-6
        A = np.genfromtxt("HD185_amps.txt").T

        # generate sine waves
        ys = show_sine(xs2, fs*2*np.pi, A)

        med = max(ys)  # MASSIVE HACK, FIXME!
        ys /= med
        ys *= 5

        plt.clf()
        plt.plot(xs, ys, color=ocols.blue)
        plt.xlabel('$\mathrm{Time~(s)}$')
        plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
        plt.savefig('HDsine_wave_gen')

        # save the synthetic light curve
        np.savetxt("HD185_rvs.txt", np.transpose((xs, ys)))
    return ys

def HD_synth_load():
    xs, ys = np.genfromtxt("HD185_rvs.txt").T
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

    # save the synthetic light curve
    np.savetxt("HD185_lc.txt", np.transpose((xs, ys)))
    return ys

# these functions use Kepler data to produce sine waves
def kepler_data(kid, rv=True):

    x, flux, flux_err = load_data(kid)

    # load kids, masses and temperatures
    k, teff, t_err, m, m_err = \
            np.genfromtxt("%s/data/AMP_subgiants.txt" % DIR, skip_header=1).T
    l = k==kid
    teff, t_err = teff[l], t_err[l]
    print k[l], teff, t_err

    # calculate rv
#     med = np.median(flux)
#     flux /= med

    y, y_err, dlL = flux_rv(flux, flux_err, teff, t_err)
    y_err = np.ones_like(y)*2.

    # load top frequencies and fit sines
    fs, f_err = np.genfromtxt("%s/data/%s_freqs.txt" % (DIR, kid)).T
    fs *= 1e-6
    f_err *= 1e-6

    plt.clf()
    plt.errorbar(x, y, yerr=y_err, alpha=.3, **reb)
    plt.savefig("%s_test" % kid)

    if rv:  # return either rvs or flux
        return x, y, y_err, fs, f_err
    return x, flux, flux_err, fs, f_err

# this code uses the top frequencies of kepler stars to generate sines
# it fits to the first ndays in the light curve
def kepler_sine_fit(kid, ndays):

    x, y, y_err, fs, f_err = kepler_data(kid)

    l = x < x[0]+ndays
    ys, A = fit_sine_err(x[l], y[l], y_err[l], fs*2*np.pi)
    print "trained on = ", max(x[l])-min(x[l]), "days"

    np.savetxt("%s_amps_test.txt" % kid, A)

    med = max(ys)  # MASSIVE HACK, FIXME!
    ys /= med
    ys *= 5

    return ys, A

# this function uses the amplitudes found above to generate sine time series
def kepler_sine_synth(kid, xs, ndays, train=False, fit=False):

    xs2 = xs*24*3600  # convert to seconds

    # (do this the first time)
    if train == True:
        ys, A = kepler_sine_fit(kid, ndays)

    # fit sine waves to the data
    if fit == True:
        ys, A = kepler_sine_fit(ndays)

    # just use amplitudes already calculated to generate RV curve
    elif fit == False:
        #  load amplitudes and freqs
        fs, f_err = np.genfromtxt("%s/data/%s_freqs.txt" % (DIR, kid)).T
        fs *= 1e-6
        f_err *= 1e-6
        A = np.genfromtxt("%s_amps.txt" % kid).T

        # generate sine waves
        ys = show_sine(xs2, fs*2*np.pi, A)

    med = max(ys)  # MASSIVE HACK, FIXME!
    ys /= med
    ys *= 5

    plt.clf()
    plt.plot(xs, ys, color=ocols.blue)
    plt.xlabel('$\mathrm{Time~(s)}$')
    plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
    plt.savefig('%s_sine_wave_gen' % kid)

    # save the synthetic light curve
    np.savetxt("%s_rvs_test.txt" % kid, np.transpose((xs, ys)))
    return ys

def kepler_synth_load(kid):
    xs, ys = np.genfromtxt("%s_rvs.txt" % kid).T
    return ys

if __name__ == "__main__":

    kid = "3424541"
    xs = np.linspace(0, 10, 1000)
    ndays = 10
    kepler_sine_synth(kid, xs, ndays, train=True, fit=True)
