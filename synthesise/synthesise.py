import numpy as np
import matplotlib.pyplot as plt
import pyfits
import sc_target
from rc_params import plot_params
reb, fbt = plot_params()
from regularised_sine_fitting import fit_sine_reg
from scipy.signal import butter, lfilter, lombscargle
from sin_tests import show_sine

# load light curve. return x, y, yerr
def load_data(fname):
    kDIR = "/Users/angusr/.kplr/data/lightcurves/%s" % fname.zfill(9)

    # try 2 different quarters (not all stars have all qs!)
    try:
        hdulist = pyfits.open('%s/kplr%s-2010296114515_slc.fits' %
                              (kDIR, fname.zfill(9)))
    except:
        "IOError:"
    try:
        hdulist = pyfits.open('%s/kplr%s-2009231120729_slc.fits' %
                              (kDIR, fname.zfill(9)))
    except:
        "IOError:"
    tbdata = hdulist[1].data
    t = np.ascontiguousarray(tbdata["TIME"], dtype=np.float64)
    flux = tbdata["PDCSAP_FLUX"]
    flux_err = tbdata["PDCSAP_FLUX_ERR"]
    q = tbdata["SAP_QUALITY"]
    n = np.isfinite(t)*np.isfinite(flux)*np.isfinite(flux_err)*(q==0)
    x, flux, flux_err = t[n], flux[n], flux_err[n]
    return x, flux, flux_err

# convert flux to rv
def convert_to_rv(flux, flux_err, teff, t_err):
    rv, rv_err, dlL = sc_target.flux_rv(flux, flux_err, teff, t_err)

    # remove weird errorbars
    l = dlL==0
    rv_err[l] = np.mean(rv_err[np.isfinite(rv_err)])
    l = rv_err > 10
    rv_err[l] = np.mean(rv_err[rv_err<10])
    return rv, rv_err

# high-pass filter
def highpass_filter(data, f_crit, order=5):
    b, a = butter(order, f_crit, btype='high', analog=False)
    y = lfilter(b, a, data)
    return y

# load frequencies
def load_freqs(fname):
    l, freqs = np.genfromtxt("%s/data/%s_freqs.txt" % (DIR, fname),
                          skip_header=1).T
    return freqs/1e6  # convert to Hz

if __name__ == "__main__":

    DIR = "/Users/angusr/Python/Subgiants"

    # load kids, masses and temperatures
    kid, teff, t_err, m, m_err = \
            np.genfromtxt("%s/data/AMP_subgiants.txt" % DIR, skip_header=1).T

    for i in range(len(kid)):
        print kid[i]

        # load fits file
        fname = str(int(kid[i]))
        x, y, yerr = load_data(fname)

        x -= x[0]

        # cut off at 10 days
        l = x < 10.
        x, y, yerr = x[l], y[l], yerr[l]
#         l = 1000
#         x, y, yerr = x[:l], y[:l], yerr[:l]

        # convert to rv
        rv, rv_err = convert_to_rv(y, yerr, teff[i], t_err[i])

        # load frequencies
        freqs = load_freqs(fname)  # Hz

        # high pass filter min(freqs) = 0.0005, max(freqs) = 0.001
        # f_crit = 0.0005
        # filtered_rv = highpass_filter(rv, f_crit)

        # fit regularised sine waves. for now fit all of them?
        reg = 10
        ys, A = fit_sine_reg(x*24*3600, rv, rv_err, freqs[:1]*2*np.pi, reg)

        plt.clf()
        plt.plot(x, ys)
        plt.savefig("%s" % str(int(kid[i])))

        # next, check the hd185 synthetic light curve x values
        # save synthetic light curve in the same way maybe it doesn't
        # matter because I'm interpolating anyway?
        # for now just do 10 days worth again
        np.savetxt("%s/synthesise/%s_rvs.txt" % (DIR, kid[i]),
                   np.transpose((x, ys)))

#         plt.clf()
#         plt.plot(x*24*3600, ys, "b")
#         plt.plot(x*24*3600, rv, "k.")
#         plt.savefig("%srv" % kid[i])
#         f = np.linspace(.00001, 0.005, 1000)
#         pgram = lombscargle(x*24*3600, ys, f*2*np.pi)
#         plt.clf()
#         plt.plot(f, pgram)
#         plt.show()
#         raw_input('enter')
