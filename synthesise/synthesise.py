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

# fit frequencies. for now use all of them?
def fit_freqs(t, rv, rv_err, freqs, reg):
    ts = t*24*3600  # convert to seconds
    ys, A = fit_sine_reg(ts, rv, rv_err, freqs*2*np.pi, reg)
    return ts, ys, A

if __name__ == "__main__":

    DIR = "/Users/angusr/Python/Subgiants"

    # load kids, masses and temperatures
    kid, teff, t_err, m, m_err = \
            np.genfromtxt("%s/data/AMP_subgiants.txt" % DIR, skip_header=1).T

    for i in range(len(kid)):
        print kid[i]

        print "load fits file"
        fname = str(int(kid[i]))
        x, y, yerr = load_data(fname)

#         print "truncate"
#         l = 10000
#         x, y, yerr = x[:l], y[:l], yerr[:l]

        print "convert to rv"
        rv, rv_err = convert_to_rv(y, yerr, teff[i], t_err[i])

        # convert ys to float64
        rv64 = np.empty((len(rv)))
        for i in range(len(rv)):
            rv64[i] = rv[i].astype('float64')

#         # compute periodogram
#         f = np.logspace(.01, 2, num=1000)
#         pgram = lombscargle(x, rv64, f*2*np.pi)
#         # load freqs
#         freqs = load_freqs(fname)
#         plt.clf()
#         plt.subplot(2, 1, 1)
#         for freq in freqs:
#             plt.axvline(freq*24*3600, color="r")
#         plt.plot(f, pgram)
#         plt.xlabel("freq")

        # high pass filter
        freqs = load_freqs(fname)
        # min(freqs) = 0.0005, max(freqs) = 0.001
        f_crit = 0.0005
        filtered_rv = highpass_filter(rv64, f_crit)

#         pgram = lombscargle(x, filtered_rv, f*2*np.pi)
#         plt.subplot(2, 1, 2)
#         for freq in freqs:
#             plt.axvline(freq*24*3600, color="r")
#         plt.plot(f, pgram)
#         plt.savefig('filter_test')

        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(x, rv, "k.")
        plt.subplot(2, 1, 2)
        plt.plot(x, filtered_rv, "k.")
        plt.savefig("test")

        print "fit regularised sine waves"
        reg = 1000
        ts, ys, A = fit_freqs(x, filtered_rv, rv_err, freqs[:1]*2*np.pi, reg)
#         ts, ys, A = fit_freqs(x, rv, rv_err, freqs*2*np.pi, reg)

#         plt.clf()
#         plt.plot(x, filtered_rv, "k.")
#         plt.plot(x, ys)
# #         plt.xlim(x[0], x[500])
#         plt.savefig("test")

        fs = np.linspace(0.00000001, max(freqs), 1000)
        pgram = lombscargle(ts, ys, fs*2*np.pi)

        plt.clf()
        for i in freqs:
            plt.axvline(i, color="r")
        plt.plot(fs, pgram)
#         plt.ylim(0, 10000)
        plt.ylim(0, 1)
        plt.savefig("frequency_detection")
        assert 0
