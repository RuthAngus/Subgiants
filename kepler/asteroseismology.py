# Attempting to perform asteroseismology on kepler data
import numpy as np
import matplotlib.pyplot as plt
import pyfits
import glob
import scipy.signal as sps
from astero_modelling import gen_freqs_nm_dn
import kplr
client = kplr.API()
from astero_modelling import model
from sin_tests import fit_sine_err
from rc_params import plot_params
reb, fbt = plot_params()
from colours import plot_colours
ocols = plot_colours()

def download(KID):
    # load data, median normalise and join together
    lc_files = glob.glob("/Users/angusr/.kplr/data/lightcurves/%s/*_slc.fits"
                         % str(KID).zfill(9))
    # check to see whether you've already downloaded data...
    print len(lc_files), ' files'
    if len(lc_files) == 0:
        print 'downloading..'
        star = client.star("%s" % KID)
        lcs = star.get_light_curves(short_cadence=True, fetch=True)
        lc_files = glob.glob("/Users/angusr/.kplr/data/lightcurves/%s/*_slc.fits"
                             % str(KID).zfill(9))
    print lc_files

def load_data(KID):

#     download(KID)  # check to see if data are downloaded

    # load data, median normalise and join together
    lc_files = glob.glob("/Users/angusr/.kplr/data/lightcurves/%s/*_slc.fits"
                         % str(KID).zfill(9))
    for i, lc_file in enumerate(lc_files):
        hdulist = pyfits.open(lc_file)
        tbdata = hdulist[1].data
        t = np.ascontiguousarray(tbdata["TIME"], dtype=np.float64)
        flux = tbdata["PDCSAP_FLUX"]
        flux_err = tbdata["PDCSAP_FLUX_ERR"]
        q = tbdata["SAP_QUALITY"]
        n = np.isfinite(t)*np.isfinite(flux)*np.isfinite(flux_err)*(q==0)
        t, flux, flux_err = t[n], flux[n], flux_err[n]
        fluxmed = np.median(flux)
        flux /= fluxmed
        flux_err /= fluxmed
        if i == 0:
            x, y, yerr = t, flux, flux_err
            break
        else:
            x = np.concatenate((x, t))
            y = np.concatenate((y, flux))
            yerr = np.concatenate((yerr, flux_err))

    # convert t to seconds
    x *= 24.*60.*60.

    # convert ys to float64
    y2 = np.empty((len(y)))
    for i in range(len(y)):
        y2[i] = y[i].astype('float64')

    return x, y2, yerr

def ls(x, y, yerr):
    # the frequency array
    fs = np.linspace(1500e-6, 3000e-6, 1000) # Hz
    ws = 2*np.pi*fs  # lombscargle uses angular frequencies
    return ws, sps.lombscargle(x, y, ws)

if __name__ == "__main__":

    KID = '9098294'
    x, y, yerr = load_data(KID)
    l = 5000
    x, y, yerr = x[:l], y[:l], yerr[:l]

    nm, dn = 2233e-6, 108.8e-6  # uHz
    fs = gen_freqs_nm_dn(nm, dn, 5)

    # fit sine waves
    ws2 = fs*2*np.pi
    ys, A = fit_sine_err(x, y, yerr, ws2)

    # calculate amplitudes
    A_even = A[::2]
    A_odd = A[1::2]
    amps1 = np.sqrt((A_even[:-1] + A_odd)**2)
    amps4 = np.sqrt(A_even[:-1]**2) + np.sqrt(A_odd)**2
    amps2 = np.sqrt((A_even[:-1] / max(A_even[:-1]))**2)
    amps3 = np.sqrt((A_odd / max(A_odd))**2)
    amps1, amps2, amps3 = amps1/max(amps1), amps2/max(amps2), amps3/max(amps3)
    amps4 /= max(amps4)

    # compute and plot periodogram
    ws, pgram = ls(x, y, yerr)
    plt.clf()
    plt.subplot(2, 1, 1)
#     plt.plot(x, ys, color=ocols.blue)
    plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
    plt.subplot(2, 1, 2)
    plt.xlabel('$\mu Hz$')
    plt.plot(ws/(2*np.pi)*1e6, pgram, color=ocols.blue)
    for i, f in enumerate(fs):
        plt.axvline(f*1e6, ymax=amps1[i], color=ocols.orange)
#         plt.axvline(f*1e6, ymax=amps2[i], color=ocols.pink)
#         plt.axvline(f*1e6, ymax=amps3[i], color=ocols.orange)
#         plt.axvline(f*1e6, ymax=amps4[i], color=ocols.orange)
    plt.savefig('KIC%s' % KID)
    print 'KIC%s.png' % KID

