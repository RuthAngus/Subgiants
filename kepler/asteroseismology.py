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

    download(KID)  # check to see if data are downloaded

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
    nu_maxHz = nu_max(m, r, teff)*1e-6
    fs = np.arange(1300e-6, 1700e-6, 1e-7) # Hz
    fs = np.linspace(5, 45, 10000) # c/d
    ws = 2*np.pi*fs  # lombscargle uses angular frequencies
    return sps.lombscargle(x, y2, ws)

if __name__ == "__main__":

    KID = '9098294'
    x, y, yerr = load_data(KID)
#     pgram = ls(x, y, yerr)

    nm, dn = 2233, 108.8  # uHz
    fs = gen_freqs_nm_dn(nm, dn, 10)
    ws = fs*2*np.pi

    plt.clf()
    plt.subplot(2, 1, 1)
    plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
    plt.subplot(2, 1, 2)
    plt.xlabel('$\mu Hz$')
#     plt.plot(ws/(2*np.pi)*1e6, pgram)
    for f in fs:
        plt.axvline(f, color='r')
    plt.savefig('KIC%s' % KID)
    print 'KIC%s.png' % KID
    plt.show()

    ys, A = fit_sine_err(x, y, yerr, ws)
    print A
