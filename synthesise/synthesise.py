import numpy as np
import matplotlib.pyplot as plt
import pyfits
import sc_target
from rc_params import plot_params
reb, fbt = plot_params()

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

# load frequencies
def load_freqs(fname):
    l, freqs = np.genfromtxt("%s/data/%s_freqs.txt" % (DIR, fname),
                          skip_header=1).T
    return freqs

# fit frequencies
# def fit_freqs(t, rv, rv_err, freqs):
#     return x, y, yerr

if __name__ == "__main__":

    DIR = "/Users/angusr/Python/Subgiants"

    # load kids, masses and temperatures
    kid, teff, t_err, m, m_err = \
            np.genfromtxt("%s/data/AMP_subgiants.txt" % DIR, skip_header=1).T

    for i in range(len(kid)):

        # load fits file
        fname = str(int(kid[i]))
        x, y, yerr = load_data(fname)

        # convert to rv
        rv, rv_err = convert_to_rv(y, yerr, teff[i], t_err[i])

        freqs = load_freqs(fname)
        print freqs
