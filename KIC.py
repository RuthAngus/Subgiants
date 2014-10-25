import numpy as np
import matplotlib.pyplot as plt
import pyfits
from scipy.signal.spectral import lombscargle
import glob
from scaling_relations import nu_max, delta_nu
from astero import astero
astr = astero()

from rc_params import plot_params
from colors import plot_colors
plot_params = plot_params()
ocols = plot_colors()

# KIC 5955122
# mass = 1.28 + 0.11 - 0.15
# radius = 2.14 + 0.07 - 0.11
# teff = 5930 +/- 55 or 5966 +/- 162
m, r, t = 1.28, 2.14, 5930
print 'KIC 5955122'
print 'nu_max = ', nu_max(m, r, t)*1e3, 'uHz'
print 'dnu = ', delta_nu(m, r), 'uHz'

print 'KIC 2010607'
m, r, t = 1.37, 2.42, 6361
print 'nu_max = ', nu_max(m, r, t)*1e3, 'uHz'
print 'dnu = ', delta_nu(m, r), 'uHz'

print 'KIC 1725815'
m, r, t = 1.43, 2.02, 6550
print 'nu_max = ', nu_max(m, r, t)*1e3, 'uHz'
print 'dnu = ', delta_nu(m, r), 'uHz'

print 'KIC 1435467'
m, r, t = 1.16, 1.63, 6433
print 'nu_max = ', nu_max(m, r, t)*1e3, 'uHz'
print 'dnu = ', delta_nu(m, r), 'uHz'

print 'KIC 2309595'
m, r, teff = 1.14, 2.40, 5238

i = 0
KID, m, r, teff = int(astr.iKID[i]), astr.im[i], \
        astr.ir[i], astr.iteff[i]

print 'KIC ', KID
print m, r, teff
print 'nu_max = ', nu_max(m, r, teff)*1e3, 'uHz'
print 'dnu = ', delta_nu(m, r), 'uHeffz'
raw_input('enter')

KID = '00%s' % KID

# load data, median normalise and join together
lc_files = glob.glob("/Users/angusr/.kplr/data/lightcurves/%s/*_slc.fits" % KID)
print len(lc_files), ' files'
for i, lc_file in enumerate(lc_files):
    hdulist = pyfits.open(lc_file)
    tbdata = hdulist[1].data
    t = tbdata["TIME"]
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

# the frequency array
nu_maxHz = nu_max(m, r, teff)*1e-6
fs = np.arange(nu_maxHz-200e-6, nu_maxHz+200e-6, 1e-7) # Hz
ws = 2*np.pi*fs  # lombscargle uses angular frequencies

# convert ys to float64
y2 = np.empty((len(y)))
for i in range(len(y)):
    y2[i] = y[i].astype('float64')

y = y2
pgram = lombscargle(x, y2, ws)

print np.var(y)
plt.clf()
plt.subplot(2, 1, 1)
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
plt.subplot(2, 1, 2)
plt.xlabel('$\mu Hz$')
plt.plot(ws/(2*np.pi)*1e6, pgram)
plt.savefig('KIC%s' % KID)
