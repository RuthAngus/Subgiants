import numpy as np
import matplotlib.pyplot as plt
import pyfits
import glob
import kplr
from scipy.signal.spectral import lombscargle
from scaling_relations import nu_max, delta_nu
from astero import astero
from rc_params import plot_params
from colors import plot_colors

client = kplr.API()
plot_params = plot_params()
ocols = plot_colors()
astr = astero()

i = 6
KID, m, r, teff = int(astr.iKID[i]), astr.im[i], \
        astr.ir[i], astr.iteff[i]

print 'KIC ', KID
print m, r, teff
print 'nu_max = ', nu_max(m, r, teff)*1e3, 'uHz'
print 'dnu = ', delta_nu(m, r), 'uHz'

# load data, median normalise and join together
lc_files = glob.glob("/Users/angusr/.kplr/data/lightcurves/%s/*_slc.fits"
                     % str(KID).zfill(9))
# check to see whether you've already downloaded data...
print len(lc_files), ' files'
if len(lc_files) == 0:
    star = client.star("%s" % KID)
    lcs = star.get_light_curves(short_cadence=True, fetch=True)
    lc_files = glob.glob("/Users/angusr/.kplr/data/lightcurves/%s/*_slc.fits"
                         % str(KID).zfill(9))
print lc_files

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
print 'KIC%s.png' % KID
