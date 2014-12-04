import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps
import scaling_relations as sr
from colours import plot_colours
ocols = plot_colours()

def flux_rv(y, teff):
    dlL = (y - np.median(y)) * 1e6  # ppm
    return dlL/17.7 * teff/5777., dlL

DIR = '/Users/angusr/Python/Subgiants'
x, y = np.genfromtxt('%s/data/hd185351.q16sc.ts' % DIR).T

# convert to seconds
l = len(x)
x, y = x[:l]*24*3600, y[:l]

# physical parameters
teff = 5045
logg = 3.27
nm = sr.nu_max_alt(logg, teff)/1e3
dn = sr.delta_nu_alt(.5)/1e6  # total guess
print 'nu_max = ', nm*1e6, 'uHz'
print 'delta nu guess = ', dn*1e6, 'uHz'

# calculate periodogram
fs = np.linspace(100e-6, 2000e-6, 1000)
ws = 2*np.pi*fs
pgram = sps.lombscargle(x, y, ws)

# convert flux to rvs
rv, dlL = flux_rv(y, teff)

# plot
plt.clf()
plt.subplot(3, 1, 1)
plt.plot(x, y, 'k.')
plt.subplot(3, 1, 2)
plt.plot(fs, pgram, color=ocols.blue)
plt.axvline(nm, color=ocols.orange)
plt.axvline(nm+dn, color=ocols.orange)
plt.axvline(nm+2*dn, color=ocols.orange)
plt.axvline(nm+3*dn, color=ocols.orange)
plt.subplot(3, 1, 4)
plt.plot(x, rv, 'k.')
plt.show()
