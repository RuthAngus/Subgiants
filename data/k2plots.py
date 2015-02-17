import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
reb, fbt = plot_params()
from colours import plot_colours
ocols = plot_colours()
import scaling_relations as sr
import sin_tests as st
from astero_modelling import gen_freqs
import scipy.signal as sps

m = 1.29
r = 4.903
teff = 4973
nm = sr.nu_max(m, r, teff)/1e3  # Hz
dn = sr.delta_nu(m, r)/1e6  # Hz

DIR = "/Users/angusr/Python/Subgiants/"
x, y, yerr = np.genfromtxt('%s/data/142091_data1.txt' % DIR, skip_header=5).T
x = (x - x[0])*24*3600

nf = 1
freqs = np.arange(nm-nf*dn, nm+nf*dn, dn)
ws = 2*np.pi*freqs
print len(ws), "frequencies used"
ys, A = st.fit_sine_err(x, y, yerr, ws)

plt.clf()
plt.errorbar(x/60., y, yerr=yerr, **reb)
plt.plot(x/60., ys, color=ocols.blue)
plt.xlabel('$\mathrm{Time~(minutes)}$')
plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
plt.xlim(min(x/60.), max(x/60.))
plt.axhline(-15, color='.6', linestyle='--')
plt.errorbar(x/60., y-ys-15, yerr=yerr, **reb)
plt.ylim(-23, 10)
plt.savefig('%s/paper/142091_1.pdf' % DIR)

print nm, "nm"
fs = np.linspace(.00000001, .001, 10000)
pgram = sps.lombscargle(x, ys, fs*2*np.pi)
plt.clf()
plt.plot(fs, pgram)
plt.savefig("pgram")

rms = np.sqrt(np.mean((y-ys)**2))
print 'rms = ', rms
print np.std(y-ys)

x, y, yerr = np.genfromtxt('%s/data/142091_data2.txt' % DIR, skip_header=5).T
x = (x - x[0])*24*3600
x, y, yerr = x[:-7], y[:-7], yerr[:-7]

ys, A = st.fit_sine_err(x, y, yerr, ws)

plt.clf()
plt.errorbar(x/60., y, yerr=yerr, **reb)
plt.plot(x/60., ys, color=ocols.blue)
plt.xlabel('$\mathrm{Time~(minutes)}$')
plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
plt.axhline(-15, color='.6', linestyle='--')
plt.errorbar(x/60., y-ys-15, yerr=yerr, **reb)
plt.ylim(-23, 10)
plt.savefig('%s/paper/142091_2.pdf' % DIR)
rms = np.sqrt(np.mean((y-ys)**2))
print 'rms = ', rms
print np.std(y-ys)
