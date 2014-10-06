import numpy as np
import matplotlib.pyplot as plt
from exoplanet_hosts import exo_hosts
from rc_params import plot_params
from astropy import constants
mearth = constants.M_jup/constants.M_earth
print constants.M_jup
print constants.M_earth
print mearth
mearth = 317.8

stars = exo_hosts()

l = (stars.mstar>0) * (stars.msini>0)
mstar, msini = stars.mstar[l], stars.msini[l]

plt.clf()
plt.plot(mstar, msini * mearth, 'k.')
plt.xlabel('$M_{star} (M_{\odot})$')
plt.ylabel('$M\sini (M_{Earth})$')
plt.axvspan(1.2, 2.5, facecolor='r', alpha=.2)
plt.ylim(0, 100)
plt.savefig('mstar_msini')
