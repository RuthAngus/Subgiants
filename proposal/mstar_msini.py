import numpy as np
import matplotlib.pyplot as plt
from exoplanet_hosts import exo_hosts
from rc_params import plot_params
reb, fbt = plot_params()
from colours import plot_colours
cols = plot_colours()
from astropy import constants
import csv

mearth = constants.M_jup/constants.M_earth
mearth = 317.8

with open('/Users/angusr/Python/Subgiants/data/exoplanet_catalogue.txt',
          'rb') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    mp, mp_err, ms, ms_err  = [], [], [], []
    for row in data:
        mp.append(row[0])
        mp_err.append(row[3])
        ms.append(row[4])
        ms_err.append(row[7])
mp.remove(mp[0])
mp.remove(mp[0])
mp_err.remove(mp_err[0])
mp_err.remove(mp_err[0])
ms.remove(ms[0])
ms.remove(ms[0])
ms_err.remove(ms_err[0])
ms_err.remove(ms_err[0])

mp = np.array(mp)
mp_err = np.array(mp_err)
ms = np.array(ms)
ms_err = np.array(ms_err)

x = np.ones(len(mp))
for i in range(len(mp)):
    if len(mp[i])==0:
       x[i] = 0
    if len(mp_err[i])==0:
       x[i] = 0
    if len(ms[i])==0:
       x[i] = 0
    if len(ms_err[i])==0:
       x[i] = 0
l = x==1
mp = np.array(map(float, mp[l]))
ms = np.array(map(float, ms[l]))
mp_err = np.array(map(float, mp_err[l]))
ms_err = np.array(map(float, ms_err[l]))

l = (ms>0) * (mp>0)
ms, ms_err, mp, mp_err = ms[l], ms_err[l], mp[l], mp_err[l]

plt.clf()
plt.errorbar(ms, mp, yerr=ms_err, xerr=mp_err, **reb)
plt.xlabel('$M_{star} (M_{\odot})$')
plt.ylabel('$M\sin i (M_{Earth})$')
# plt.axvspan(1.2, 2.5, facecolor='r', alpha=.2)
# plt.axvspan(1.338, 1.753, facecolor='r', alpha=.2)
plt.axvspan(1.5, 5, facecolor=cols.blue, alpha=.2, edgecolor=None)
plt.xlim(0, 5)
plt.savefig('/Users/angusr/Python/Subgiants/paper/mstar_msini.pdf')
