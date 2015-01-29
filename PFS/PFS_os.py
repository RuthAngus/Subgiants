import numpy as np
from PFS_selection import read_file
from scaling_relations import delta_nu, nu_max
from scaling_synth import amps_and_freqs, fake_rvs
from observing_strategy import os
from selection import exp_time
from PFS_selection import select

names, ra, dec, M, vmag, R, T = select()
l = (150 < ra) * (ra < 152)
names, ra, dec, M, vmag, R, T = names[l], ra[l], dec[l], M[l], vmag[l], R[l], \
        T[l]

# calculate dn and nm
dn = delta_nu(M, R)
nm = nu_max(M, R, T)

# load example x values
DIR = "/Users/angusr/Python/Subgiants"
x, y = np.genfromtxt("%s/synthesise/3424541_rvs.txt" % DIR).T
dt = x[1] - x[0]
ndays = 10
xs = np.arange(0, ndays, dt)

# # generate fake rvs
# for i, name in enumerate(names):
#     amps_and_freqs(name, dn[i], nm[i]*1e3)
#     fake_rvs(xs, name, ndays)

# calculate exposure times
PFS = (8., 88.)  # 88 secs on a 8 mag star for S/N = 188
Vg, expg = PFS
exptime = exp_time(vmag, Vg, expg)

nmins = 2  # interval between observations in minutes
ndays = 10  # number of nights observed. if this is greater than the no.
# of nights in the simulated data, the rms will increase
ntests = 50  # number of changes in interval
nsamples = [2, 3, 5]  # number of observations per night
nsim = 1  # number of simulations
exptime = 100
stype = "sine"

for sample in nsamples:
    os(nmins, ndays, ntests, sample, nsim, exptime, names, stype)
