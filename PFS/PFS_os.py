import numpy as np
from PFS_selection import read_file
from scaling_synth import amps_and_freqs, fake_rvs
from observing_strategy import os
from selection import exp_time
from PFS_selection import select
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel
import scaling_relations as sc

names, ra, dec, M, vmag, R, T = select()
l = (150 < ra) * (ra < 152)
names, ra, dec, M, vmag, R, T = names[l], ra[l], dec[l], M[l], vmag[l], R[l], \
        T[l]

# calculate dn and nm
dn = sc.delta_nu(M, R)
nm = sc.nu_max(M, R, T)

# load example x values
DIR = "/Users/angusr/Python/Subgiants"
x, y = np.genfromtxt("%s/synthesise/3424541_rvs.txt" % DIR).T
dt = x[1] - x[0]
ndays = 10
xs = np.arange(0, ndays, dt)

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
stype = "GP"

# generate fake rvs
for i, name in enumerate(names):
    if stype=="sine":
        amps_and_freqs(name, dn[i], nm[i]*1e3)
        fake_rvs(xs, name, ndays)
    elif stype=="GP":
        nm = 24*3600*sc.nu_max(M, R, T)*1e-3  # in days-1
        theta = [5, 1, 1, 1./nm[i]]
        k = theta[0] * ExpSquaredKernel(theta[1]) \
                * ExpSine2Kernel(theta[2], theta[3])
        gp = george.GP(k)
        xs = np.linspace(x[0], x[-1], 1000)
        yerr = np.ones_like(xs)*.1
        gp.compute(xs, yerr)
        ys = gp.sample()
        np.savetxt("%s/injections/%s_%s_GP.txt" % (DIR, name, ndays),
                np.transpose((xs, ys)))

for sample in nsamples:
    os(nmins, ndays, ntests, sample, nsim, exptime, names, stype)
