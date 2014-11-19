import numpy as np
import matplotlib.pyplot as pl
import triangle
import h5py
import acor
import sys

fname = sys.argv[1]

with h5py.File("samples_%s" %fname, "r") as f:
    samples = f["samples"][:, 50:, :]
nwalkers, n, ndim = samples.shape
flatchain = samples.reshape((-1, ndim))
mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                  zip(*np.percentile(flatchain, [16, 50, 84], axis=0)))
mres = np.array(mcmc_result)[:, 0]
print 'mcmc_result = ', mres
np.savetxt("parameters%s.txt" %fname, np.array(mcmc_result))

pl.clf()
# fig_labels = ['$\log~g$', '$\rho$', '$T_{eff}$']
fig_labels = ['logg', 'rho', 'Teff']
fig = triangle.corner(flatchain, truths=mres, labels=fig_labels)
fig.savefig("triangle%s.png" %fname)

print("Plotting traces")
pl.figure()
for i in range(ndim):
    pl.clf()
    pl.plot(samples[:, :, i].T, 'k-', alpha=0.3)
    pl.savefig("%s%s.png" %(i, fname))
