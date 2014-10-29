import numpy as np
import matplotlib.pyplot as plt
from scaling_relations import delta_nu, nu_max
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel, WhiteKernel
from BGdata import BetaGem
BG = BetaGem()
from colors import plot_colors
ocols = plot_colors()

# Load Beta Gem data
t, rv, rv_err = BG.rvHJD-BG.rvHJD[0], BG.rv/np.median(BG.rv), \
        BG.rv_err/np.median(BG.rv)
nm_days = BG.nm*3600*24
P = 1./nm_days

# train your GP
theta = [1.**2, 2.**2, 1., P]
k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[3])
gp = george.GP(k)
gp.compute(t, rv_err)
gp.optimize(t, rv, rv_err)
xs = np.linspace(min(t), max(t), 1000)
mu, cov = gp.predict(rv, xs)
std = np.sqrt(np.diag(cov))

plt.clf()
plt.errorbar(t, rv, yerr=rv_err, fmt='k.', capsize=0, ecolor='.8')
plt.plot(xs, mu, ocols.blue)
plt.fill_between(mu+std, mu-std, color=ocols.blue, alpha=.5)
plt.savefig('BGGP_training')
