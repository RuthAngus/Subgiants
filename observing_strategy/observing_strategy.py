import numpy as np
import matplotlib.pyplot as plt
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel, WhiteKernel
from colours import plot_colours
ocols = plot_colours()
from rc_params import plot_params
reb, fbk = plot_params()
from BGdata import BetaGem
BG = BetaGem()

# Load Beta Gem data
t, rv, rv_err = BG.rvHJD-BG.rvHJD[0], BG.rv/np.median(BG.rv), \
        BG.rv_err/np.median(BG.rv)
nm_days = BG.nm*3600*24
P = 1./nm_days
l = 38
# l = len(t)
t, rv, rv_err = t[:l], rv[:l], rv_err[:l]

# train GP on Beta Gem (without white noise)
theta = [15., 5., 1., P]
k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[3])
gp = george.GP(k)
gp.compute(t, rv_err)
results = gp.optimize(t, rv, rv_err)[0]
print np.exp(results)

# plot BG data and result
xs = np.linspace(min(t), max(t), 1000)
mu, cov = gp.predict(rv, xs)
plt.clf()
plt.errorbar(t, rv, yerr=rv_err, **reb)
plt.plot(xs, mu, color=ocols.blue)
plt.xlabel('$\mathrm{Time~(days)}$')
plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
std = np.sqrt(np.diag(cov))
plt.fill_between(xs, mu+std, mu-std, **fbk)
plt.savefig('BG_trained')

# Generate time series with a GP
xs = np.linspace(0, 10, 48*10) # one point every 1/2 hour
yerr = np.ones_like(xs)*.01

# theta = [15., 5., 1., 2./24.]
theta = results
k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[3])
gp = george.GP(k)
gp.compute(xs, yerr)
np.random.seed(1234)
s = gp.sample(xs)
xs1, xs2, xs3 = xs[::48], xs[2::48], xs[4::48]
ss1, ss2, ss3 = s[::48], s[2::48], s[4::48]
yerrs1, yerrs2, yerrs3 = yerr[::48], yerr[2::48], yerr[4::48]
xss = np.concatenate((xs1, xs2, xs3))
ss = np.concatenate((ss1, ss2, ss3))
yerrs = np.concatenate((yerrs3, yerrs2, yerrs3))
print np.std(ss)

plt.clf()
plt.plot(xs, s, color=ocols.blue)
plt.errorbar(xss, ss, yerr=yerrs, **reb)
plt.xlabel('$\mathrm{Time~(days)}$')
plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
plt.show()
