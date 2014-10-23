import numpy as np
from BGdata import BetaGem
import matplotlib.pyplot as plt
from fit_rv import fit_rv
from scipy.signal import periodogram, lombscargle
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel, WhiteKernel, CosineKernel
from scipy.misc import derivative
from rc_params import plot_params
from colors import plot_colors
ocol = plot_colors()

def deriv(mu, xs):
    return (mu[1:] - mu[:-1]) / (xs[1:] - xs[:-1])

def predict(xs):
    return gp.predict(y, xs)[0]

# load flux data
BG = BetaGem()
x = BG.fHJD
y, yerr = BG.flux-np.median(BG.flux), BG.flux_err
xs = np.linspace(min(x), max(x), 200)

# compute periodogram
plt.clf()
freqs = np.linspace(.1, 30, 1000)
pgram = lombscargle(x, y, freqs)
plt.plot(1./freqs, pgram)
# plt.plot(xs, mu, color=ocol.blue)
l = pgram==max(pgram)
w = freqs[l]
print w, 'angular freq'
P = 2*np.pi / w
print
P = 1./(8.29606681747e-05)

# # plot sine fit
# plt.clf()
# plt.ylabel("flux")
# plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
# sinefit = np.sin(w*xs)
# plt.plot(xs, sinefit, 'r-')
# plt.show()
# plt.savefig('BGflux')

# theta = [.16, .1 ** 2, .5, .01, .5]
theta = [.16, .1 ** 2, .5, .01, P]

# GP
k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[4])
k += WhiteKernel(theta[3])
gp = george.GP(k)
gp.compute(x, yerr)
print "initial likelihood = ", gp.lnlikelihood(y, quiet=True)
# result = gp.optimize(x, y, yerr)
# print 'result', result

# predicting
mu = gp.predict(y, xs)[0]
derivs = deriv(mu, xs)

# plot initial guess
plt.clf()
plt.subplot(2,1,1)
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
plt.plot(xs, mu, color=ocol.orange)
plt.subplot(2,1,2)
plt.plot(xs[1:-2], derivs[1:-1], color=ocol.blue)
# plt.show()
# plt.savefig('init_rv_fit')

# load rv data
rvx, rv, rverr = BG.rvHJD, BG.rv, BG.rv_err

# # plot rv derivs
# plt.clf()
# plt.subplot(3,1,1)
# plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
# plt.plot(xs, mu, color=ocol.orange)
# plt.ylabel("Flux")

plt.subplot(2,1,1)
plt.ylabel("Flux")
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
plt.plot(xs, mu, color=ocol.orange)
plt.xlim(x[0], x[0]+.2)

l = 38
rvx, rv, rverr = rvx[:l], rv[:l], rverr[:l]
k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[4])
k += WhiteKernel(theta[3])
gp = george.GP(k)
gp.compute(rvx, rverr)
rvxs = np.linspace(min(rvx), max(rvx), 100)
mu = gp.predict(rv, rvxs)[0]

plt.subplot(2,1,2)
plt.ylabel("RV (m/s)")
plt.xlabel("Time")
plt.errorbar(rvx, rv, yerr=rverr, fmt='k.',
             capsize=0, ecolor='.8')
result = gp.optimize(rvx, rv, rverr)
print result[0][-1]
mu = gp.predict(rv, rvxs)[0]
plt.plot(rvxs, mu, ocol.pink)
plt.xlim(rvx[0], rvx[0]+0.2)
plt.savefig('GPRV')
plt.show()
