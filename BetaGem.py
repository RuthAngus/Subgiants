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

# load flux data
BG = BetaGem()
x = BG.fHJD - BG.fHJD[0]
y, yerr = BG.flux-np.median(BG.flux), BG.flux_err
x = np.arange(len(y))
xs = np.linspace(min(x), max(x), 100)

theta = [.16, .1 ** 2, .5, .01, .5]

k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[4])
k += WhiteKernel(theta[3])
gp = george.GP(k)
gp.compute(x, yerr)
print "initial likelihood = ", gp.lnlikelihood(y, quiet=True)

def predict(xs):
    return gp.predict(y, xs)[0]

mu = predict(xs)
derivs = derivative(predict, xs)

# plot initial guess
plt.clf()
plt.subplot(2,1,1)
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
plt.plot(xs, mu, color=ocol.orange)
plt.subplot(2,1,2)
plt.plot(xs, derivs, color=ocol.blue)
plt.show()
# plt.savefig('init_rv_fit')

# compute periodogram
plt.clf()
freqs = np.linspace(.1, 30, 1000)
pgram = lombscargle(x, y, freqs)
plt.plot(1./freqs, pgram)
plt.plot(xs, mu, color=ocol.blue)
l = pgram==max(pgram)
w = freqs[l]
print w, 'angular freq'

# plot sine fit
plt.clf()
plt.ylabel("flux")
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
sinefit = np.sin(w*xs)
plt.plot(xs, sinefit, 'r-')
# plt.show()
# plt.savefig('BGflux')

# load rv data
rvx, rv, rverr = BG.rvHJD - BG.rvHJD[0], BG.rv, BG.rv_err
# diff = rvx[1:] - rvx[:-1]  # find gaps
# l = diff>1e-1

plt.clf()
plt.ylabel("RV (m/s)")
plt.xlabel("Time")
plt.errorbar(rvx, rv, yerr=rverr, fmt='k.',
             capsize=0, ecolor='.8')
# for i in range(len(diff[l])):  # plot gaps
#     plt.axvline(rvx[l][i], color='r')
plt.xlim(.8, 1.1)
plt.show()

# plt.clf()
# plt.ylabel("RV (m/s)")
# plt.xlabel("Time")
# x1 = rvx[0]
# ngaps = len(diff[l])
# for i in range(1, len(diff[l])):  # plot gaps
#     x2 = rvx[l][i]
#     plt.errorbar(rvx[x1:x2], rv[x1:x2], yerr=rverr[x1:x2], fmt='k.',
#                  capsize=0, ecolor='.8')
#     x1 = x2
#     plt.show()
#     raw_input('enter')
