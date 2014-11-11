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
interval = 48  # 96 for 1/4, 48 for 1/2, 24 for hr, 12 for 2 hr
ndays = 10
xs = np.linspace(0, ndays, interval*ndays) # one point every 1/2 hour
yerr = np.ones_like(xs)*.01

# Compute GP prior sample
theta = results
k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[3])
gp = george.GP(k)
gp.compute(xs, yerr)
# np.random.seed(1234)
# np.random.seed(123)
# np.random.seed(12)
np.random.seed(1)

l = interval  # sample every day
f = interval/24  # sample every hour

samples = gp.sample(xs, 10)

for i, s in enumerate(samples):
    # sample every 1/2, 1, 1.5, 2, hours
    ms = np.arange(.5, 10, .5)  # sampling interval in hours
    stds = []
    for m in ms:

        xs1, xs2, xs3 = xs[::l], xs[m::l], xs[2*m::l]
        ss1, ss2, ss3 = s[::l], s[m::l], s[2*m::l]
        yerrs1, yerrs2, yerrs3 = yerr[::l], yerr[m::l], yerr[2*m::l]

        xss = np.vstack((xs1, xs2, xs3))
        ss = np.vstack((ss1, ss2, ss3))
        yerrs = np.vstack((yerrs1, yerrs2, yerrs3))
        xmean = np.mean(xss, axis=0)
        ymean = np.mean(ss, axis=0)

        stds.append(np.std(ymean))

#         plt.clf()
#         plt.plot(xs, s, color=ocols.blue)
#         plt.errorbar(xs1, ss1, yerr=yerrs1, **reb)
#         plt.errorbar(xs2, ss2, yerr=yerrs2, **reb)
#         plt.errorbar(xs3, ss3, yerr=yerrs3, **reb)
#         plt.plot(xmean, ymean, 'o', color=ocols.pink)
#         plt.xlabel('$\mathrm{Time~(days)}$')
#         plt.ylabel('$\mathrm{RV~(ms}^{-1}\mathrm{)}$')
    #     plt.show()

    sdts = np.array(stds)
    plt.clf()
    plt.subplot(2, 1, 1)
    plt.plot(xs, s, color=ocols.blue)
    plt.xlabel("$\mathrm{Time~(days)}$")
    plt.ylabel("$\mathrm{RV~(ms}^{-1}\mathrm{)}$")
    plt.subplot(2, 1, 2)
    plt.plot(ms, stds, color=ocols.orange,
             label="$%s$" % ms[stds==min(stds)][0])
    plt.xlabel("$\mathrm{Time~between~samples~(hours)}$")
    plt.ylabel("$\mathrm{RMS}$")
    plt.subplots_adjust(hspace=.3)
    plt.legend()
    plt.savefig('GP_%s' % i)
