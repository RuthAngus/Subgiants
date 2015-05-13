import numpy as np
import matplotlib.pyplot as plt
import george
from george.kernels import ExpSine2Kernel, ExpSquaredKernel, CosineKernel
import scaling_relations as sc
import scipy.signal as sps

D = "/Users/angusr/Python"
ndays = 10
xs = np.linspace(0, ndays, 2000)
yerr = np.ones_like(xs)*1.

# load Luan's subgiants
data = np.genfromtxt("%s/Subgiants/proposal/sample_luan.out" % D,
                     skip_header=1, dtype=str).T
kids = data[0]
data = np.genfromtxt("%s/Subgiants/proposal/sample_luan.out" % D,
                     skip_header=1).T
T, T_err = data[1], data[2]
M, M_err = data[9], data[10]
R, R_err = data[13], data[14]

# mode lifetime scaling relation
tau = 3.2 * (T/5780.)**-4  # in days
nm = 24*3600*sc.nu_max(M, R, T)*1e-3  # in days-1

for i in range(len(kids)):
    print kids[i], i, "of", len(kids)
    k = 5.*ExpSquaredKernel(tau[i])*CosineKernel(1./nm[i])
    gp = george.GP(k)
    gp.compute(xs, yerr)
    ys = gp.sample()
    np.savetxt("%s/Subgiants/injections/%s_%s_GP.txt" % (D, kids[i], ndays),
               np.transpose((xs, ys)))
