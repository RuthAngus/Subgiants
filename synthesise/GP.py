import numpy as np
import matplotlib.pyplot as plt
import george
from george.kernels import ExpSine2Kernel, ExpSquaredKernel
import scaling_relations as sc

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

nm = 24*3600*sc.nu_max(M, R, T)*1e-3  # in days-1

for i in range(len(kids)):
    print kids[i], i, "of", len(kids)

    theta = [5, 1, 1, 1./nm[i]]
    k = theta[0]*ExpSquaredKernel(theta[1])*ExpSine2Kernel(theta[2], theta[3])
    gp = george.GP(k)
    gp.compute(xs, yerr)
    ys = gp.sample()

    np.savetxt("%s/Subgiants/injections/%s_%s_GP.txt" % (D, kids[i], ndays),
               np.transpose((xs, ys)))
