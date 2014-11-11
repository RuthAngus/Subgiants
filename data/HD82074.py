import numpy as np
import matplotlib.pyplot as plt
import idlsave
from rc_params import plot_params
from colors import plot_colors
params = plot_params()
ocols = plot_colors()

class HD(object):

    def __init__(self):

        DIR = '/Users/angusr/Python/Subgiants'
        data = idlsave.read("%s/data/vsthd82074.dat" % DIR)

        t0 = 16700
        t = data.cf3['JD'] - t0
        rv = data.cf3['MNVEL']
        rv_err = data.cf3['ERRVEL']

        xlims = [(1.6*24, 1.85*24), (85, 95), (110, 120),
                 (135, 145), (155, 170), (205, 215), (230, 240),
                 (255, 265), (275, 290)]  # hours

        for i in range(len(xlims)):
            l = (xlims[i][0]/24. < t) * (t < xlims[i][1]/24.)
            if i==0:
                t0, rv0, rv_err0 = t[l], rv[l], rv_err[l]
            elif i==1:
                t1, rv1, rv_err1 = t[l], rv[l], rv_err[l]
            elif i==2:
                t2, rv2, rv_err2 = t[l], rv[l], rv_err[l]
            elif i==3:
                t3, rv3, rv_err3 = t[l], rv[l], rv_err[l]
            elif i==4:
                t4, rv4, rv_err4 = t[l], rv[l], rv_err[l]
            elif i==5:
                t5, rv5, rv_err5 = t[l], rv[l], rv_err[l]
            elif i==6:
                t6, rv6, rv_err6 = t[l], rv[l], rv_err[l]
            elif i==7:
                t7, rv7, rv_err7 = t[l], rv[l], rv_err[l]
            elif i==8:
                t8, rv8, rv_err8 = t[l], rv[l], rv_err[l]

        self.t = t
        self.rv = rv
        self.rv_err = rv_err
        self.t0 = t0
        self.rv0 = rv0
        self.rv_err0 = rv_err0
        self.t1 = t1
        self.rv1 = rv1
        self.rv_err1 = rv_err1
        self.t2 = t2
        self.rv2 = rv2
        self.rv_err2 = rv_err2
        self.t3 = t3
        self.rv3 = rv3
        self.rv_err3 = rv_err3
        self.t4 = t4
        self.rv4 = rv4
        self.rv_err4 = rv_err4
        self.t5 = t5
        self.rv5 = rv5
        self.rv_err5 = rv_err5
        self.t6 = t6
        self.rv6 = rv6
        self.rv_err6 = rv_err6
        self.t7 = t7
        self.rv7 = rv7
        self.rv_err7 = rv_err7
        self.t8 = t8
        self.rv8 = rv8
        self.rv_err8 = rv_err8

if __name__ == "__main__":
    plt.clf()
    plt.errorbar(t*24, rv, yerr=rv_err, ecolor='.8', capsize=0, fmt='k.')
    plt.xlabel('$\mathrm{JD-%s}$' % t0)
    plt.ylabel('$\mathrm{RV}~(ms^{-1})$')
    # plt.xlim(1.6*24, 1.85*24)
    # plt.xlim(85, 95)
    # plt.xlim(110, 120)
    # plt.xlim(135, 145)
    # plt.xlim(155, 170)
    plt.xlim(205, 215)
    # plt.xlim(170, 300)
    plt.ylim(-30, 30)
    plt.savefig("%s/figures/HD82074" % DIR)
    plt.show()
