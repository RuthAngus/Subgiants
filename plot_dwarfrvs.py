import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
import idlsave
import glob
from rc_params import plot_params
from colors import plot_colors
params = plot_params()
ocols = plot_colors()

# load star names
DIR = '/Users/angusr/Python/Subgiants'
stars = np.genfromtxt('%s/data/star_names.txt' % DIR)
rv_files = glob.glob('%s/data/vst*.dat' % DIR)

for star in stars:
    data = idlsave.read("%s/data/vst%s.dat" % (DIR, int(star)))

    t = data.cf3['JD']
    rv = data.cf3['MNVEL']
    rv_err = data.cf3['ERRVEL']

    plt.clf()
    plt.errorbar(t, rv, yerr=rv_err, capsize=0, ecolor='.8', fmt='k.')
    plt.xlabel("$\mathrm{JD}$")
    plt.ylabel("$\mathrm{RV~(ms^{-1})}$")
    plt.savefig('%s/figures/%s' % (DIR, int(star)))
