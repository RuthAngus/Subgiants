import numpy as np
import matplotlib.pyplot as plt
import idlsave
from rc_params import plot_params
from colors import plot_colors
params = plot_params()
ocols = plot_colors()

# load star names
DIR = '/Users/angusr/Python/Subgiants'
stars = np.genfromtxt('%s/data/star_names.txt' % DIR)
rv_files = glob.glob('%s/data/vst*.dat' % DIR)

plt.clf()
for i, star in enumerate(stars):
    data = idlsave.read("%s/data/vst%s.dat" % (DIR, int(star)))

    t = data.cf3['JD']
    rv = data.cf3['MNVEL']
    rv_err = data.cf3['ERRVEL']

    plt.subplot(len(stars), 1, i+1)
    plt.errorbar(t, rv, yerr=rv_err, capsize=0, ecolor='.8', fmt='k.')

plt.savefig("%s/figures/all_stars" % DIR)
