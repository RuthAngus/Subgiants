import numpy as np
from scaling_relations import delta_nu, nu_max

class BetaGem(object):

    def __init__(self):

        # load rv data
        data = np.genfromtxt("/Users/angusr/Python/Subgiants/data/BetaGem.txt",
                             skip_header=43, skip_footer=1659-845,
                             invalid_raise=False).T
        rvHJD, rv, rv_err = data

        # load photometry
        data = np.genfromtxt("/Users/angusr/Python/Subgiants/data/BetaGem.txt",
                             skip_header=856, invalid_raise=False).T
        fHJD, flux, flux_err = data

        self.rvHJD = rvHJD
        self.fHJD = fHJD
        self.rv = rv
        self.rv_err = rv_err
        self.flux = flux
        self.flux_err = flux_err
        self.m = 1.91
        self.m_err = 0.09
        self.r = 8.8
        self.r_err = 0.1
        self.teff = 4750
        self.teff_err = 150
        self.nm = nu_max(1.91, 8.8, 4750)/1e3
        self.dnu = delta_nu(1.91, 8.8)/1e6
