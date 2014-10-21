import numpy as np

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
