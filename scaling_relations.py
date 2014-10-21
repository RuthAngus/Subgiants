import numpy as np

def nu_max(M, R, Teff):
    # M and R in Solar units, Teff in K
    # returns mHz
    return M / (R**2 * np.sqrt(Teff/5777.)) * 3.05

def delta_nu(M, R):
    # M and R in Solar units, Teff in K
    # returns uHz
    return np.sqrt(M) * R**(-3./2.) * 134.9
