import numpy as np

def nu_max(M, R, Teff):
    # M and R in Solar units, Teff in K
    # returns mHz
    return M / (R**2 * np.sqrt(Teff/5777.)) * 3.05

def delta_nu(M, R):
    # M and R in Solar units, Teff in K
    # returns uHz
    return np.sqrt(M) * R**(-3./2.) * 134.9

def nu_max_alt(g, Teff):
    # returns mHz
    # solar g = 27.94
    # g in solar units
    return g * np.sqrt(Teff/5777.) * 3.05

def delta_nu_alt(rho):
    # rho in Solar units, Teff in K
    # returns uHz
    return np.sqrt(rho) * 134.9
