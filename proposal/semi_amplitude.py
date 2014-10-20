import numpy as np
from orbit import rv_K
from astropy import constants

def Kepler3(a, Ms):
    return np.sqrt( 4*np.pi**2 / (constants.G*Ms) * a**3 )

def mass_function(m1, m2, incl):
    return m1 * (np.sin(incl))**2 / (m1+m2)**3

def my_semi_amplitude(P, m1, m2, incl):
    return (2*np.pi*G/P)**(1/3.) * m2*np.sin(incl) / m1**(2/3.)

m1 = 1.5 * constants.M_sun
m2 = 1 * constants.M_jup
incl = np.pi/2.
a = 0.5 * constants.au
P = Kepler3(a, m1)
P_days = P/60./60./24.
G = constants.G

print my_semi_amplitude(P, m1, m2, incl), 'Jupiter orbiting 1.5 M_sun at .5 AU'

m2 = .6 * constants.M_jup
print my_semi_amplitude(P, m1, m2, incl), 'Saturn orbiting 1.5 M_sun at .5 AU'

m2 = 17.15 * constants.M_earth
print my_semi_amplitude(P, m1, m2, incl), 'Neptune orbiting 1.5 M_sun at .5 AU'

a = 0.3 * constants.au
P = Kepler3(a, m1)
m2 = 17.15 * constants.M_earth
print my_semi_amplitude(P, m1, m2, incl), 'Neptune orbiting 1.5 M_sun at .3 AU'

a = 0.1 * constants.au
P = Kepler3(a, m1)
m2 = 10 * constants.M_earth
print my_semi_amplitude(P, m1, m2, incl), '10 M_earth orbiting 1.5 M_sun at .1 AU'

