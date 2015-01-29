import idlsave
from astropy.coordinates import SkyCoord
from astropy import units as u
from selection import exp_time
import numpy as np

def read_file(filename):
    s = idlsave.read(filename)
    data = s.sme
    return data.name, data.mass, data.vmag, data.ra, data.dec, \
            data.radius_iso, data.teff

def select():
    name, m, v, ra, dec, R, T = read_file("spocsiv.dat")
    ramin = SkyCoord('11 40 31 +00 00 00', 'icrs', unit=(u.hourangle, u.deg))
    ramax = SkyCoord('14 56 44 +00 00 00', 'icrs', unit=(u.hourangle, u.deg))
#     ramin, ramax = 175, 224

    # Dec < +10
    # RA = observable in Feb
    # 6.5 < V  < 8.5
    # Range of masses from 1.0 Msun to 1.8 Msun

    # RA_sun = 20:40:31 - 21:56:44
    # RA_midnight_UT = 8:40:31 - 9:56:44
    # timezone = 4
    # RA_midnight_chile = 12:40:31 - 13:56:44

#     ramin, ramax = 115, 165
#     l = (dec<10) * (6.5<v) * (v<8.5) * (1<m) * (m<1.8) \
#             * (ramin<ra) * (ra<ramax)

    ramin, ramax = 90, 165
    l = (dec<20) * (6.5<v) * (v<8.5) * (1<m) * (m<1.8) \
            * (ramin<ra) * (ra<ramax)

    return name[l], ra[l], dec[l], m[l], v[l], R[l], T[l]

if __name__ == "__main__":

    name, ra, dec, m, v, R, T = select()

    # calculate exposure times
    PFS = (8., 88.)  # 88 secs on a 8 mag star for S/N = 188
    Vg, expg = PFS
    exptime = exp_time(v, Vg, expg)

    etimes = []
    for i in range(len(name)):
#         print name[i], ra[i], dec[i], v[i]
        print "M = ", m[i], "R = ", R[i], "T = ", T[i]
#         print "exptime = ", exptime[i], "\n"
        etimes.append(exptime[i])

    np.savetxt("exptimes.txt", np.array(etimes))

# JAN 27, 2014       20:40:31       -18:20: 4
# JAN 28, 2014       20:44:38       -18: 4:19
# JAN 29, 2014       20:48:46       -17:48:10
# JAN 30, 2014       20:52:53       -17:31:43
# JAN 31, 2014       20:56:58       -17:15: 1
# FEB  1, 2014       21: 1: 3       -16:57:56
# FEB  2, 2014       21: 5: 7       -16:40:33
# FEB  3, 2014       21: 9:10       -16:22:58
# FEB  4, 2014       21:13:12       -16: 5: 0
# FEB  5, 2014       21:17:14       -15:46:47
# FEB  6, 2014       21:21:14       -15:28:21
# FEB  7, 2014       21:25:14       -15: 9:35
# FEB  8, 2014       21:29:14       -14:50:34
# FEB  9, 2014       21:33:11       -14:31:23
# FEB 10, 2014       21:37: 9       -14:11:52
# FEB 11, 2014       21:41: 6       -13:52: 7
# FEB 12, 2014       21:45: 1       -13:32:13
# FEB 13, 2014       21:48:56       -13:12: 1
# FEB 14, 2014       21:52:51       -12:51:35
# FEB 15, 2014       21:56:44       -12:31: 3
