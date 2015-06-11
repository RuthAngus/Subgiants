import idlsave
from astropy.coordinates import SkyCoord
from astropy import units as u
from selection import exp_time
import numpy as np
import scaling_relations as sr

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

    ramin, ramax = 190, 230
    l = (dec<30) * (6.5<v) * (v<8.5) * (1<m) * (m<1.8) \
            * (ramin<ra) * (ra<ramax)

    # August run
    ramin, ramax = 300, 340
    l = (dec<30) * (6.5<v) * (v<8.5) * (1<m) * (m<1.8) \
            * (ramin<ra) * (ra<ramax)

    print len(name), "targets"
    return name[l], ra[l], dec[l], m[l], v[l], R[l], T[l]

if __name__ == "__main__":

    name, ra, dec, m, v, R, T = select()

    print "hd185 nm = ", sr.nu_max(1.99, 5.35, 5016)*1e3
    print "hd185 dn = ", sr.delta_nu(1.99, 5.35)

    # calculate exposure times
    PFS = (8., 88.)  # 88 secs on a 8 mag star for S/N = 188
    Vg, expg = PFS
    exptime = exp_time(v, Vg, expg)

    etimes = []
    for i in range(len(name)):
        print name[i], ra[i], dec[i], v[i]
        print "M = ", m[i], "R = ", R[i], "T = ", T[i]
        nm = sr.nu_max(m[i], R[i], T[i])*1e3
        print "nm = ", nm
        print "dn = ", sr.delta_nu(m[i], R[i])
        print "1/2 period = ", 1./(nm * 1e-6) / 60. / 2., "mins"
        print "exptime = ", exptime[i], "\n"
        etimes.append(exptime[i])

    print len(name), "targets"
    np.savetxt("exptimes.txt", np.array(etimes))
