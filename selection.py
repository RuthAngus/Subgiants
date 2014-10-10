import numpy as np
from hipparcos import hipp
import matplotlib.pyplot as plt
from rc_params import plot_params
from colors import plot_colors

# RA of the Sun = 0 on the 21st of March, 6 on 31st June, 12 on the 21st September, 18 on 21st Dec
# RA at midnight = 12 on 21st march,  on 31st June, 0 on 21st sept,
# Therefore RA ~ 8 at midnight in Jan and 22 in July.

def flux_ratio(dm):
    return 2.5 ** dm

def exp_time(V, Vg, expg):

    # CHIRON
    # to get 1 m/s precision you need S/N of 200
    # for a V = 6.5 star that's a 300s exposure

    # PFS
    # to get 2 m/s precision you need S/N of ~200 (188)
    # which is ~88 secs for a V = 8 star

    dm = Vg - V
    e_time = np.zeros_like(V)
    l1 = dm >= 0
    e_time[l1] = abs(expg / flux_ratio(dm[l1]))
    l2 = dm < 0
    e_time[l2] = abs(expg / flux_ratio(dm[l2]))
    return e_time

def cuts(ID, M_v, B_V, V, RA, dec):
    RA_sec, dec_as =  RA, dec # RA_dec_conv(RA, dec)
    l1 = (1.8 < M_v) * (M_v < 3.) * (.8 < B_V) * (B_V < 1.1) * (V < Vlim)
    l2 = (M_v[l1] < 2.) * (B_V[l1] > .8)
    l2 = l2==False
    l3 = (RA_mns-RA_lims < RA_sec[l1][l2]) * (RA_sec[l1][l2] < RA_mns+RA_lims) * \
           (lat_as-lat_lim_as < dec_as[l1][l2]) * (dec_as[l1][l2] < lat_as+lat_lim_as)
    return ID[l1][l2][l3], M_v[l1][l2][l3], B_V[l1][l2][l3], V[l1][l2][l3], \
            RA_sec[l1][l2][l3], dec_as[l1][l2][l3]

def m2M(m, plx_as):  # plx in arcsec
    return m + 5*(1 + np.log10(plx_as/1000.))

def RA_dec_conv(RA, dec):
    RA_sec = RA[2, :] + RA[1, :] * 60 + RA[0, :]*60*60
    dec_as = dec[2, :] + dec[1, :] * 60 + dec[0, :]*60*60
    return RA_sec, dec_as

def RA_dec_r(RA_sec, dec_as):
    RA = np.zeros((3, len(RA_sec)))
    dec = np.zeros((3, len(dec_as)))

    RA[0, :] = [int(ra/60./60.) for ra in RA_sec]
    RA[1, :] = [int(ra/60.) for ra in RA_sec]
    RA[2, :] = RA_sec % 60.

    dec[0, :] = [int(ra/60./60.) for ra in dec_as]
    dec[1, :] = [int(ra/60.) for ra in dec_as]
    dec[2, :] = dec_as % 60.

#     RA[1, :] = [int(RA_sec[i] - RA[0, :][i]*60*60)/60./60. \
#             for i in range(len(RA_sec))]
#     RA[2, :] = [int(RA_sec[i] - RA[0, :][i]*60*60 - \
#             RA[1, :][i]*60)/60./60. for i in range(len(RA_sec))]

#     dec[0, :] = [round(d/60./60.) for d in dec_as]
#     dec[1, :] = [int(d/60.) for d in dec_as]
#     dec[2, :] = dec_as % 60.

#     dec[1, :] = [int(dec_as[i] - dec[0, :][i]*60*60)/60./60. \
#             for i in range(len(dec_as))]
#     dec[2, :] = [int(dec_as[i] - dec[0, :][i]*60*60 - \
#             dec[1, :][i]*60)/60./60. for i in range(len(dec_as))]

    return RA, dec

if __name__ == "__main__":

    # load data
    h = hipp()
    V = h.Vmag
    B_V = h.B_V
    plx = h.plx
    hipno = h.hipno
    M_v = m2M(V, plx)  # calc abs mag

    # GLOBAL VARIABLES
    # Coordinates of Las Campanas Observatory:
    # 29.0146 S, 70.6926 W
    lat, lat_lim = -29.0146, 30.
    lat_as, lat_lim_as = lat*60*60, lat_lim*60*60
    RA_mn, RA_lim = 20., 4.
    RA_mns, RA_lims = RA_mn*60*60, RA_lim*60*60
    Vlim = 8.5

    ID, M_v, B_V, V, RA_sec, dec_as = \
            cuts(np.arange(len(M_v)), M_v, B_V, V, h.RA, h.dec)

    print len(RA_sec), 'targets selected'
    plt.clf()
    plt.plot(RA_sec/60./60., dec_as/60./60., 'k.')
    plt.axhline(lat, color='r')
    plt.axhline(lat+lat_lim, color='k', linestyle='--')
    plt.axhline(lat-lat_lim, color='k', linestyle='--')
    plt.axvline(RA_mn, color='b')
    plt.axvline(RA_mn+RA_lim, color='k', linestyle='--')
    plt.axvline(RA_mn-RA_lim, color='k', linestyle='--')
    plt.xlim(0,24)
    plt.ylim(-90, 90)
    plt.xlabel('RA (hours)')
    plt.ylabel('dec (degrees)')
    plt.savefig('RA_dec')

    # make magnitude plot
    plt.clf()
    plt.plot(B_V, M_v, 'k.')
    plt.savefig('Magnitude')

    Vg = 8.  # 6.5 for CHIRON
    expg = 88. # 300. for CHIRON

    # load Luan's list of Subgiants
    ID = np.genfromtxt('subgiant_list.txt', skip_header=2, dtype=str, usecols=1).T

    data = np.genfromtxt('subgiant_list.txt', skip_header=2).T
    RA_sub = data[4:7, :]
    dec_sub = data[7:10, :]
    plx_sub, B_sub, V_sub = data[10], data[11], data[12]
    B_V_sub = B_sub - V_sub

    Mv_sub = m2M(V_sub, plx_sub)

    ID, Mv_sub, B_V_sub, V_sub, RA_sub_sec, dec_sub_as = \
            cuts(ID, Mv_sub, B_V_sub, V_sub, RA_sub, dec_sub)

    print len(Mv_sub), 'subgiants selected'
    for i in range(len(ID)):
        print ID[i].upper(), RA_sub, dec_sub
