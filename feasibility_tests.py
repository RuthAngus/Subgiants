import numpy as np
import matplotlib.pyplot as plt
from scaling_relations import nu_max, delta_nu
import astero
data = astero.astero()

# load astero data
m = data.bm
m_errp = data.bm_errp
m_errm = data.bm_errm
r = data.br
r_errp = data.br_errp
r_errm = data.br_errm
logg = data.blogg
logg_errp = data.blogg_errp
logg_errm = data.blogg_errm
teff = data.bteff
teff_err = data.bteff_err

mlower = 1.8
l = m > mlower
ms, ms_errp, ms_errm, rs, rs_errp, rs_errm, teffs, teffs_err = \
        m[l], m_errp[l], m_errm[l], r[l], r_errp[l], r_errm[l], \
        teff[l], teff_err[l]

# print ms, "+", ms_errp, "-", ms_errm
# print rs, "+", rs_errp, "-", rs_errm
# print teffs, "+/-", teffs_err

nm = nu_max(ms, rs, teffs)
nm_errp = nu_max(ms+ms_errp, rs+rs_errp, teffs+teffs_err)
nm_errm = nu_max(ms-ms_errm, rs-rs_errm, teffs-teffs_err)

dn = delta_nu(ms, rs)
dn_errp = delta_nu(ms+ms_errp, rs+rs_errp)
dn_errm = delta_nu(ms-ms_errm, rs-rs_errm)

# print nm, '+', nm_errp-nm, '-', nm-nm_errm
# print dn, '+', dn_errp-dn, '-', dn-dn_errm

print 'Beta Gem'
BGm = 1.91
BGm_err = 0.09
BGr = 8.8
BGr_err = 0.1
BGteff = 4750
BGteff_err = 150
print nu_max(BGm, BGr, BGteff)
print delta_nu(BGm, BGr)
