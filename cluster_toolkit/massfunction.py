"""Halo mass function.
"""

import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np

def dndM_at_M(M, k, P, Omega_m, d=1.97, e=1.0, f=0.51, g=1.228):
    """Tinker et al. 2008 appendix C mass function at a given mass. 
    Default behavior is for :math:`M_{200m}` mass definition.

    Args:
        M (float or array like): Mass in Msun/h.
        k (array like): Wavenumbers of the matter power spectrum in h/Mpc comoving.
        P_lin (array like): Linear matter power spectrum in (Mpc/h)^3 comoving.
        Omega_m (float): Matter density fraction.
        d (float; optional): First Tinker parameter. Default is 1.97.
        e (float; optional): Second Tinker parameter. Default is 1.
        f (float; optional): Third Tinker parameter. Default is 0.51.
        g (float; optional): Fourth Tinker parameter. Default is 1.228.

    """
    if type(M) is list or type(M) is np.ndarray:
        dndM = np.zeros_like(M)
        cluster_toolkit._lib.dndM_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), Omega_m, d, e, f, g, _dcast(dndM))
        return dndM
    else:
        return cluster_toolkit._lib.dndM_at_M(M, _dcast(k), _dcast(P), len(k), Omega_m, d, e, f, g)

def _calc_dndM_at_M(M, k, P, om, dndM, d=1.97, e=1.0, f=0.51, g=1.228):
    return cluster_toolkit._lib.dndM_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), om, d, e, f, g, _dcast(dndM))

def _n_in_bin(Mlo, Mhi, Marr, dndM):
    return cluster_toolkit._lib.n_in_bin(_dcast(Marr), _dcast(dndM), len(Marr), Mlo, Mhi)

def _n_in_bins(edges, Marr, dndM):
    N = np.zeros(len(edges)-1)
    cluster_toolkit._lib.n_in_bins(_dcast(Marr), _dcast(dndM), len(Marr), _dcast(edges), len(edges), _dcast(N))
    return N
