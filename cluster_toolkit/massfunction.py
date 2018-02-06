"""Halo mass function.
"""

import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np

def dndM_at_M(M, k, P, Omega_m, d=1.97, e=1.0, f=0.51, g=1.228):
    """Tinker et al. 2008 appendix C mass function at a given mass. 
    Default behavior is for :math:`M_{200m}` mass definition.

    NOTE: by default, this function is only valid at :math:`z=0`. For use
    at higher redshifts either recompute the parameters yourself, or
    wait for this behavior to be patched.

    Args:
        M (float or array like): Mass in Msun/h.
        k (array like): Wavenumbers of the matter power spectrum in h/Mpc comoving.
        P_lin (array like): Linear matter power spectrum in (Mpc/h)^3 comoving.
        Omega_m (float): Matter density fraction.
        d (float; optional): First Tinker parameter. Default is 1.97.
        e (float; optional): Second Tinker parameter. Default is 1.
        f (float; optional): Third Tinker parameter. Default is 0.51.
        g (float; optional): Fourth Tinker parameter. Default is 1.228.

    Returns:
        float or array like: Mass function :math:`dn/dM`.

    """
    if type(M) is list or type(M) is np.ndarray:
        dndM = np.zeros_like(M)
        cluster_toolkit._lib.dndM_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), Omega_m, d, e, f, g, _dcast(dndM))
        return dndM
    else:
        return cluster_toolkit._lib.dndM_at_M(M, _dcast(k), _dcast(P), len(k), Omega_m, d, e, f, g)

def G_at_M(M, k, P, Omega_m, d=1.97, e=1.0, f=0.51, g=1.228):
    """Tinker et al. 2008 appendix C multiplicity funciton G(M) as 
    a function of mass. Default behavior is for :math:`M_{200m}` mass 
    definition.

    Args:
        M (float or array like): Mass in Msun/h.
        k (array like): Wavenumbers of the matter power spectrum in h/Mpc comoving.
        P_lin (array like): Linear matter power spectrum in (Mpc/h)^3 comoving.
        Omega_m (float): Matter density fraction.
        d (float; optional): First Tinker parameter. Default is 1.97.
        e (float; optional): Second Tinker parameter. Default is 1.
        f (float; optional): Third Tinker parameter. Default is 0.51.
        g (float; optional): Fourth Tinker parameter. Default is 1.228.

    Returns:
        float or array like: Halo multiplicity :math:`G(M)`.
    """
    if type(M) is list or type(M) is np.ndarray:
        G = np.zeros_like(M)
        cluster_toolkit._lib.G_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), Omega_m, d, e, f, g, _dcast(G))
    else:
        return cluster_toolkit._lib.G_at_M(M, _dcast(k), _dcast(P), len(k), Omega_m, d, e, f, g)

def G_at_sigma(sigma, d=1.97, e=1.0, f=0.51, g=1.228):
    """Tinker et al. 2008 appendix C multiplicity funciton G(sigma) as 
    a function of sigma.

    NOTE: by default, this function is only valid at :math:`z=0`. For use
    at higher redshifts either recompute the parameters yourself, or
    wait for this behavior to be patched.

    Args:
        sigma (float or array like): RMS variance of the matter density field.
        d (float; optional): First Tinker parameter. Default is 1.97.
        e (float; optional): Second Tinker parameter. Default is 1.
        f (float; optional): Third Tinker parameter. Default is 0.51.
        g (float; optional): Fourth Tinker parameter. Default is 1.228.

    Returns:
        float or array like: Halo multiplicity :math:`G(\sigma)`.
    """
    if type(sigma) is list or type(sigma) is np.ndarray:
        G = np.zeros_like(sigma)
        cluster_toolkit._lib.G_at_sigma_arr(_dcast(sigma), len(sigma), d, e, f, g, _dcast(G))
        return G
    else:
        return cluster_toolkit._lib.G_at_sigma(sigma, d, e, f, g)

def n_in_bin(Mlo, Mhi, Marr, dndM):
    """Tinker et al. 2008 appendix C binned mass function.

    Args:
        Mlo (float): Lower mass edge.
        Mhi (float): Upper mass edge.
        Marr (array like): Array of locations that dndM has been evaluated at.
        dndM (array like): Array of dndM.

    Returns:
       float: number density of halos in the mass bin.

    """
    return cluster_toolkit._lib.n_in_bin(Mlo, Mhi, _dcast(Marr), _dcast(dndM), len(Marr))

def n_in_bins(edges, Marr, dndM):
    """Tinker et al. 2008 appendix C binned mass function.

    Args:
        edges (array like): Edges of the mass bins.
        Marr (array like): Array of locations that dndM has been evaluated at.
        dndM (array like): Array of dndM.

    Returns:
       numpy.ndarray: number density of halos in the mass bins. Length is :code:`len(edges)-1`.

    """
    n = np.zeros(len(edges)-1)
    cluster_toolkit._lib.n_in_bins(_dcast(edges), len(edges), _dcast(Marr), _dcast(dndM), len(Marr), _dcast(n))
    return n

def _dndM_sigma2_precomputed(M, sigma2, sigma2_top, sigma2_bot, Omega_m, d=1.97, e=1.0, f=0.51, g=1.228):
    dndM = np.zeros_like(M)
    cluster_toolkit._lib.dndM_sigma2_precomputed(_dcast(M), _dcast(sigma2), _dcast(sigma2_top), _dcast(sigma2_bot), len(M), Omega_m, d, e, f, g, _dcast(dndM))
    return dndM
