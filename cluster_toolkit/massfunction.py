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
    M = np.asarray(M)
    scalar_input = False
    if M.ndim == 0:
        M = M[None] #makes M 1D
        scalar_input = True
    if M.ndim > 1:
        raise Exception("M cannot be a >1D array.")
    M = np.require(M, dtype=np.float64, requirements=["C"])
    k = np.require(k, dtype=np.float64, requirements=["C"])
    P = np.require(P, dtype=np.float64, requirements=["C"])

    dndM = np.zeros_like(M)
    cluster_toolkit._lib.dndM_at_M_arr(_dcast(M), len(M), _dcast(k),
                                       _dcast(P), len(k), Omega_m,
                                       d, e, f, g, _dcast(dndM))
    if scalar_input:
        return np.squeeze(dndM)
    return dndM

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
    M = np.asarray(M)
    scalar_input = False
    if M.ndim == 0:
        M = M[None] #makes M 1D
        scalar_input = True
    if M.ndim > 1:
        raise Exception("M cannot be a >1D array.")
    M = np.require(M, dtype=np.float64, requirements=["C"])
    k = np.require(k, dtype=np.float64, requirements=["C"])
    P = np.require(P, dtype=np.float64, requirements=["C"])

    G = np.zeros_like(M)
    cluster_toolkit._lib.G_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), Omega_m, d, e, f, g, _dcast(G))
    if scalar_input:
        return np.squeeze(G)
    return G

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
    sigma = np.asarray(sigma)
    scalar_input = False
    if sigma.ndim == 0:
        sigma = sigma[None]
        scalar_input = True
    if sigma.ndim > 1:
        raise Exception("sigma cannot be a >1D array.")

    G = np.zeros_like(sigma)
    cluster_toolkit._lib.G_at_sigma_arr(_dcast(sigma), len(sigma), d, e, f, g, _dcast(G))
    if scalar_input:
        return np.squeeze(G)
    return G

def n_in_bins(edges, Marr, dndM):
    """Tinker et al. 2008 appendix C binned mass function.

    Args:
        edges (array like): Edges of the mass bins.
        Marr (array like): Array of locations that dndM has been evaluated at.
        dndM (array like): Array of dndM.

    Returns:
       numpy.ndarray: number density of halos in the mass bins. Length is :code:`len(edges)-1`.

    """
    edges = np.asarray(edges)
    if edges.ndim == 0:
        edges = edges[None]
    if edges.ndim > 1:
        raise Exception("edges cannot be a >1D array.")

    n = np.zeros(len(edges)-1)
    cluster_toolkit._lib.n_in_bins(_dcast(edges), len(edges), _dcast(Marr), _dcast(dndM), len(Marr), _dcast(n))
    return n

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
    return np.squeeze(n_in_bins([Mlo, Mhi], Marr, dndM))

def _dndM_sigma2_precomputed(M, sigma2, dsigma2dM, Omega_m, d=1.97, e=1.0, f=0.51, g=1.228):
    dndM = np.zeros_like(M)
    cluster_toolkit._lib.dndM_sigma2_precomputed(_dcast(M), _dcast(sigma2), _dcast(dsigma2dM), len(M), Omega_m, d, e, f, g, _dcast(dndM))
    return dndM
