"""Correlation functions for matter and halos.

"""
import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np

def xi_nfw_at_R(R, M, c, om, delta=200):
    """NFW halo profile correlation function.

    Args:
        R (float or array like): 3d distances from halo center in Mpc/h comoving
        M (float): Mass in Msun/h
        c (float): Concentration
        om (float): Omega_matter, matter fraction of the density
        delta (int; optional): Overdensity, default is 200

    Returns:
        float or array like: NFW halo profile.

    """
    if type(R) is list or type(R) is np.ndarray:
        xi = np.zeros_like(R)
        cluster_toolkit._lib.calc_xi_nfw(_dcast(R), len(R), M, c, delta, om, _dcast(xi))
        return xi
    else:
        return cluster_toolkit._lib.xi_nfw_at_R(R, M, c, delta, om)

def xi_einasto_at_R(R, M, rs, alpha, om, delta=200, rhos=-1.):
    """Einasto halo profile.

    Args:
        R (float or array like): 3d distances from halo center in Mpc/h comoving
        M (float): Mass in Msun/h; not used if rhos is specified
        rhos (float): Scale density in Msun h^2/Mpc^3 comoving; optional
        rs (float): Scale radius
        alpha (float): Profile exponent
        om (float): Omega_matter, matter fraction of the density
        delta (int): Overdensity, default is 200

    Returns:
        float or array like: Einasto halo profile.

    """
    if type(R) is list or type(R) is np.ndarray:
        xi = np.zeros_like(R)
        cluster_toolkit._lib.calc_xi_einasto(_dcast(R), len(R), M, rhos, rs, alpha, delta, om, _dcast(xi))
        return xi
    else:
        return cluster_toolkit._lib.xi_einasto_at_R(R, M, rhos, rs, alpha, delta, om)

def xi_mm_at_R(R, k, P, N=200, step=0.005):
    """Matter-matter correlation function.

    Args:
        R (float or array like): 3d distances from halo center in Mpc/h comoving
        k (array like): Wavenumbers of power spectrum in h/Mpc comoving
        P (array like): Matter power spectrum in (Mpc/h)^3 comoving
        N (int; optional): Quadrature step count, default is 200
        step (float; optional): Quadrature step size, default is 5e-3

    Returns:
        float or array like: Matter-matter correlation function

    """
    if type(R) is list or type(R) is np.ndarray:
        xi = np.zeros_like(R)
        cluster_toolkit._lib.calc_xi_mm(_dcast(R), len(R), _dcast(k), _dcast(P), len(k), _dcast(xi), N, step)
        return xi
    return cluster_toolkit._lib.xi_mm_at_R(R, _dcast(k), _dcast(P), len(k), N, step)

def xi_2halo(bias, xi_mm):
    """2-halo term in halo-matter correlation function

    Args:
        bias (float): Halo bias
        xi_mm (float or array like): Matter-matter correlation function

    Returns:
        float or array like: 2-halo term in halo-matter correlation function

    """
    NR = len(xi_mm)
    xi = np.zeros_like(xi_mm)
    cluster_toolkit._lib.calc_xi_2halo(NR, bias, _dcast(xi_mm), _dcast(xi))
    return xi

def xi_hm(xi_1halo, xi_2halo):
    """Halo-matter correlation function

    Args:
        xi_1halo (float or array like): 1-halo term
        xi_2halo (float or array like, same size as xi_1halo): 2-halo term

    Returns:
        float or array like: Halo-matter correlation function

    """

    NR = len(xi_1halo)
    xi = np.zeros_like(xi_1halo)
    cluster_toolkit._lib.calc_xi_hm(NR, _dcast(xi_1halo), _dcast(xi_2halo), _dcast(xi))
    return xi

def _calc_xi_nfw(R, M, c, om, xi, delta=200):
    """Direct call to the vectorized version of xi_nfw(R).

    """
    cluster_toolkit._lib.calc_xi_nfw(_dcast(R), len(R), M, c, delta, om, _dcast(xi))
    return

def _calc_xi_mm(R, k, P, xi, N=200, step=0.005):
    """Direct call to the vectorized version of xi_mm(R).

    """

    cluster_toolkit._lib.calc_xi_mm(_dcast(R), len(R), _dcast(k), _dcast(P), len(k), _dcast(xi), N, step)
    return

def _calc_xi_2halo(bias, xi_mm, xi_2halo):
    """Direct call to the vectorized version of xi_2halo(R).

    """
    NR = len(xi_mm)
    cluster_toolkit._lib.calc_xi_2halo(NR, bias, _dcast(xi_mm), _dcast(xi_2halo))
    return

def _calc_xi_hm(xi_1halo, xi_2halo, xi_hm):
    """Direct call to the vectorized version of xi_hm(R).

    """

    NR = len(xi_1halo)
    cluster_toolkit._lib.calc_xi_hm(NR, _dcast(xi_1halo), _dcast(xi_2halo), _dcast(xi_hm))
    return
