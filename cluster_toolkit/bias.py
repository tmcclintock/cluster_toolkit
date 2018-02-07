"""Halo bias.

"""
import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np
from numpy import ascontiguousarray as ACA

def bias_at_M(M, k, P, Omega_m, delta=200):
    """Tinker et al. 2010 bais at mass M [Msun/h].

    Args:
        M (float or array like): Mass in Msun/h.
        k (array like): Wavenumbers of power spectrum in h/Mpc comoving.
        P (array like): Power spectrum in (Mpc/h)^3 comoving.
        Omega_m (float): Matter density fraction.
        delta (int; optional): Overdensity, default is 200.

    Returns:
        float or array like: Halo bias.

    """
    if type(M) is list or type(M) is np.ndarray:
        bias = np.zeros_like(M)
        M = ACA(M, dtype=np.float64)
        k = ACA(k, dtype=np.float64)
        P = ACA(P, dtype=np.float64)
        cluster_toolkit._lib.bias_at_M_arr(_dcast(M), len(M), delta, _dcast(k), _dcast(P), len(k), Omega_m, _dcast(bias))
        return bias
    else:
        return cluster_toolkit._lib.bias_at_M(M, delta, _dcast(k), _dcast(P), len(k), Omega_m)

def bias_at_R(R, k, P, delta=200):
    """Tinker 2010 bais at mass M [Msun/h] corresponding to radius R [Mpc/h comoving].

    Args:
        R (float or array like): Lagrangian radius in Mpc/h comoving.
        k (array like): Wavenumbers of power spectrum in h/Mpc comoving.
        P (array like): Power spectrum in (Mpc/h)^3 comoving.
        delta (int; optional): Overdensity, default is 200.

    Returns:
        float or array like: Halo bias.

    """
    if type(R) is list or type(R) is np.ndarray:
        bias = np.zeros_like(R)
        cluster_toolkit._lib.bias_at_R_arr(_dcast(R), len(R), delta, _dcast(k), _dcast(P), len(k), _dcast(bias))
        return bias
    else:
        return cluster_toolkit._lib.bias_at_R(R, delta, _dcast(k), _dcast(P), len(k))
    
def bias_at_nu(nu, delta=200):
    """Tinker 2010 bais at peak height nu.

    Args:
        nu (float or array like): Peak height.
        delta (int; optional): Overdensity, default is 200.

    Returns:
        float or array like: Halo bias.

    """
    if type(nu) is list or type(nu) is np.ndarray:
        bias = np.zeros_like(nu)
        cluster_toolkit._lib.bias_at_nu_arr(_dcast(nu), len(nu), delta, _dcast(bias))
        return bias
    else:
        return cluster_toolkit._lib.bias_at_nu(nu, delta)

def sigma2_at_M(M, k, P, Omega_m):
    """RMS variance in top hat sphere of lagrangian radius R [Mpc/h comoving] corresponding to a mass M [Msun/h] of linear power spectrum.

    Args:
        M (float or array like): Mass in Msun/h.
        k (array like): Wavenumbers of power spectrum in h/Mpc comoving.
        P (array like): Power spectrum in (Mpc/h)^3 comoving.
        Omega_m (float): Omega_matter, matter density fraction.

    Returns:
        float or array like: RMS variance of top hat sphere.

    """
    if type(M) is list or type(M) is np.ndarray:
        s2 = np.zeros_like(M)
        cluster_toolkit._lib.sigma2_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), Omega_m, _dcast(s2))
        return s2
    else:
        return cluster_toolkit._lib.sigma2_at_M(M, _dcast(k), _dcast(P), len(k), Omega_m)

def sigma2_at_R(R, k, P):
    """RMS variance in top hat sphere of radius R [Mpc/h comoving] of linear power spectrum.

    Args:
        R (float or array like): Radius in Mpc/h comoving.
        k (array like): Wavenumbers of power spectrum in h/Mpc comoving.
        P (array like): Power spectrum in (Mpc/h)^3 comoving.

    Returns:
        float or array like: RMS variance of a top hat sphere.

    """
    if type(R) is list or type(R) is np.ndarray:
        s2 = np.zeros_like(R)
        cluster_toolkit._lib.sigma2_at_R_arr(_dcast(R), len(R), _dcast(k), _dcast(P), len(k), _dcast(s2))
        return s2
    else:
        return cluster_toolkit._lib.sigma2_at_R(R, _dcast(k), _dcast(P), len(k))

def nu_at_M(M, k, P, Omega_m):
    """Peak height of top hat sphere of lagrangian radius R [Mpc/h comoving] corresponding to a mass M [Msun/h] of linear power spectrum.

    Args:
        M (float or array like): Mass in Msun/h.
        k (array like): Wavenumbers of power spectrum in h/Mpc comoving.
        P (array like): Power spectrum in (Mpc/h)^3 comoving.
        Omega_m (float): Omega_matter, matter density fraction.

    Returns:
        nu (float or array like): Peak height.

    """
    if type(M) is list or type(M) is np.ndarray:
        nu = np.zeros_like(M)
        cluster_toolkit._lib.nu_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), Omega_m, _dcast(nu))
        return nu
    else:
        return cluster_toolkit._lib.nu_at_M(M, _dcast(k), _dcast(P), len(k), Omega_m)

def nu_at_R(R, k, P):
    """Peak height of top hat sphere of radius R [Mpc/h comoving] of linear power spectrum.

    Args:
        R (float or array like): Radius in Mpc/h comoving.
        k (array like): Wavenumbers of power spectrum in h/Mpc comoving.
        P (array like): Power spectrum in (Mpc/h)^3 comoving.

    Returns:
        float or array like: Peak height.

    """
    if type(R) is list or type(R) is np.ndarray:
        nu = np.zeros_like(R)
        cluster_toolkit._lib.nu_at_R_arr(_dcast(R), len(R), _dcast(k), _dcast(P), len(k), _dcast(nu))
        return nu
    else:
        return cluster_toolkit._lib.nu_at_R(R, _dcast(k), _dcast(P), len(k))

def _calc_sigma2_at_R(R, k, P, s2):
    """Direct call to vectorized version of RMS variance in top hat sphere of radius R [Mpc/h comoving] of linear power spectrum.

    """
    cluster_toolkit._lib.sigma2_at_R_arr(_dcast(R), len(R), _dcast(k), _dcast(P), len(k), _dcast(s2))
    return

def _calc_sigma2_at_M(M, k, P, Omega_m, s2):
    """Direct call to vectorized version of RMS variance in top hat sphere of lagrangian radius R [Mpc/h comoving] corresponding to a mass M [Msun/h] of linear power spectrum.

    """
    cluster_toolkit._lib.sigma2_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), Omega_m, _dcast(s2))
    return


def _calc_nu_at_R(R, k, P, nu):
    """Direct call to vectorized version of peak height of R.

    """
    cluster_toolkit._lib.nu_at_R_arr(_dcast(R), len(R), _dcast(k), _dcast(P), len(k), _dcast(nu))
    return

def _calc_nu_at_M(M, k, P, Omega_m, nu):
    """Direct call to vectorized version of peak height of M.

    """
    cluster_toolkit._lib.nu_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), Omega_m, _dcast(nu))
    return

def _calc_bias_at_R(R, k, P, bias, delta=200):
    """Direct call to vectorized version of Tinker 2008 bias at R.

    """
    cluster_toolkit._lib.bias_at_R_arr(_dcast(R), len(R), delta, _dcast(k), _dcast(P), len(k), _dcast(bias))
    return

def _calc_bias_at_M(M, k, P, Omega_m, bias, delta=200):
    """Direct call to vectorized version of Tinker 2008 bias at M.

    """
    cluster_toolkit._lib.bias_at_M_arr(_dcast(M), len(M), delta, _dcast(k), _dcast(P), len(k), Omega_m, _dcast(bias))
    return

def _bias_at_nu_FREEPARAMS(nu,A,a,B,b,C,c, delta=200):
    """A special function used only for quickly computing best fit parameters
    for the halo bias models.
    """
    bias = np.zeros_like(nu)
    cluster_toolkit._lib.bias_at_nu_arr_FREEPARAMS(_dcast(nu), len(nu), delta, A,a,B,b,C,c, _dcast(bias))
    return bias
