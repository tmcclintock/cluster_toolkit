import cluster_toolkit
from ctypes import c_double, c_int, POINTER
import numpy as np

def _dcast(x):
    return cluster_toolkit._ffi.cast('double*', x.ctypes.data)

def sigma2_at_R(R, k, P):
    """RMS variance in top hat sphere of radius R [Mpc/h comoving] of linear power spectrum

    Args:
    R (float or array like): Radius in Mpc/h comoving
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving

    Returns:
    sigma2 (float or array like): RMS variance

    """
    if type(R) is list or np.ndarray:
        s2 = np.zeros_like(R)
        cluster_toolkit._lib.sigma2_at_R_arr(_dcast(R), len(R), _dcast(k), _dcast(P), len(k), _dcast(s2))
        return s2
    else:
        return cluster_toolkit._lib.sigma2_at_R(R, _dcast(k), _dcast(P), len(k))

def sigma2_at_M(M, k, P, om):
    """RMS variance in top hat sphere of lagrangian radius R [Mpc/h comoving] corresponding to a mass M [Msun/h] of linear power spectrum

    Args:
    M (float or array like): Mass in Msun/h
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving
    om (float): Omega_matter, matter density fraction

    Returns:
    sigma2 (float or array like): RMS variance

    """
    if type(M) is list or np.ndarray:
        s2 = np.zeros_like(M)
        cluster_toolkit._lib.sigma2_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), om, _dcast(s2))
        return s2
    else:
        return cluster_toolkit._lib.sigma2_at_M(M, _dcast(k), _dcast(P), len(k), om)

def _calc_sigma2_at_R(R, k, P, s2):
    """Direct call to vectorized version of RMS variance in top hat sphere of radius R [Mpc/h comoving] of linear power spectrum

    Args:
    R (array like): Radius in Mpc/h comoving
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving
    sigma2 (array like): RMS variance, populated with the result

    """
    cluster_toolkit._lib.sigma2_at_R_arr(_dcast(R), len(R), _dcast(k), _dcast(P), len(k), _dcast(s2))
    return

def _calc_sigma2_at_M(M, k, P, om, s2):
    """Direct call to vectorized version of RMS variance in top hat sphere of lagrangian radius R [Mpc/h comoving] corresponding to a mass M [Msun/h] of linear power spectrum

    Args:
    M (array like): Mass in Msun/h
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving
    om (float): Omega_matter, matter density fraction
    sigma2 (array like): RMS variance, populated with the result

    """
    cluster_toolkit._lib.sigma2_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), om, _dcast(s2))
    return

def nu_at_R(R, k, P):
    """Peak height of top hat sphere of radius R [Mpc/h comoving] of linear power spectrum

    Args:
    R (float or array like): Radius in Mpc/h comoving
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving

    Returns:
    nu (float or array like): Peak height

    """
    if type(R) is list or np.ndarray:
        nu = np.zeros_like(R)
        cluster_toolkit._lib.nu_at_R_arr(_dcast(R), len(R), _dcast(k), _dcast(P), len(k), _dcast(nu))
        return nu
    else:
        return cluster_toolkit._lib.nu_at_R(R, _dcast(k), _dcast(P), len(k))

def nu_at_M(M, k, P, om):
    """Peak height of top hat sphere of lagrangian radius R [Mpc/h comoving] corresponding to a mass M [Msun/h] of linear power spectrum

    Args:
    M (float or array like): Mass in Msun/h
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving
    om (float): Omega_matter, matter density fraction

    Returns:
    nu (float or array like): Peak height

    """
    if type(M) is list or np.ndarray:
        nu = np.zeros_like(M)
        cluster_toolkit._lib.nu_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), om, _dcast(nu))
        return nu
    else:
        return cluster_toolkit._lib.nu_at_M(M, _dcast(k), _dcast(P), len(k), om)

def _calc_nu_at_R(R, k, P, nu):
    """Direct call to vectorized version of peak height of R

    Args:
    R (array like): Radius in Mpc/h comoving
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving
    nu (array like): RMS variance, populated with the result

    """
    cluster_toolkit._lib.nu_at_R_arr(_dcast(R), len(R), _dcast(k), _dcast(P), len(k), _dcast(nu))
    return

def _calc_nu_at_M(M, k, P, om, nu):
    """Direct call to vectorized version of peak height of M

    Args:
    R (array like): Radius in Mpc/h comoving
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving
    om (float): Matter density fraction
    nu (array like): RMS variance, populated with the result

    """
    cluster_toolkit._lib.nu_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), om, _dcast(nu))
    return

def bias_at_nu(nu, delta=200):
    """Tinker 2008 bais at peak height nu

    Args:
    nu (float or array like): Peak height
    delta (int; optional): Overdensity, default is 200

    Returns:
    bias (float or array like): Bias

    """
    if type(nu) is list or np.ndarray:
        bias = np.zeros_like(nu)
        cluster_toolkit._lib.bias_at_nu_arr(_dcast(nu), len(nu), delta, _dcast(k), _dcast(P), len(k), _dcast(bias))
        return bias
    else:
        return cluserwl._lib.bias_at_nu(nu, delta)

def bias_at_R(R, k, P, delta=200):
    """Tinker 2008 bais at mass M [Msun/h] corresponding to radius R [Mpc/h comoving]

    Args:
    R (float or array like): Lagrangian radius in Mpc/h comoving
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving
    delta (int; optional): Overdensity, default is 200

    Returns:
    bias (float or array like): Bias

    """
    if type(R) is list or np.ndarray:
        bias = np.zeros_like(R)
        cluster_toolkit._lib.bias_at_R_arr(_dcast(R), len(R), delta, _dcast(k), _dcast(P), len(k), _dcast(bias))
        return bias
    else:
        return cluster_toolkit._lib.bias_at_R(R, delta, _dcast(k), _dcast(P), len(k))

def bias_at_M(M, k, P, om, delta=200):
    """Tinker 2008 bais at mass M [Msun/h]

    Args:
    M (float or array like): Mass in Msun/h
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving
    om (float): Matter density fraction
    delta (int; optional): Overdensity, default is 200

    Returns:
    bias (float or array like): Bias

    """
    if type(M) is list or type(M) is np.ndarray:
        bias = np.zeros_like(M)
        cluster_toolkit._lib.bias_at_M_arr(_dcast(M), len(M), delta, _dcast(k), _dcast(P), len(k), om, _dcast(bias))
        return bias
    else:
        return cluster_toolkit._lib.bias_at_M(M, delta, _dcast(k), _dcast(P), len(k), om)

def _calc_bias_at_R(R, k, P, bias, delta=200):
    """Direct call to vectorized version of Tinker 2008 bias at R

    Args:
    R (array like): Radius in Mpc/h comoving
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving
    om (float): Matter density fraction
    bias (array like): Bias, populated with the result
    """
    cluster_toolkit._lib.bias_at_R_arr(_dcast(R), len(R), delta, _dcast(k), _dcast(P), len(k), _dcast(bias))
    return

def _calc_bias_at_M(M, k, P, om, bias, delta=200):
    """Direct call to vectorized version of Tinker 2008 bias at M

    Args:
    M (array like): Mass in Msun/h
    k (array like): Wavenumbers of power spectrum in h/Mpc comoving
    P (array like): Power spectrum in (Mpc/h)^3 comoving
    om (float): Matter density fraction
    bias (array like): Bias, populated with the result
    """
    cluster_toolkit._lib.bias_at_M_arr(_dcast(M), len(M), delta, _dcast(k), _dcast(P), len(k), om, _dcast(bias))
    return
