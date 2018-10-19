"""Integrals of the power spectrum. This includes RMS variance of the density field, sigma2, as well as peak neight, nu. These were previously implemented in the bias module, but have been migrated to here.

"""
import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np
from numpy import ascontiguousarray as ACA

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

def dsigma2dM_at_M(M, k, P, Omega_m):
    """Derivative w.r.t. mass of RMS variance in top hat sphere of lagrangian radius R [Mpc/h comoving] corresponding to a mass M [Msun/h] of linear power spectrum.

    Args:
        M (float or array like): Mass in Msun/h.
        k (array like): Wavenumbers of power spectrum in h/Mpc comoving.
        P (array like): Power spectrum in (Mpc/h)^3 comoving.
        Omega_m (float): Omega_matter, matter density fraction.

    Returns:
        float or array like: d/dM of RMS variance of top hat sphere.

    """
    if type(M) is list or type(M) is np.ndarray:
        ds2dM = np.zeros_like(M)
        cluster_toolkit._lib.dsigma2dM_at_M_arr(_dcast(M), len(M), _dcast(k), _dcast(P), len(k), Omega_m, _dcast(ds2dM))
        return ds2dM
    else:
        return cluster_toolkit._lib.dsigma2dM_at_M(M, _dcast(k), _dcast(P), len(k), Omega_m)

    return 0

    
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
