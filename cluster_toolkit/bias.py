"""Halo bias.

"""
import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np
from .peak_height import *

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
    M = np.asarray(M)
    scalar_input = False
    if M.ndim == 0:
        M = M[None] #makes r 1D
        scalar_input = True
    if M.ndim > 1:
        raise Exception("M cannot be a >1D array.")
    M = np.require(M, dtype=np.float64, requirements=["C"])
    k = np.require(k, dtype=np.float64, requirements=["C"])
    P = np.require(P, dtype=np.float64, requirements=["C"])
    bias = np.zeros_like(M)
    cluster_toolkit._lib.bias_at_M_arr(_dcast(M), len(M), delta, _dcast(k), _dcast(P), len(k), Omega_m, _dcast(bias))
    if scalar_input:
        return np.squeeze(bias)
    return bias

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
    R = np.asarray(R)
    scalar_input = False
    if R.ndim == 0:
        R = R[None] #makes r 1D
        scalar_input = True
    if R.ndim > 1:
        raise Exception("R cannot be a >1D array.")
    R = np.require(R, dtype=np.float64, requirements=["C"])
    k = np.require(k, dtype=np.float64, requirements=["C"])
    P = np.require(P, dtype=np.float64, requirements=["C"])
    bias = np.zeros_like(R)
    cluster_toolkit._lib.bias_at_R_arr(_dcast(R), len(R), delta, _dcast(k), _dcast(P), len(k), _dcast(bias))
    if scalar_input:
        return np.squeeze(bias)
    return bias
    
def bias_at_nu(nu, delta=200):
    """Tinker 2010 bais at peak height nu.

    Args:
        nu (float or array like): Peak height.
        delta (int; optional): Overdensity, default is 200.

    Returns:
        float or array like: Halo bias.

    """
    nu = np.asarray(nu)
    scalar_input = False
    if nu.ndim == 0:
        nu = nu[None] #makes nu 1D
        scalar_input = True
    if nu.ndim > 1:
        raise Exception("nu cannot be a >1D array.")

    bias = np.zeros_like(nu)
    cluster_toolkit._lib.bias_at_nu_arr(_dcast(nu), len(nu), delta, _dcast(bias))
    if scalar_input:
        return np.squeeze(bias)
    return bias

def _bias_at_nu_FREEPARAMS(nu,A,a,B,b,C,c, delta=200):
    """A special function used only for quickly computing best fit parameters
    for the halo bias models.
    """
    bias = np.zeros_like(nu)
    cluster_toolkit._lib.bias_at_nu_arr_FREEPARAMS(_dcast(nu), len(nu), delta, A,a,B,b,C,c, _dcast(bias))
    return bias
