"""Galaxy cluster shear and magnification profiles also known as DeltaSigma and Sigma, respectively.

"""
import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np

def Sigma_nfw_at_R(R, mass, concentration, Omega_m, delta=200):
    """Surface mass density of an NFW profile [Msun h/pc^2 comoving].

    Args:
        R (float or array like): Projected radii Mpc/h comoving.
        mass (float): Halo mass Msun/h.
        concentration (float): concentration.
        Omega_m (float): Matter density fraction.
        delta (int; optional): Overdensity, default is 200.

    Returns:
        float or array like

    """
    R = np.asarray(R)
    scalar_input = False
    if R.ndim == 0:
        R = R[None] #makes R 1D
        scalar_input = True
    if R.ndim > 1:
        raise Exception("R cannot be a >1D array.")

    Sigma = np.zeros_like(R)
    cluster_toolkit._lib.Sigma_nfw_at_R_arr(_dcast(R), len(R), mass,
                                            concentration, delta,
                                            Omega_m, _dcast(Sigma))
    if scalar_input:
        return np.squeeze(Sigma)
    return Sigma

def Sigma_at_R(R, Rxi, xi, mass, concentration, Omega_m, delta=200):
    """Surface mass density given some 3d profile [Msun h/pc^2 comoving].

    Args:
        R (float or array like): Projected radii Mpc/h comoving.
        Rxi (array like): 3D radii of xi_hm Mpc/h comoving.
        xi_hm (array like): Halo matter correlation function.
        mass (float): Halo mass Msun/h.
        concentration (float): concentration.
        Omega_m (float): Matter density fraction.
        delta (int; optional): Overdensity, default is 200.

    Returns:
        float or array like: Surface mass density Msun h/pc^2 comoving.

    """
    R = np.asarray(R)
    scalar_input = False
    if R.ndim == 0:
        R = R[None] #makes R 1D
        scalar_input = True
    if R.ndim > 1:
        raise Exception("R cannot be a >1D array.")
    if np.min(R) < np.min(Rxi):
        raise Exception("Minimum R for Sigma(R) must be >= than min(r) of xi(r).")
    if np.max(R) > np.max(Rxi):
        raise Exception("Maximum R for Sigma(R) must be <= than max(r) of xi(r).")

    Sigma = np.zeros_like(R)
    cluster_toolkit._lib.Sigma_at_R_full_arr(_dcast(R), len(R), _dcast(Rxi),
                                             _dcast(xi), len(Rxi), mass, concentration,
                                             delta, Omega_m, _dcast(Sigma))
    if scalar_input:
        return np.squeeze(Sigma)
    return Sigma

def DeltaSigma_at_R(R, Rs, Sigma, mass, concentration, Omega_m, delta=200):
    """Excess surface mass density given Sigma [Msun h/pc^2 comoving].

    Args:
        R (float or array like): Projected radii Mpc/h comoving.
        Rs (array like): Projected radii of Sigma, the surface mass density.
        Sigma (array like): Surface mass density.
        mass (float): Halo mass Msun/h.
        concentration (float): concentration.
        Omega_m (float): Matter density fraction.
        delta (int; optional): Overdensity, default is 200.

    Returns:
        float or array like: Excess surface mass density Msun h/pc^2 comoving.

    """
    R = np.asarray(R)
    scalar_input = False
    if R.ndim == 0:
        R = R[None] #makes R 1D
        scalar_input = True
    if R.ndim > 1:
        raise Exception("R cannot be a >1D array.")

    if np.min(R) < np.min(Rs):
        raise Exception("Minimum R for DeltaSigma(R) must be "+
                        ">= than min(R) of Sigma(R).")
    if np.max(R) > np.max(Rs):
        raise Exception("Maximum R for DeltaSigma(R) must be "+
                        "<= than max(R) of Sigma(R).")

    DeltaSigma = np.zeros_like(R)
    cluster_toolkit._lib.DeltaSigma_at_R_arr(_dcast(R), len(R), _dcast(Rs),
                                             _dcast(Sigma), len(Rs), mass,
                                             concentration, delta, Omega_m,
                                             _dcast(DeltaSigma))
    if scalar_input:
        return np.squeeze(DeltaSigma)
    return DeltaSigma
