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
        float or array like: Surface mass density Msun h/pc^2 comoving.

    """
    if type(R) is list or type(R) is np.ndarray:
        Sigma = np.zeros_like(R)
        cluster_toolkit._lib.Sigma_nfw_at_R_arr(_dcast(R), len(R), mass, concentration, delta, Omega_m, _dcast(Sigma))
        return Sigma
    else:
        return cluster_toolkit._lib.Sigma_nfw_at_R(R, mass, concentration, delta, Omega_m)

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
    if np.min(R) < np.min(Rxi):
        raise Exception("Minimum R for Sigma(R) must be >= than min(r) of xi(r).")
    if np.max(R) > np.max(Rxi):
        raise Exception("Maximum R for Sigma(R) must be <= than max(r) of xi(r).")
    if type(R) is list or type(R) is np.ndarray:
        Sigma = np.zeros_like(R)
        cluster_toolkit._lib.Sigma_at_R_full_arr(_dcast(R), len(R), _dcast(Rxi), _dcast(xi), len(Rxi), mass, concentration, delta, Omega_m, _dcast(Sigma))
        return Sigma
    else:
        return cluster_toolkit._lib.Sigma_at_R_full(R, _dcast(Rxi), _dcast(xi), len(Rxi), mass, concentration, delta, Omega_m)


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
    if np.min(R) <= np.min(Rs):
        raise Exception("Minimum R for DeltaSigma(R) must be > than min(R) of Sigma(R).")
    if np.max(R) >= np.max(Rs):
        raise Exception("Maximum R for DeltaSigma(R) must be < than max(R) of Sigma(R).")
    if type(R) is list or type(R) is np.ndarray:
        DeltaSigma = np.zeros_like(R)
        cluster_toolkit._lib.DeltaSigma_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(Rs), mass, concentration, delta, Omega_m, _dcast(DeltaSigma))
        return DeltaSigma
    else:
        return cluster_toolkit._lib.DeltaSigma_at_R(R, _dcast(Rs), _dcast(Sigma), len(Rs), mass, concentration, delta, Omega_m)

def _calc_DeltaSigma_at_R(R, Rs, Sigma, mass, concentration, Omega_m, DeltaSigma, delta=200):
    """Direct call to vectorized excess surface mass density given Sigma [Msun h/pc^2 comoving]

    """
    return cluster_toolkit._lib.DeltaSigma_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(Rs), mass, concentration, delta, Omega_m, _dcast(DeltaSigma))

def _calc_Sigma_nfw_at_R(R, mass, concentration, Omega_m, Sigma, delta=200):
    """Direct call to vectorized surface mass density of NFW.

    """
    return cluster_toolkit._lib.Sigma_nfw_at_R_arr(_dcast(R), len(R), mass, concentration, delta, Omega_m, _dcast(Sigma))

def _calc_Sigma_at_R(R, Rxi, xi, mass, concentration, Omega_m, Sigma, delta=200):
    """Direct call to vectorized surface mass density given xi_hm

    """
    return cluster_toolkit._lib.Sigma_at_R_full_arr(_dcast(R), len(R), _dcast(Rxi), _dcast(xi), len(Rxi), mass, concentration, delta, Omega_m, _dcast(Sigma))
