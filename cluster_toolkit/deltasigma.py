"""Galaxy cluster shear and magnification profiles also known as DeltaSigma and Sigma, respectively.

"""
import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np

def Sigma_nfw_at_R(R, M, conc, om, delta=200):
    """Surface mass density of an NFW profile [Msun h/pc^2 comoving].

    Args:
        R (float or array like): Projected radii Mpc/h comoving.
        M (float): Halo mass Msun/h.
        conc (float): concentration.
        om (float): Matter density fraction.
        delta (int; optional): Overdensity, default is 200.

    Returns:
        float or array like: Surface mass density Msun h/pc^2 comoving.

    """
    if type(R) is list or type(R) is np.ndarray:
        Sigma = np.zeros_like(R)
        cluster_toolkit._lib.Sigma_nfw_at_R_arr(_dcast(R), len(R), M, conc, delta, om, _dcast(Sigma))
        return Sigma
    else:
        return cluster_toolkit._lib.Sigma_nfw_at_R(R, M, conc, delta, om)

def Sigma_at_R(R, Rxi, xi, M, conc, om, delta=200):
    """Surface mass density given some 3d profile [Msun h/pc^2 comoving].

    Args:
        R (float or array like): Projected radii Mpc/h comoving.
        Rxi (array like): 3D radii of xi_hm Mpc/h comoving.
        xi_hm (array like): Halo matter correlation function.
        M (float): Halo mass Msun/h.
        conc (float): concentration.
        om (float): Matter density fraction.
        delta (int; optional): Overdensity, default is 200.

    Returns:
        float or array like: Surface mass density Msun h/pc^2 comoving.

    """
    if type(R) is list or type(R) is np.ndarray:
        Sigma = np.zeros_like(R)
        cluster_toolkit._lib.Sigma_at_R_full_arr(_dcast(R), len(R), _dcast(Rxi), _dcast(xi), len(Rxi), M, conc, delta, om, _dcast(Sigma))
        return Sigma
    else:
        return cluster_toolkit._lib.Sigma_at_R_full(R, _dcast(Rxi), _dcast(xi), len(Rxi), M, conc, delta, om)


def DeltaSigma_at_R(R, Rs, Sigma, M, conc, om, delta=200):
    """Excess surface mass density given Sigma [Msun h/pc^2 comoving].

    Args:
        R (float or array like): Projected radii Mpc/h comoving.
        Rs (array like): Projected radii of Sigma, the surface mass density.
        Sigma (array like): Surface mass density.
        M (float): Halo mass Msun/h.
        conc (float): concentration.
        om (float): Matter density fraction.
        delta (int; optional): Overdensity, default is 200.

    Returns:
        float or array like: Excess surface mass density Msun h/pc^2 comoving.

    """
    if type(R) is list or type(R) is np.ndarray:
        DeltaSigma = np.zeros_like(R)
        cluster_toolkit._lib.DeltaSigma_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, _dcast(DeltaSigma))
        return DeltaSigma
    else:
        return cluster_toolkit._lib.DeltaSigma_at_R(R, _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om)

def _calc_DeltaSigma_at_R(R, Rs, Sigma, M, conc, om, DeltaSigma, delta=200):
    """Direct call to vectorized excess surface mass density given Sigma [Msun h/pc^2 comoving]

    """
    return cluster_toolkit._lib.DeltaSigma_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, _dcast(DeltaSigma))

def _calc_Sigma_nfw_at_R(R, M, conc, om, Sigma, delta=200):
    """Direct call to vectorized surface mass density of NFW.

    """
    return cluster_toolkit._lib.Sigma_nfw_at_R_arr(_dcast(R), len(R), M, conc, delta, om, _dcast(Sigma))

def _calc_Sigma_at_R(R, Rxi, xi, M, conc, om, Sigma, delta=200):
    """Direct call to vectorized surface mass density given xi_hm

    """
    return cluster_toolkit._lib.Sigma_at_R_full_arr(_dcast(R), len(R), _dcast(Rxi), _dcast(xi), len(Rxi), M, conc, delta, om, _dcast(Sigma))
