"""Galaxy cluster density profiles.

"""
import cluster_toolkit
from ctypes import c_double, c_int, POINTER
import numpy as np

def _dcast(x):
    if type(x) is list: x = np.array(x)
    return cluster_toolkit._ffi.cast('double*', x.ctypes.data)

def rho_nfw_at_R(R, M, c, om, delta=200):
    """NFW halo density profile.

    Args:
        R (float or array like): 3d distances from halo center in Mpc/h comoving.
        M (float): Mass in Msun/h.
        c (float): Concentration.
        om (float): Omega_matter, matter fraction of the density.
        delta (int; optional): Overdensity, default is 200.

    Returns:
        float or array like: NFW halo density profile in Msun h^2/Mpc^3 comoving.

    """
    if type(R) is list or type(R) is np.ndarray:
        rho = np.zeros_like(R)
        cluster_toolkit._lib.calc_rho_nfw(_dcast(R), len(R), M, c, delta, om, _dcast(rho))
        return rho
    else:
        return cluster_toolkit._lib.rho_nfw_at_R(R, M, c, delta, om)

def rho_einasto_at_R(R, M, rs, alpha, om, delta=200, rhos=-1.):
    """Einasto halo density profile.

    Args:
        R (float or array like): 3d distances from halo center in Mpc/h comoving.
        M (float): Mass in Msun/h; not used if rhos is specified.
        rhos (float): Scale density in Msun h^2/Mpc^3 comoving; optional.
        rs (float): Scale radius.
        alpha (float): Profile exponent.
        om (float): Omega_matter, matter fraction of the density.
        delta (int): Overdensity, default is 200.

    Returns:
        float or array like: Einasto halo density profile in Msun h^2/Mpc^3 comoving.

    """
    if type(R) is list or type(R) is np.ndarray:
        rho = np.zeros_like(R)
        cluster_toolkit._lib.calc_rho_einasto(_dcast(R), len(R), M, rhos, rs, alpha, delta, om, _dcast(rho))
        return rho
    else:
        return cluster_toolkit._lib.rho_einasto_at_R(R, M, rhos, rs, alpha, delta, om)

__all__ = [
    'rho_nfw_at_R',
    'rho_einasto_at_R'
]
