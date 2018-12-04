"""
Derivatives of halo profiles. Used to plot splashback results.
"""
import cluster_toolkit as ct
from cluster_toolkit import _dcast
import numpy as np
from numpy import ascontiguousarray as ACA

def drho_nfw_dr_at_R(Radii, Mass, conc, Omega_m, delta=200):
    """Derivative of the NFW halo density profile.

    Args:
        Radii (float or array like): 3d distances from halo center in Mpc/h comoving
        Mass (float): Mass in Msun/h
        conc (float): Concentration
        Omega_m (float): Matter fraction of the density
        delta (int; optional): Overdensity, default is 200

    Returns:
        float or array like: derivative of the NFW profile.

    """
    if type(Radii) is list or type(Radii) is np.ndarray:
        drhodr = np.zeros_like(Radii)
        ct._lib.drho_nfw_dr_at_R_arr(_dcast(R), len(R), Mass, conc,
                                     delta, Omega_m, _dcast(drhodr))
        return xi
    else:
        return cluster_toolkit._lib.drho_nfw_dr_at_R(Radii, Mass, conc,
                                                     delta, Omega_m)

