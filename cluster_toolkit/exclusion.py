"""Correlation functions with halo exclusion.
"""
import cluster_toolkit
from cluster_toolkit import _dcast as dc
import numpy as np

def xi_hm_exclusion_at_r(radii, Mass, conc,
                         rt, beta, Ma, ca, Mb, cb,
                         bias, ximm, Omega_m, delta=200):
    """Halo-matter correlation function with halo exclusion incorporated.

    Args:
        radii (float or array-like): Radii of the profile in Mpc/h
        Mass (float): Mass in Msun/h
        conc (float): concentration of the 1-halo profile
        rt (float): truncation radius in Mpc/h
        beta (float): width of the truncation distribution (erfc) in Mpc/h
        Ma (float): Mass of ...
        ca (float): concentration of ...
        Mb (float): Mass of ...
        cb (float): concentration of ...
        bias (float): halo bias at large scales
        ximm (float): matter correlation function
        Omega_m (float): Matter density fraction
        delta (float): halo overdensity; default is 200

    Returns:
        float or array-like: profile at radii

    """
    if type(radii) is list or type(radii) is np.ndarray:
        xihm = np.zeros_like(radii)
        cluster_toolkit._lib.xihm_exclusion_at_r_arr(dc(radii), len(radii),
                                                     Mass, conc, rt, beta, Ma, ca,
                                                     Mb, cb, bias, dc(ximm), delta,
                                                     Omega_m, dc(xihm))
        return xihm
    else:
        r = np.array(radii)
        xihm = np.zeros_like(r)
        cluster_toolkit._lib.xihm_exclusion_at_r_arr(dc(r), len(r), Mass, conc, rt,
                                                     beta, Ma, ca, Mb, cb,
                                                     bias, dc(ximm), delta, Omega_m,
                                                     dc(xihm))
        return xihm
    

def theta_at_r(radii, rt, beta):
    """Truncation function.
    
    Args:
        radii (float or array-like): Radii of the profile in Mpc/h
        rt (float): truncation radius in Mpc/h
        beta (float): width of the truncation distribution (erfc) in Mpc/h

    Returns:
        float or array-like: Truncation function

    """
    if type(radii) is list or type(radii) is np.ndarray:
        th = np.zeros_like(radii)
        cluster_toolkit._lib.theta_erfc_at_r_arr(dc(radii), len(radii), rt, beta, dc(th))
        return th
    else:
        return cluster_toolkit._lib.theta_erfc_at_r(radii, rt, beta)
