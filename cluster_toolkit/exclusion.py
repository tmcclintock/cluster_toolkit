"""Correlation functions with halo exclusion.
"""
import cluster_toolkit
from cluster_toolkit import _dcast as dc
import numpy as np

def xi_hm_exclusion_at_r(radii, Mass, conc,
                         rt, beta, Ma, ca, Mb, cb,
                         bias, ximm, Omega_m, delta=200,
                         exclusion_scheme=0):
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
        ximm (float or array-like): matter correlation function. Must have same shape as the radii.
        Omega_m (float): Matter density fraction
        delta (float): halo overdensity; default is 200
        exclusion_scheme (int): halo exclusion scheme; either 0,1,2; default is 0

    Returns:
        float or array-like: exclusion profile at each radii

    """
    if exclusion_scheme not in [0,1,2]:
        raise Exception("Halo exclusion_scheme must be in [0,1,2].")
    lib =  cluster_toolkit._lib
    if type(radii) is list or type(radii) is np.ndarray:
        xihm = np.zeros_like(radii)
        lib.xihm_exclusion_at_r_arr(dc(radii), len(radii),
                                    Mass, conc, rt, beta, Ma, ca,
                                    Mb, cb, bias, dc(ximm), delta,
                                    Omega_m, exclusion_scheme, dc(xihm))
        return xihm
    else:
        r = np.array(radii)
        xihm = np.zeros_like(r)
        if len(r) != len(np.atleast_1d(ximm)):
            raise Exception("xi_hm_exclusion: Supply single xi_mm value for a single radii.")
        lib.xihm_exclusion_at_r_arr(dc(r), len(r), Mass, conc, rt,
                                    beta, Ma, ca, Mb, cb,
                                    bias, dc(np.array(ximm)), delta, Omega_m,
                                    exclusion_scheme, dc(xihm))
        return xihm

def xi_1h_exclusion_at_r(radii, Mass, conc,
                         rt, beta, Omega_m, delta=200):
    """Halo-matter correlation function with halo exclusion incorporated,
    but just the 1-halo term.

    Args:
        radii (float or array-like): Radii of the profile in Mpc/h
        Mass (float): Mass in Msun/h
        conc (float): concentration of the 1-halo profile
        rt (float): truncation radius in Mpc/h
        beta (float): width of the truncation distribution (erfc) in Mpc/h
        Omega_m (float): Matter density fraction
        delta (float): halo overdensity; default is 200

    Returns:
        float or array-like: 1-halo of the exclusion profile at each radii

    """
    lib =  cluster_toolkit._lib
    if type(radii) is list or type(radii) is np.ndarray:
        xi1h = np.zeros_like(radii)
        lib.xi_1h_at_r_arr(dc(radii), len(radii),
                           Mass, conc, rt, beta, delta,
                           Omega_m, dc(xi1h))
        return xi1h
    else:
        return lib.xi_1h_at_r(radii, Mass, conc, rt, beta, delta, Omega_m)

def xi_2h_exclusion_at_r(radii, bias, ximm):
    """Halo-matter correlation function with halo exclusion incorporated,
    but just the 2-halo term.

    Args:
        radii (float or array-like): Radii of the profile in Mpc/h
        bias (float): halo bias at large scales
        ximm (float or array-like): matter correlation function. Must have same shape as the radii.

    Returns:
        float or array-like: 2-halo of the exclusion profile at each radii

    """

    lib =  cluster_toolkit._lib
    if type(radii) is list or type(radii) is np.ndarray:
        xi2h = np.zeros_like(radii)
        lib.xi_2h_at_r_arr(dc(radii), len(radii),
                           bias, dc(ximm), dc(xi2h))
        return xi2h
    else:
        r = np.array(radii)
        xi2h = np.zeros_like(r)
        if len(r) != len(np.atleast_1d(ximm)):
            raise Exception("xi_2h_exclusion: Supply single xi_mm value for a single radii.")
        return lib.xi_2h_at_r(dc(r), len(r), bias, dc(np.array(ximm)))

def xi_correction_at_r(radii, Mass, rt, Ma, ca, Mb, cb, bias, ximm,
                       Omega_m, delta=200, exclusion_scheme=0):
    """Halo-matter correlation function with halo exclusion incorporated.

    Args:
        radii (float or array-like): Radii of the profile in Mpc/h
        Mass (float): Mass in Msun/h
        rt (float): truncation radius in Mpc/h
        Ma (float): Mass of ...
        ca (float): concentration of ...
        Mb (float): Mass of ...
        cb (float): concentration of ...
        bias (float): halo bias at large scales
        ximm (float or array-like): matter correlation function. Must have same shape as the radii.
        Omega_m (float): Matter density fraction
        delta (float): halo overdensity; default is 200
        exclusion_scheme (int): halo exclusion scheme; either 0,1,2; default is 0

    Returns:
        float or array-like: correction term for the exclusion profile at each radii

    """

    if exclusion_scheme not in [0,1,2]:
        raise Exception("Halo exclusion_scheme must be in [0,1,2].")
    lib =  cluster_toolkit._lib
    if type(radii) is list or type(radii) is np.ndarray:
        xic = np.zeros_like(radii)
        lib.xi_correction_at_r_arr(dc(radii), len(radii), Mass, rt, Ma, ca, Mb, cb,
                                   bias, dc(ximm), delta, Omega_m, exclusion_scheme, dc(xic))
        return xic
    else:
        r = np.array(radii)
        xic = np.zeros_like(r)
        if len(r) != len(np.atleast_1d(ximm)):
            raise Exception("xi_correction_exclusion: Supply single xi_mm value for a single radii.")
        return lib.xi_correction_at_r_arr(dc(r), len(r), Mass, rt, Ma, ca, Mb, cb,
                                          bias, dc(ximm), delta, Omega_m, exclusion_scheme, dc(xic))

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
