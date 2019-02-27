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
    radii = np.asarray(radii)
    ximm = np.asarray(ximm)
    scalar_input = False
    if radii.ndim == 0:
        radii = radii[None] #makes r 1D
        scalar_input = True
    if radii.ndim > 1:
        raise Exception("radii cannot be a >1D array.")

    if len(radii) != len(ximm):
        raise Exception("len(r) must equal len(ximm)")

    if exclusion_scheme not in [0,1,2]:
        raise Exception("Halo exclusion_scheme must be in [0,1,2].")
    lib =  cluster_toolkit._lib
    xihm = np.zeros_like(radii)
    lib.xihm_exclusion_at_r_arr(dc(radii), len(radii),
                                Mass, conc, rt, beta, Ma, ca,
                                Mb, cb, bias, dc(ximm), delta,
                                Omega_m, exclusion_scheme, dc(xihm))
    if scalar_input:
        return np.squeeze(xihm)
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
    radii = np.asarray(radii)
    scalar_input = False
    if radii.ndim == 0:
        radii = radii[None] #makes r 1D
        scalar_input = True
    if radii.ndim > 1:
        raise Exception("radii cannot be a >1D array.")

    lib =  cluster_toolkit._lib
    xi1h = np.zeros_like(radii)
    lib.xi_1h_at_r_arr(dc(radii), len(radii),
                       Mass, conc, rt, beta, delta,
                       Omega_m, dc(xi1h))
    if scalar_input:
        return np.squeeze(xi1h)
    return xi1h

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
    radii = np.asarray(radii)
    ximm = np.asarray(ximm)
    scalar_input = False
    if radii.ndim == 0:
        radii = radii[None] #makes r 1D
        scalar_input = True
    if radii.ndim > 1:
        raise Exception("radii cannot be a >1D array.")

    if len(radii) != len(ximm):
        raise Exception("len(r) must equal len(ximm)")

    lib =  cluster_toolkit._lib
    xi2h = np.zeros_like(radii)
    lib.xi_2h_at_r_arr(dc(radii), len(radii),
                       bias, dc(ximm), dc(xi2h))
    if scalar_input:
        return np.squeeze(xi2h)
    return xi2h

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
    radii = np.asarray(radii)
    ximm = np.asarray(ximm)
    scalar_input = False
    if radii.ndim == 0:
        radii = radii[None] #makes r 1D
        scalar_input = True
    if radii.ndim > 1:
        raise Exception("radii cannot be a >1D array.")

    if len(radii) != len(ximm):
        raise Exception("len(r) must equal len(ximm)")

    if exclusion_scheme not in [0,1,2]:
        raise Exception("Halo exclusion_scheme must be in [0,1,2].")
    lib =  cluster_toolkit._lib
    xic = np.zeros_like(radii)
    lib.xi_correction_at_r_arr(dc(radii), len(radii), Mass, rt, Ma, ca, Mb, cb,
                               bias, dc(ximm), delta, Omega_m, exclusion_scheme, dc(xic))
    if scalar_input:
        return np.squeeze(xic)
    return xic

def theta_at_r(radii, rt, beta):
    """Truncation function.
    
    Args:
        radii (float or array-like): Radii of the profile in Mpc/h
        rt (float): truncation radius in Mpc/h
        beta (float): width of the truncation distribution (erfc) in Mpc/h

    Returns:
        float or array-like: Truncation function

    """
    radii = np.asarray(radii)
    scalar_input = False
    if radii.ndim == 0:
        radii = radii[None] #makes r 1D
        scalar_input = True
    if radii.ndim > 1:
        raise Exception("radii cannot be a >1D array.")

    th = np.zeros_like(radii)
    cluster_toolkit._lib.theta_erfc_at_r_arr(dc(radii), len(radii), rt, beta, dc(th))
    if scalar_input:
        return np.squeeze(th)
    return th
