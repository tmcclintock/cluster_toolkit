import cluster_toolkit
from ctypes import c_double, c_int, POINTER
import numpy as np

def _dcast(x):
    return cluster_toolkit._ffi.cast('double*', x.ctypes.data)

def Sigma_mis_single_at_R(R, Rs, Sigma, M, conc, om, Rmis, delta=200):
    """Miscentered surface mass density [Msun h/pc^2 comoving] of a profile miscentered by an amount Rmis Mpc/h comoving. Units are Msun h/pc^2 comoving

    Args:
    R (float or array like): Projected radii Mpc/h comoving
    Rs (array like): Projected radii of the centered surface mass density profile
    Sigma (float or array like): Surface mass density Msun h/pc^2 comoving
    M (float): Halo mass Msun/h
    conc (float): concentration
    om (float): Matter density fraction
    Rmis (float): Miscentered distance in Mpc/h comoving
    delta (int; optional): Overdensity, default is 200

    Returns:
    Sigma_mis (float or array like): Miscentered projected surface mass density

    """
    if type(R) is list or type(R) is np.ndarray:
        Sigma_mis = np.zeros_like(R)
        cluster_toolkit._lib.Sigma_mis_single_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, _dcast(Sigma_mis))
        return Sigma_mis
    else:
        return cluster_toolkit._lib.Sigma_mis_single_at_R(R, _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, Rmis)

def _calc_Sigma_mis_single_at_R(R, Rs, Sigma, M, conc, om, Rmis, Sigma_mis, delta=200):
    """Direct call to vectorized miscentered sufrace mass density profile

    Args:
    R (float or array like): Projected radii Mpc/h comoving
    Rs (array like): Projected radii of the centered surface mass density profile
    Sigma (float or array like): Surface mass density Msun h/pc^2 comoving
    M (float): Halo mass Msun/h
    conc (float): concentration
    om (float): Matter density fraction
    Rmis (float): Miscentered distance in Mpc/h comoving
    delta (int; optional): Overdensity, default is 200
    Sigma_mis (float or array like): Miscentered projected surface mass density, populated with the result
    """
    return cluster_toolkit._lib.Sigma_mis_single_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, _dcast(Sigma_mis))

def Sigma_mis_at_R(R, Rs, Sigma, M, conc, om, Rmis, delta=200, kernel="gaussian"):
    """Miscentered surface mass density [Msun h/pc^2 comoving] convolved with a distribution for Rmis. Units are Msun h/pc^2 comoving

    Args:
    R (float or array like): Projected radii Mpc/h comoving
    Rs (array like): Projected radii of the centered surface mass density profile
    Sigma (float or array like): Surface mass density Msun h/pc^2 comoving
    M (float): Halo mass Msun/h
    conc (float): concentration
    om (float): Matter density fraction
    Rmis (float): Miscentered distance in Mpc/h comoving
    delta (int; optional): Overdensity, default is 200
    kernel (string; optional): Kernal for convolution, default is gaussian, options are gaussian or exponential

    Returns:
    Sigma_mis (float or array like): Miscentered projected surface mass density

    """
    if kernel == "gaussian": integrand_switch = 0
    elif kernel == "exponential": integrand_switch = 1
    if type(R) is list or type(R) is np.ndarray:
        Sigma_mis = np.zeros_like(R)
        cluster_toolkit._lib.Sigma_mis_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, integrand_switch, _dcast(Sigma_mis))
        return Sigma_mis
    else:
        return cluster_toolkit._lib.Sigma_mis_at_R(R, _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, integrand_switch)

def _calc_Sigma_mis_at_R(R, Rs, Sigma, M, conc, om, Rmis, Sigma_mis, delta=200, kernel="gaussian"):
    """Direct call to vectorized miscentered surface density profile with convolution. Units are Msun h/pc^2 comoving

    Args:
    R (float or array like): Projected radii Mpc/h comoving
    Rs (array like): Projected radii of the centered surface mass density profile
    Sigma (float or array like): Surface mass density Msun h/pc^2 comoving
    M (float): Halo mass Msun/h
    conc (float): concentration
    om (float): Matter density fraction
    Rmis (float): Miscentered distance in Mpc/h comoving
    delta (int; optional): Overdensity, default is 200
    kernel (string; optional): Kernal for convolution, default is gaussian, options are gaussian or exponential
    Sigma_mis (float or array like): Miscentered projected surface mass density, populated with the result

    """
    if kernel == "gaussian": integrand_switch = 0
    elif kernel == "exponential": integrand_switch = 1
    return cluster_toolkit._lib.Sigma_mis_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, integrand_switch, _dcast(Sigma_mis))

def DeltaSigma_mis_at_R(R, Rs, Sigma):
    """Miscentered excess surface mass density profile at R. Units are Msun h/pc^2 comoving

    Args:
    R (float or array like): Projected radii to evaluate profile
    Rs (array like): Projected radii of miscentered Sigma profile
    Sigma (array like): Miscentered Sigma profile

    Returns:
    DeltaSigma_mis (float array like): Miscentered excess surface mass 

    """
    if type(R) is list or type(R) is np.ndarray:
        DeltaSigma_mis = np.zeros_like(R)
        cluster_toolkit._lib.DeltaSigma_mis_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(R), _dcast(DeltaSigma_mis))
        return DeltaSigma_mis
    else:
        return cluster_toolkit._lib.DeltaSigma_mis_at_R(R, _dcast(Rs), _dcast(Sigma), len(Rs))

def _calc_DeltaSigma_mis_at_R(R, Rs, Sigma, DeltaSigma_mis):
    """Direct call to vectorized miscentered excess surface mass density profile at R. Units are Msun h/pc^2 comoving

    Args:
    R (float or array like): Projected radii to evaluate profile
    Rs (array like): Projected radii of miscentered Sigma profile
    Sigma (array like): Miscentered Sigma profile
    DeltaSigma_mis (float array like): Miscentered excess surface mass, populated with the result

    """
    return cluster_toolkit._lib.DeltaSigma_mis_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(R), _dcast(DeltaSigma_mis))
