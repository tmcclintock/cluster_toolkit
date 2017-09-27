import clusterwl
from ctypes import c_double, c_int, POINTER
import numpy as np

def _dcast(x):
    return clusterwl._ffi.cast('double*', x.ctypes.data)

def Sigma_mis_single_at_R(R, Rs, Sigma, M, conc, om, Rmis, delta=200):
    """Miscentered surface mass density [Msun h/pc^2 comoving] of a profile miscentered by an amount Rmis Mpc/h comoving

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
        clusterwl._lib.Sigma_mis_single_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, _dcast(Sigma_mis))
        return Sigma_mis
    else:
        return clusterwl._lib.Sigma_mis_single_at_R(R, _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, Rmis)

def _calc_Sigma_mis_single_at_R(R, Rs, Sigma, M, conc, om, Rmis, Sigma_mis, delta=200):
    return clusterwl._lib.Sigma_mis_single_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, _dcast(Sigma_mis))

def Sigma_mis_at_R(R, Rs, Sigma, M, conc, om, Rmis, delta=200, kernel="gaussian"):
    if kernel == "gaussian": integrand_switch = 0
    elif kernel == "exponential": integrand_switch = 1
    return clusterwl._lib.Sigma_mis_at_R(R, _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, integrand_switch)

def calc_Sigma_mis_at_R(R, Rs, Sigma, M, conc, om, Rmis, Sigma_mis, delta=200, kernel="gaussian"):
    if kernel == "gaussian": integrand_switch = 0
    elif kernel == "exponential": integrand_switch = 1
    return clusterwl._lib.Sigma_mis_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, integrand_switch, _dcast(Sigma_mis))

def DeltaSigma_mis_at_R(R, Rs, Sigma):
    return clusterwl._lib.DeltaSigma_mis_at_R(R, _dcast(Rs), _dcast(Sigma), len(Rs))

def calc_DeltaSigma_mis_at_R(R, Rs, Sigma, DeltaSigma_mis):
    return clusterwl._lib.DeltaSigma_mis_at_R_arr(_dcast(R), len(R), _dcast(Rs), _dcast(Sigma), len(R), _dcast(DeltaSigma_mis))
