import clusterwl
from ctypes import c_double, c_int, POINTER

def dcast(x):
    return clusterwl._ffi.cast('double*', x.ctypes.data)

def Sigma_mis_single_at_R(R, Rs, Sigma, M, conc, om, Rmis, delta=200):
    return clusterwl._lib.Sigma_mis_single_at_R(R, dcast(Rs), dcast(Sigma), len(Rs), M, conc, delta, om, Rmis)

def calc_Sigma_mis_single_at_R(R, Rs, Sigma, M, conc, om, Rmis, Sigma_mis, delta=200):
    return clusterwl._lib.Sigma_mis_single_at_R_arr(dcast(R), len(R), dcast(Rs), dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, dcast(Sigma_mis))

def Sigma_mis_at_R(R, Rs, Sigma, M, conc, om, Rmis, delta=200, kernel="gaussian"):
    if kernel == "gaussian": integrand_switch = 0
    elif kernel == "exponential": integrand_switch = 1
    return clusterwl._lib.Sigma_mis_at_R(R, dcast(Rs), dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, integrand_switch)

def calc_Sigma_mis_at_R(R, Rs, Sigma, M, conc, om, Rmis, Sigma_mis, delta=200, kernel="gaussian"):
    if kernel == "gaussian": integrand_switch = 0
    elif kernel == "exponential": integrand_switch = 1
    return clusterwl._lib.Sigma_mis_at_R_arr(dcast(R), len(R), dcast(Rs), dcast(Sigma), len(Rs), M, conc, delta, om, Rmis, integrand_switch, dcast(Sigma_mis))

def DeltaSigma_mis_at_R(R, Rs, Sigma):
    return clusterwl._lib.DeltaSigma_mis_at_R(R, dcast(Rs), dcast(Sigma), len(Rs))

def calc_DeltaSigma_mis_at_R(R, Rs, Sigma, DeltaSigma_mis):
    return clusterwl._lib.DeltaSigma_mis_at_R_arr(dcast(R), len(R), dcast(Rs), dcast(Sigma), len(R), dcast(DeltaSigma_mis))
