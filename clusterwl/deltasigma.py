import clusterwl
from ctypes import c_double, c_int, POINTER

def dcast(x):
    return clusterwl._ffi.cast('double*', x.ctypes.data)

def Sigma_nfw_at_R(R, M, conc, om, delta=200):
    return clusterwl._lib.Sigma_nfw_at_R(R, M, conc, delta, om)

def calc_Sigma_nfw_at_R(R, M, conc, om, Sigma, delta=200):
    return clusterwl._lib.Sigma_nfw_at_R_arr(dcast(R), len(R), M, conc, delta, om, dcast(Sigma))

def Sigma_at_R(R, Rxi, xi, M, conc, om, delta=200):
    return clusterwl._lib.Sigma_at_R_full(R, dcast(Rxi), dcast(xi), len(Rxi), M, conc, delta, om)

def calc_Sigma_at_R(R, Rxi, xi, M, conc, om, Sigma, delta=200):
    return clusterwl._lib.Sigma_at_R_full_arr(dcast(R), len(R), dcast(Rxi), dcast(xi), len(Rxi), M, conc, delta, om, dcast(Sigma))

def DeltaSigma_at_R(R, Rs, Sigma, M, conc, om, delta=200):
    return clusterwl._lib.DeltaSigma_at_R(R, dcast(Rs), dcast(Sigma), len(Rs), M, conc, delta, om)

def calc_DeltaSigma_at_R(R, Rs, Sigma, M, conc, om, DeltaSigma, delta=200):
    return clusterwl._lib.DeltaSigma_at_R_arr(dcast(R), len(R), dcast(Rs), dcast(Sigma), len(Rs), M, conc, delta, om, dcast(DeltaSigma))
