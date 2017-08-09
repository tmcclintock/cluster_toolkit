import clusterwl
from ctypes import c_double, c_int, POINTER

def dcast(x):
    return clusterwl._ffi.cast('double*', x.ctypes.data)

def Sigma_at_R(R, Rxi, xi, M, conc, om, delta=200):
    return clusterwl._lib.Sigma_at_R(R, dcast(Rxi), dcast(xi), len(Rxi), M, conc, delta, om)

def calc_Sigma_at_R(R, Rxi, xi, M, conc, om, Sigma, delta=200):
    clusterwl._lib.Sigma_at_R_arr(dcast(R), len(R), dcast(Rxi), dcast(xi), len(Rxi), M, conc, delta, om, dcast(Sigma))
    return
