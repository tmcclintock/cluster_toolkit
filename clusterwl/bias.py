import clusterwl
from ctypes import c_double, c_int, POINTER

def dcast(x):
    return clusterwl._ffi.cast('double*', x.ctypes.data)

def sigma2_at_R(R, k, P):
    return clusterwl._lib.sigma2_at_R(R, dcast(k), dcast(P), len(k))

def sigma2_at_M(M, k, P, om):
    return clusterwl._lib.sigma2_at_M(M, dcast(k), dcast(P), len(k), om)

def calc_sigma2_at_R(R, k, P, s2):
    clusterwl._lib.sigma2_at_R_arr(dcast(R), len(R), dcast(k), dcast(P), len(k), dcast(s2))
    return

def calc_sigma2_at_M(M, k, P, om, s2):
    clusterwl._lib.sigma2_at_M_arr(dcast(R), len(R), dcast(k), dcast(P), len(k), om, dcast(s2))
    return

def nu_at_R(R, k, P):
    return clusterwl._lib.nu_at_R(R, dcast(k), dcast(P), len(k))

def nu_at_M(M, k, P, om):
    return clusterwl._lib.nu_at_M(M, dcast(k), dcast(P), len(k), om)

def calc_nu_at_R(R, k, P, nu):
    clusterwl._lib.nu_at_R_arr(dcast(R), len(R), dcast(k), dcast(P), len(k), dcast(nu))
    return

def calc_nu_at_M(M, k, P, om, nu):
    clusterwl._lib.nu_at_M_arr(dcast(M), len(M), dcast(k), dcast(P), len(k), om, dcast(nu))
    return

def bias_at_nu(nu, delta=200):
    return cluserwl._lib.bias_at_nu(nu, delta)

def bias_at_R(R, k, P, delta=200):
    return clusterwl._lib.bias_at_R(R, delta, dcast(k), dcast(P), len(k))

def bias_at_M(M, k, P, om, delta=200):
    return clusterwl._lib.bias_at_M(M, delta, dcast(k), dcast(P), len(k), om)

def calc_bias_at_R(R, k, P, bias, delta=200):
    clusterwl._lib.bias_at_R_arr(dcast(R), len(R), delta, dcast(k), dcast(P), len(k), dcast(bias))
    return

def calc_bias_at_M(M, k, P, om, bias, delta=200):
    clusterwl._lib.bias_at_M_arr(dcast(M), len(M), delta, dcast(k), dcast(P), len(k), om, dcast(bias))
    return
