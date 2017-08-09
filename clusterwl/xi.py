import clusterwl
from ctypes import c_double, c_int, POINTER

def dcast(x):
    return clusterwl._ffi.cast('double*', x.ctypes.data)

def calc_xi_nfw(R, M, c, delta, om, xi):
    NR = len(R)
    clusterwl._lib.calc_xi_nfw(dcast(R), NR, M, c, delta, om, dcast(xi))
    return

def calc_xi_mm(R, k, P, xi, N=200, step=0.005):
    NR = len(R)
    Nk = len(k)
    clusterwl._lib.calc_xi_mm(dcast(R), NR, dcast(k), dcast(P), Nk, dcast(xi), N, step)
    return

def calc_xi_2halo(bias, xi_mm, xi_2halo):
    NR = len(xi_mm)
    clusterwl._lib.calc_xi_2halo(NR, bias, dcast(xi_mm), dcast(xi_2halo))
    return

def calc_xi_hm(xi_1halo, xi_2halo, xi_hm):
    NR = len(xi_1halo)
    clusterwl._lib.calc_xi_hm(NR, dcast(xi_1halo), dcast(xi_2halo), dcast(xi_hm))
    return
