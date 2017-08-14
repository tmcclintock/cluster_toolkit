import clusterwl
from ctypes import c_double, c_int, POINTER
import numpy as np

def dcast(x):
    return clusterwl._ffi.cast('double*', x.ctypes.data)

def xi_nfw_at_R(R, M, c, om, delta=200):
    if type(R) is list or np.ndarray:
        xi = np.zeros_like(R)
        clusterwl._lib.calc_xi_nfw(dcast(R), len(R), M, c, delta, om, dcast(xi))
        return xi
    else:
        return clusterwl._lib.xi_nfw_at_R(R, M, c, delta, om)

def xi_mm_at_R(R, k, P, N=200, step=0.005):
    if type(R) is list or np.ndarray:
        xi = np.zeros_like(R)
        clusterwl._lib.calc_xi_mm(dcast(R), len(R), dcast(k), dcast(P), len(k), dcast(xi), N, step)
        return xi
    return clusterwl._lib.xi_mm_at_R(R, dcast(k), dcast(P), len(k), N, step)

def xi_2halo(bias, xi_mm):
    NR = len(xi_mm)
    xi = np.zeros_like(xi_mm)
    clusterwl._lib.calc_xi_2halo(NR, bias, dcast(xi_mm), dcast(xi))
    return xi

def xi_hm(xi_1halo, xi_2halo):
    NR = len(xi_1halo)
    xi = np.zeros_like(xi_1halo)
    clusterwl._lib.calc_xi_hm(NR, dcast(xi_1halo), dcast(xi_2halo), dcast(xi))
    return xi

def calc_xi_nfw(R, M, c, om, xi, delta=200):
    clusterwl._lib.calc_xi_nfw(dcast(R), len(R), M, c, delta, om, dcast(xi))
    return

def calc_xi_mm(R, k, P, xi, N=200, step=0.005):
    clusterwl._lib.calc_xi_mm(dcast(R), len(R), dcast(k), dcast(P), len(k), dcast(xi), N, step)
    return

def calc_xi_2halo(bias, xi_mm, xi_2halo):
    NR = len(xi_mm)
    clusterwl._lib.calc_xi_2halo(NR, bias, dcast(xi_mm), dcast(xi_2halo))
    return

def calc_xi_hm(xi_1halo, xi_2halo, xi_hm):
    NR = len(xi_1halo)
    clusterwl._lib.calc_xi_hm(NR, dcast(xi_1halo), dcast(xi_2halo), dcast(xi_hm))
    return
