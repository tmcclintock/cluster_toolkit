import clusterwl
import numpy as np
from ctypes import c_double, c_int, POINTER

def dcast(x):
    return clusterwl._ffi.cast('double*', x.ctypes.data)

def dndM_at_M(M, k, P, om, d=1.97, e=1.0, f=0.51, g=1.228):
    return clusterwl._lib.dndM_at_M(M, dcast(k), dcast(P), len(k), om, d, e, f, g)

def calc_dndM_at_M(M, k, P, om, dndM, d=1.97, e=1.0, f=0.51, g=1.228):
    return clusterwl._lib.dndM_at_M_arr(dcast(M), len(M), dcast(k), dcast(P), len(k), om, d, e, f, g, dcast(dndM))

def N_in_bin(Mlo, Mhi, volume, Marr, dndM):
    return clusterwl._lib.N_in_bin(dcast(Marr), dcast(dndM), len(Marr), volume, Mlo, Mhi)

def N_in_bins(edges, volume, Marr, dndM):
    N = np.zeros(len(edges)-1)
    clusterwl._lib.N_in_bins(dcast(Marr), dcast(dndM), len(Marr), volume, dcast(edges), len(edges), dcast(N))
    return N
