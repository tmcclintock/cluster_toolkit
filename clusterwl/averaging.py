import clusterwl
from ctypes import c_double, c_int, POINTER

def dcast(x):
    return clusterwl._ffi.cast('double*', x.ctypes.data)

def average_profile_in_bin(Rlow, Rhigh, R, prof):
    return clusterwl._lib.average_profile_in_bin(Rlow, Rhigh, dcat(R), len(R), dcast(prof))

def average_profile_in_bins(Redges, R, prof, ave_prof):
    return clusterwl._lib.average_profile_in_bins(dcast(Redges), len(Redges), dcast(R), len(R), dcast(prof), dcast(ave_prof))
