import cluster_toolkit
from ctypes import c_double, c_int
import numpy as np

def _dcast(x):
    if type(x) is list: x = np.array(x)
    return cluster_toolkit._ffi.cast('double*', x.ctypes.data)
