"""
cluster_toolkit is a module for computing galaxy cluster models.
"""

import cffi
import glob
import os
import numpy as np

__author__ = "Tom McClintock <tmcclintock@email.arizona.edu>"

cluster_toolkit_dir = os.path.dirname(__file__)
include_dir = os.path.join(cluster_toolkit_dir,'include')
lib_file = os.path.join(cluster_toolkit_dir,'_cluster_toolkit.so')

_ffi = cffi.FFI()
for file_name in glob.glob(os.path.join(include_dir,'*.h')):
    _ffi.cdef(open(file_name).read())
_lib = _ffi.dlopen(lib_file)

def _dcast(x):
    if type(x) is list: x = np.array(x)
    return _ffi.cast('double*', x.ctypes.data)

from . import averaging, bias, boostfactors, deltasigma, density, massfunction, miscentering, xi
