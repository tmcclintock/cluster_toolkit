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
    if isinstance(x, list): x = np.asarray(x, dtype=np.float64, order='C')
    return _ffi.cast('double*', x.ctypes.data)

from . import averaging, bias, boostfactors, concentration, deltasigma, density, exclusion, massfunction, miscentering, peak_height, xi
