"""cluster_toolkit

A module for computing galaxy cluster models.
"""
import os, cffi, glob

cluster_toolkit_dir = os.path.dirname(__file__)
include_dir = os.path.join(cluster_toolkit_dir,'include')
lib_file = os.path.join(cluster_toolkit_dir,'_cluster_toolkit.so')

_ffi = cffi.FFI()
for file_name in glob.glob(os.path.join(include_dir,'*.h')):
    _ffi.cdef(open(file_name).read())
_lib = _ffi.dlopen(lib_file)

from .xi import *
from .bias import *
from .deltasigma import *
from .miscentering import *
from .averaging import *
from .massfunction import *
from .boostfactors import *
from .density import *
