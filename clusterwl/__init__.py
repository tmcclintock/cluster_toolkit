"""clsuterwl

You can find docstrings for each of the package contents by inspecting help(clusterwl.PACKAGENAME). The package names include: xi, bias, deltasigma, miscentering, and boostfactors.
"""
import os, cffi, glob

clusterwl_dir = os.path.dirname(__file__)
include_dir = os.path.join(clusterwl_dir,'include')
lib_file = os.path.join(clusterwl_dir,'_clusterwl.so')

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
