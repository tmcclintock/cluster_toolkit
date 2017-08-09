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
