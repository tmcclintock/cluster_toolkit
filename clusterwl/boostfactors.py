import clusterwl
from ctypes import c_double, c_int, POINTER
import numpy as np

def _dcast(x):
    return clusterwl._ffi.cast('double*', x.ctypes.data)

def boost_nfw_at_R(R, B0, Rs):
    """NFW boost factor model.

    Args:
    R (float or array like): Distances on the sky in the same units as Rs. Mpc/h comoving suggested for consistency with other modules.
    B0 (float): NFW profile amplitude.
    Rs (float): NFW profile scale radius.

    Returns:
    boost_nfw (float or array like): NFW boost factor profile. B = (1-fcl)^-1
    """
    if type(R) is list or type(R) is np.ndarray:
        boost = np.zeros_like(R)
        print B0, Rs
        clusterwl._lib.boost_factor_nfw_at_R_arr(_dcast(R), len(R), B0, Rs, _dcast(boost))
        return boost
    else:
        return clusterwl._lib.boost_nfw_at_R(R, B0, Rs)
