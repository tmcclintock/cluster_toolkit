"""Galaxy cluster boost factors, also known as membership dilution models.

"""
import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np

def boost_nfw_at_R(R, B0, R_scale):
    """NFW boost factor model.

    Args:
        R (float or array like): Distances on the sky in the same units as R_scale. Mpc/h comoving suggested for consistency with other modules.
        B0 (float): NFW profile amplitude.
        R_scale (float): NFW profile scale radius.

    Returns:
        float or array like: NFW boost factor profile; B = (1-fcl)^-1.

    """
    if type(R) is list or type(R) is np.ndarray:
        boost = np.zeros_like(R)
        cluster_toolkit._lib.boost_nfw_at_R_arr(_dcast(R), len(R), B0, R_scale, _dcast(boost))
        return boost
    else:
        return cluster_toolkit._lib.boost_nfw_at_R(R, B0, R_scale)

def boost_powerlaw_at_R(R, B0, R_scale, alpha):
    """Power law boost factor model.

    Args:
        R (float or array like): Distances on the sky in the same units as R_scale. Mpc/h comoving suggested for consistency with other modules.
        B0 (float): Boost factor amplitude.
        R_scale (float): Power law scale radius.
        alpha (float): Power law exponent.

    Returns:
        float or array like: Power law boost factor profile; B = (1-fcl)^-1.

    """
    if type(R) is list or type(R) is np.ndarray:
        boost = np.zeros_like(R)
        cluster_toolkit._lib.boost_powerlaw_at_R_arr(_dcast(R), len(R), B0, R_scale, alpha, _dcast(boost))
        return boost
    else:
        return cluster_toolkit._lib.boost_powerlaw_at_R(R, B0, R_scale, alpha)
