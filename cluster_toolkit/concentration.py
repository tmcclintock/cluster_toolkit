"""Halo concentration.

"""
import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np

def concentration_at_M(Mass, k, P, Omega_m, delta=200, Mass_type="mean"):
    """Concentration of the NFW profile at mass M [Msun/h].
    Only implemented relation at the moment is Diemer & Kravtsov (2015).

    Note: only single concentrations at a time are allowed at the moment.
    
    Args:
        Mass (float): Mass in Msun/h.
        k (array like): Wavenumbers of power spectrum in h/Mpc comoving.
        P (array like): Power spectrum in (Mpc/h)^3 comoving.
        Omega_m (float): Matter density fraction.
        delta (int; optional): Overdensity, default is 200.
        Mass_type(string; optional); Defines either Mcrit or Mmean. Default is mean. Choose "crit" for Mcrit. Other values will raise an exception.

    Returns:
        float: NFW concentration.

    """
    if Mass_type is "mean":
        return cluster_toolkit._lib.DK15_concentration_at_Mmean(Mass, _dcast(k), _dcast(P), len(k), delta, Omega_m)
    elif Mass_type is "crit":
        return cluster_toolkit._lib.DK15_concentration_at_Mcrit(Mass, _dcast(k), _dcast(P), len(k), delta, Omega_m)
    else:
        raise Exception("ConcentrationError: must choose either 'mean' or 'crit'")


