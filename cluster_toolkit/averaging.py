"""Averaging projected cluster profiles.

"""
import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np

def average_profile_in_bins(Redges, R, prof):
    """Average profile in bins.

    Calculates the average of some projected profile in a 
    radial bins in Mpc/h comoving.

    Args:
        Redges (array like): Array of radial bin edges.
        R (array like): Radii of the profile.
        prof (array like): Projected profile.

    Returns:
        numpy.array: Average profile in bins between the edges provided.

    """
    Redges = np.asarray(Redges)
    if Redges.ndim == 0: #Is scalar
        Redges = r[None] #makes r 1D
        raise Exception("Must supply a left and right edge.")
    if Redges.ndim > 1:
        raise Exception("Redges cannot be a >1D array.")
    if np.min(Redges) < np.min(R):
        raise Exception("Minimum edge must be >= minimum R")
    if np.max(Redges) > np.max(R):
        raise Exception("Maximum edge must be <= maximum R")
    
    ave_prof = np.zeros(len(Redges)-1)
    cluster_toolkit._lib.average_profile_in_bins(_dcast(Redges), len(Redges),
                                                 _dcast(R), len(R), _dcast(prof),
                                                 _dcast(ave_prof))
    if len(Redges) == 2:#Single bin
        return np.squeeze(ave_prof)
    return ave_prof

def average_profile_in_bin(Rlow, Rhigh, R, prof):
    """Average profile in a bin.

    Calculates the average of some projected profile in a 
    radial bin in Mpc/h comoving.

    Args:
        Rlow (float): Inner radii.
        Rhigh (float): Outer radii.
        R (array like): Radii of the profile.
        prof (array like): Projected profile.

    Returns:
        float: Average profile in the radial bin, or annulus.

    """
    return average_profile_in_bins([Rlow, Rhigh], R, prof)
