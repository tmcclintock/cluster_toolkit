import cluster_toolkit
from ctypes import c_double, c_int, POINTER

def _dcast(x):
    return cluster_toolkit._ffi.cast('double*', x.ctypes.data)

def average_profile_in_bin(Rlow, Rhigh, R, prof):
    """Calculate the average of some projected profile in a radial bin in Mpc/h comoving

    Args:
    Rlow (float): Inner radii
    Rhigh (float): Outer radii
    R (array like): Radii of the profile
    prof (array like): Projected profile

    Returns:
    ave_prof (float): Average profile in the radial bin, or annulus

    """
    return cluster_toolkit._lib.average_profile_in_bin(Rlow, Rhigh, dcat(R), len(R), _dcast(prof))

def average_profile_in_bins(Redges, R, prof):
    """Calculate the average of some projected profile in a radial bins in Mpc/h comoving

    Args:
    Redges (array like): Array of radial bin edges
    R (array like): Radii of the profile
    prof (array like): Projected profile
    ave_prof (float): Average profile in the radial bins, populated with the result

    """
    ave_prof = np.zeros(len(Redges)-1)
    cluster_toolkit._lib.average_profile_in_bins(_dcast(Redges), len(Redges), _dcast(R), len(R), _dcast(prof), _dcast(ave_prof))
    return ave_prof
