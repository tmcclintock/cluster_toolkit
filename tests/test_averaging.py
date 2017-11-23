import pytest
from cluster_toolkit import averaging
from os.path import dirname, join
import numpy as np
import numpy.testing as npt

R = np.logspace(-1, 1, num=100)
def profile_func(R):
    return R**2
prof = profile_func(R) #Imaginary profile is an R^2 power law
#Any bin should have (2/(R_2^2 - R_1^2)) * (R_2^4 - R_1^4)/4

def true_ave(Rlow, Rhigh):
    return (Rhigh**4 - Rlow**4)/(Rhigh**2 - Rlow**2)/2.

def test_bin():
    Rlow  = .5
    Rhigh = 1.5
    npt.assert_almost_equal(true_ave(Rlow, Rhigh), averaging.average_profile_in_bin(Rlow, Rhigh, R, prof))

def test_bins():
    Rbins = [[.6, 1.6], [1.6, 1.8]]
    Redges = [.6, 1.6, 1.8]
    aves = np.array([true_ave(Rlo, Rhi) for Rlo,Rhi in Rbins])
    npt.assert_array_almost_equal(aves, averaging.average_profile_in_bins(Redges, R, prof))

#Regression tests below. TODO
