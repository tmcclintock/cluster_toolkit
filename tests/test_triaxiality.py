import pytest
from cluster_toolkit import triaxiality as tr
from os.path import dirname, join
import numpy as np
import numpy.testing as npt

def test_mapped_r():
    zs = np.linspace(0.01, 10)
    R = 5
    r = np.sqrt(zs**2 + R**2)

    #Test that no triaxiality results in the same thing as a spherical halo
    mapped_r = tr.mapped_radii(zs, R, i=0, phi=0, q=1, s=1)
    npt.assert_equal(r, mapped_r)

    #Test that any triaxiality results in larger mapped radii
    mapped_r = tr.mapped_radii(zs, R, i=np.pi/2, phi=np.pi, q=0.6, s=0.6)
    npt.assert_array_less(r, mapped_r)

    #Test that if the inputs are all arrays, then we get
    #an output of the correct dimensions and size
    zs = np.linspace(0.01, 10, 4)
    Rs = zs
    i = np.linspace(0, np.pi/2, 4)
    phi = np.linspace(0, np.pi, 4)
    mapped_r = tr.mapped_radii(zs, Rs, i, phi, q=0.6, s=0.6)
    npt.assert_equal(len(zs)**4, mapped_r.size)

    return
