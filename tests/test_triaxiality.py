import pytest
import cluster_toolkit as ct
from cluster_toolkit import deltasigma as ds
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

def FAILING_test_single_halo_NFW():
    R = np.logspace(-1, 2, num = 100)
    M = 1e14
    c = 5
    Omega_m = 0.3

    r = np.logspace(-2, 3, 1000)
    xi = ct.xi.xi_nfw_at_r(r, M, c, Omega_m)
    import matplotlib.pyplot as plt
    plt.loglog(r, xi)
    plt.show()
    arr1 = ds.Sigma_at_R(R, r, xi, M, c, Omega_m)
    arr2 = tr.Ellipsoidal_Sigma_nfw_single_halo(R, M, c, 0, 1, 1, Omega_m)
    plt.loglog(R, arr2)
    plt.loglog(R, arr1)
    plt.yscale("symlog")
    plt.show()

    npt.assert_equal(arr1, arr2)
    return

#test_single_halo_NFW()
