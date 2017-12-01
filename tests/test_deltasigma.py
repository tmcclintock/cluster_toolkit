import pytest
from cluster_toolkit import deltasigma as ds
from os.path import dirname, join
import numpy as np
import numpy.testing as npt

#Need some test data to use first
M = 1e14
c = 5
Om = 0.3
try:
    here = dirname(__file__)
    Rxi = np.loadtxt(here+"/data_for_testing/r3d.txt")
    xihm = np.loadtxt(here+"/data_for_testing/xi_hm.txt")
except IOError:
    Rxi = np.loadtxt("./data_for_testing/r3d.txt")
    xihm = np.loadtxt("./data_for_testing/xi_hm.txt")
R = np.logspace(-1, 2, num = 1000)

def test_Sigma_nfw():
    arrout = ds.Sigma_nfw_at_R(R, M, c, Om)
    assert len(arrout) == len(R)
    for i in range(len(R)):
        npt.assert_equal(arrout[i], ds.Sigma_nfw_at_R(R[i], M, c, Om))

def test_Sigma():
    arrout = ds.Sigma_at_R(R, Rxi, xihm, M, c, Om)
    assert len(arrout) == len(R)
    for i in range(len(R)):
        npt.assert_equal(arrout[i], ds.Sigma_at_R(R[i], Rxi, xihm, M, c, Om))

def test_DeltaSigma():
    #Try with nfw
    snfw = ds.Sigma_nfw_at_R(R, M, c, Om)
    arrout = ds.DeltaSigma_at_R(R, R, snfw, M, c, Om)
    assert len(arrout) == len(R)
    for i in range(len(R)):
        npt.assert_equal(arrout[i], ds.DeltaSigma_at_R(R[i], R, snfw, M, c, Om))
    #Now try with general Sigma(R)
    sig = ds.Sigma_at_R(R, Rxi, xihm, M, c, Om)
    arrout = ds.DeltaSigma_at_R(R, R, sig, M, c, Om)
    assert len(arrout) == len(R)
    for i in range(len(R)):
        npt.assert_equal(arrout[i], ds.DeltaSigma_at_R(R[i], R, sig, M, c, Om))

def test_Sigma_and_xi():
    #Test that changes to xi_hm below a cut only affect Sigma below that cut
    arrout = ds.Sigma_at_R(R, Rxi, xihm, M, c, Om)
    norm = arrout.copy()
    xi2 = xihm.copy()
    cut = 0.2
    inds = np.where(Rxi < cut)[0]
    xi2[inds] *= 2
    arrout2 = ds.Sigma_at_R(R, Rxi, xi2, M, c, Om)
    arrout /= norm
    arrout2 /= norm
    inds = np.where(R < cut)[0]
    npt.assert_array_less(arrout[inds], arrout2[inds])
    inds = np.where(R >= cut)[0]
    #Splines mean that we won't get perfectly equal
    npt.assert_array_almost_equal(arrout[inds], arrout2[inds], decimal=2)
