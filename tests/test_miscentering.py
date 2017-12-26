import pytest
from cluster_toolkit import miscentering as mis
from os.path import dirname, join
import numpy as np
import numpy.testing as npt

#Need some test data to use first
M = 1e14
c = 5
Om = 0.3
Rmis = 1.0
dpath = "./data_for_testing/"
R = np.loadtxt(join(dirname(__file__),dpath+"Rds.txt"))
Sigma = np.loadtxt(join(dirname(__file__),dpath+"Sigma.txt"))
DeltaSigma = np.loadtxt(join(dirname(__file__),dpath+"DeltaSigma.txt"))
Rm = np.logspace(np.log10(min(R)), np.log10(max(R)), num=100)

def test_errors():
    with pytest.raises(Exception):
        mis.Sigma_mis_single_at_R(np.min(R)*0.9, R, Sigma, M, c, Om, Rmis)
        mis.Sigma_mis_single_at_R(np.max(R)*1.1, R, Sigma, M, c, Om, Rmis)

def test_Sigma_mis():
    arrout = mis.Sigma_mis_at_R(Rm, R, Sigma, M, c, Om, Rmis)
    assert len(arrout) == len(Rm)
    for i in range(len(Rm)):
        npt.assert_equal(arrout[i], mis.Sigma_mis_at_R(Rm[i], R, Sigma, M, c, Om, Rmis))

def test_Single():
    arrout = mis.Sigma_mis_single_at_R(Rm, R, Sigma, M, c, Om, Rmis)
    assert len(arrout) == len(Rm)
    for i in range(len(Rm)):
        npt.assert_equal(arrout[i], mis.Sigma_mis_single_at_R(Rm[i], R, Sigma, M, c, Om, Rmis))

def test_DSmis():
    Smis = mis.Sigma_mis_at_R(Rm, R, Sigma, M, c, Om, Rmis)
    arrout = mis.DeltaSigma_mis_at_R(Rm, Rm, Smis)
    assert len(arrout) == len(Rm)
    for i in range(len(Rm)):
        npt.assert_equal(arrout[i], mis.DeltaSigma_mis_at_R(Rm[i], Rm, Smis))
    Smis = mis.Sigma_mis_single_at_R(Rm, R, Sigma, M, c, Om, Rmis)
    arrout = mis.DeltaSigma_mis_at_R(Rm, Rm, Smis)
    assert len(arrout) == len(Rm)
    for i in range(len(Rm)):
        npt.assert_equal(arrout[i], mis.DeltaSigma_mis_at_R(Rm[i], Rm, Smis))

def test_nomis():
    #Test what happens when Rmis=0
    Rmis = 0.0001
    lo, hi = 0, -1 #We expect the end points to be off
    dec = 3 #To this decimal place
    ones = np.ones_like(Sigma)
    Smis = mis.Sigma_mis_at_R(R, R, Sigma, M, c, Om, Rmis)/Sigma
    #npt.assert_array_almost_equal(Sigma[lo:hi], Smis[lo:hi], decimal=dec)
    npt.assert_array_almost_equal(ones[lo:hi], Smis[lo:hi], decimal=dec)
    Smis = mis.Sigma_mis_single_at_R(R, R, Sigma, M, c, Om, Rmis)/Sigma
    #npt.assert_array_almost_equal(Sigma[lo:hi], Smis[lo:hi], decimal=dec)
    npt.assert_array_almost_equal(ones[lo:hi], Smis[lo:hi], decimal=dec)
