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
dpath = dirname(__file__)+"/data_for_testing/"
R = np.loadtxt(dpath+"Rds.txt")
Sigma = np.loadtxt(dpath+"Sigma.txt")
DeltaSigma = np.loadtxt(dpath+"DeltaSigma.txt")
Rm = np.logspace(np.log10(min(R)), np.log10(max(R)), num=100)

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
