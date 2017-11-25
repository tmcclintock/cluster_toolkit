import pytest
from cluster_toolkit import deltasigma as ds
from os.path import dirname, join
import numpy as np
import numpy.testing as npt

#Need some test data to use first
M = 1e14
c = 5
Om = 0.3
here = dirname(__file__)
Rxi = np.loadtxt(here+"/data_for_testing/r3d.txt")
xihm = np.loadtxt(here+"/data_for_testing/xi_hm.txt")
R = np.logspace(-1, 2, num = 1000)

def test_Sigma():
    arrout = ds.Sigma_at_R(R, Rxi, xihm, M, c, Om)
    assert len(arrout) == len(R)
    for i in range(len(R)):
        npt.assert_almost_equal(arrout[i], ds.Sigma_at_R(R[i], Rxi, xihm, M, c, Om))
