import pytest
from cluster_toolkit import boostfactors as bfs
from os.path import dirname, join
import numpy as np
import numpy.testing as npt

R = np.logspace(-1,1.5,num=100)
B0 = .19
Rs = 1.0
alpha = -1.0

def test_nfw():
    npt.assert_array_less(1.0, bfs.boost_nfw_at_R(R, B0, Rs))

def test_powerlaw():
    npt.assert_array_less(1.0, bfs.boost_powerlaw_at_R(R, B0, Rs, alpha))

def test_single_vs_array():
    arrout = np.array([bfs.boost_nfw_at_R(Ri, B0, Rs) for Ri in R])
    npt.assert_array_equal(arrout, bfs.boost_nfw_at_R(R, B0, Rs))
    arrout = np.array([bfs.boost_powerlaw_at_R(Ri, B0, Rs, alpha) for Ri in R])
    npt.assert_array_equal(arrout, bfs.boost_powerlaw_at_R(R, B0, Rs, alpha))
