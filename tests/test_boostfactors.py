import pytest
from cluster_toolkit import boostfactors as bfs
from os.path import dirname, join
import numpy as np
import numpy.testing as npt

R = np.logspace(-1,1.5,num=100)
B0 = .19
Rs = 1.0

def test_nfw():
    npt.assert_array_less(1.0, bfs.boost_nfw_at_R(R, B0, Rs))

def test_powerlaw():
    alpha = -1.0
    npt.assert_array_less(1.0, bfs.boost_powerlaw_at_R(R, B0, Rs, alpha))
