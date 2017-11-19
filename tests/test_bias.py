import pytest
from cluster_toolkit import bias
from os.path import dirname, join
import numpy as np
import numpy.testing as npt

#Halo properties that are inputs
Mass = 1e14 #Msun/h
Omega_m = 0.3 #arbitrary
datapath = "./data_for_testing/"
klin = np.loadtxt(join(dirname(__file__),datapath+"klin.txt")) #h/Mpc; wavenumber
plin = np.loadtxt(join(dirname(__file__),datapath+"plin.txt")) #[Mpc/h]^3 linear power spectrum
M_arr = np.array([1e13, 1e14, 1e15]) #Msun/h

def test_exceptions_bias_at_M():
    with pytest.raises(TypeError):
        bias.bias_at_M()
        bias.bias_at_M(Mass, klin, plin)
        bias.bias_at_M(Mass, klin, plin, Omega_m, Omega_m)
        bias.bias_at_M("a string", klin, plin, Omega_m)

def test_outputs_bias_at_M():
    #List vs. numpy.array
    npt.assert_array_equal(bias.bias_at_M(M_arr, klin, plin, Omega_m), bias.bias_at_M(M_arr.tolist(), klin, plin, Omega_m))
    #Single value vs numpy.array
    arrout = bias.bias_at_M(M_arr, klin, plin, Omega_m)
    for i in range(len(M_arr)):
        npt.assert_equal(bias.bias_at_M(M_arr[i], klin, plin, Omega_m), arrout[i])
