import pytest
from cluster_toolkit import xi
import numpy as np
import numpy.testing as npt

R = 1.0 #Mpc/h
R_nparray = np.array([0.1, 1.0, 10.0]) #Mpc/h
R_list = R_nparray.tolist() #Mpc/h
Mass = 1e14 #Msun/h
conc = 5 #concentration; no units, for NFW
Rscale = 1.0 #Mph/h; for Einasto
alpha = 0.19 #typical value; for Einasto
Omega_m = 0.3 #arbitrary


def test_exceptions_xi_nfw_at_R():
    with pytest.raises(TypeError):
        xi.xi_nfw_at_R() #No args
        xi.xi_nfw_at_R(R, Mass, conc) #Too few args
        xi.xi_nfw_at_R(R, Mass, conc, Omega_m, Omega_m) #Too many args
        xi.xi_nfw_at_R("a string", Mass, conc, Omega_m, Omega_m) #Wrong type

def test_outputs_xi_nfw_at_R():
    #List vs. numpy.array
    npt.assert_array_equal(xi.xi_nfw_at_R(R_nparray, Mass, conc, Omega_m), xi.xi_nfw_at_R(R_list, Mass, conc, Omega_m))
    #Single value vs numpy.array
    arrout = xi.xi_nfw_at_R(R_nparray, Mass, conc, Omega_m)
    for i in range(len(R_nparray)):
        npt.assert_equal(xi.xi_nfw_at_R(R_nparray[i], Mass, conc, Omega_m), arrout[i])

def test_exceptions_xi_einasto_at_R():
    with pytest.raises(TypeError):
        xi.xi_einasto_at_R() #No args
        xi.xi_einasto_at_R(R, Mass, Rscale, alpha) #Too few args
        xi.xi_einasto_at_R(R, Mass, Rscale, alpha, Omega_m, Omega_m) #Too many args
        xi.xi_einasto_at_R("a string", Mass, Rscale, alpha, Omega_m, Omega_m) #Wrong type

def test_outputs_xi_einasto_at_R():
    #List vs. numpy.array
    npt.assert_array_equal(xi.xi_einasto_at_R(R_nparray, Mass, Rscale, alpha, Omega_m), xi.xi_einasto_at_R(R_list, Mass, Rscale, alpha, Omega_m))
    #Single value vs numpy.array
    arrout = xi.xi_einasto_at_R(R_nparray, Mass, Rscale, alpha, Omega_m)
    for i in range(len(R_nparray)):
        npt.assert_equal(xi.xi_einasto_at_R(R_nparray[i], Mass, Rscale, alpha, Omega_m), arrout[i])
