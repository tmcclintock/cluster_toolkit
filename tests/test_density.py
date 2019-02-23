import pytest
from cluster_toolkit import density, xi
import numpy as np
import numpy.testing as npt

#Halo properties that are inputs
R = 1.0 #Mpc/h
R_arr = np.array([0.1, 1.0, 10.0]) #Mpc/h
Mass = 1e14 #Msun/h
conc = 5 #concentration; no units, for NFW
Rscale = 1.0 #Mph/h; for Einasto
alpha = 0.19 #typical value; for Einasto
Omega_m = 0.3 #arbitrary
rhocrit = 2.77533742639e+11
rhom = rhocrit*Omega_m

def test_exceptions_rho_nfw_at_r():
    with pytest.raises(TypeError):
        density.rho_nfw_at_r() #No args
        density.rho_nfw_at_r(R, Mass, conc) #Too few args
        density.rho_nfw_at_r(R, Mass, conc, Omega_m, Omega_m) #Too many args
        density.rho_nfw_at_r("a string", Mass, conc, Omega_m, Omega_m) #Wrong type

def test_outputs_rho_nfw_at_r():
    #List vs. numpy.array
    npt.assert_array_equal(density.rho_nfw_at_r(R_arr, Mass, conc, Omega_m), density.rho_nfw_at_r(R_arr.tolist(), Mass, conc, Omega_m))
    #Single value vs numpy.array
    arrout = density.rho_nfw_at_r(R_arr, Mass, conc, Omega_m)
    for i in range(len(R_arr)):
        npt.assert_equal(density.rho_nfw_at_r(R_arr[i], Mass, conc, Omega_m), arrout[i])

def test_exceptions_rho_einasto_at_r():
    with pytest.raises(TypeError):
        density.rho_einasto_at_r() #No args
        density.rho_einasto_at_r(R, Mass, Rscale, alpha) #Too few args
        density.rho_einasto_at_r(R, Mass, Rscale, alpha, Omega_m, Omega_m) #Too many args
        density.rho_einasto_at_r("a string", Mass, Rscale, alpha, Omega_m, Omega_m) #Wrong type

def test_outputs_rho_einasto_at_r():
    #List vs. numpy.array
    arr1 = density.rho_einasto_at_r(R_arr, Mass, Rscale, alpha, Omega_m)
    arr2 = density.rho_einasto_at_r(R_arr.tolist(), Mass, Rscale, alpha, Omega_m)
    npt.assert_array_equal(arr1, arr2)
    #Single value vs numpy.array
    arr1 = density.rho_einasto_at_r(R_arr, Mass, Rscale, alpha, Omega_m)
    arr2 = np.array([density.rho_einasto_at_r(R_arr[i], Mass, Rscale, alpha, Omega_m) for i in range(len(R_arr))])
    npt.assert_array_equal(arr1, arr2)

def test_rho_vs_xi_nfw():
    #Note: testing fractional precision to avoid numerical issues
    arr1 = density.rho_nfw_at_r(R_arr, Mass, conc, Omega_m)
    arr2 = xi.xi_nfw_at_r(R_arr, Mass, conc, Omega_m)
    arr2 = rhom*(1+arr2)
    arr2 = arr2/arr1
    ones = np.ones_like(arr1)
    npt.assert_array_equal(ones, arr2)
    
def test_rho_vs_xi_einasto():
    arr1 = density.rho_einasto_at_r(R_arr, Mass, Rscale, alpha, Omega_m)
    arr2 = xi.xi_einasto_at_r(R_arr, Mass, Rscale, alpha, Omega_m)
    arr2 = rhom*(1+arr2)
    arr2 = arr2/arr1
    ones = np.ones_like(arr1)
    npt.assert_array_almost_equal(ones, arr2)


if __name__ == "__main__":
    #test_outputs_rho_nfw_at_r()
    #test_outputs_rho_einasto_at_r()
    test_outputs_rho_nfw_at_r()
