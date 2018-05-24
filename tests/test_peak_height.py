import pytest
from cluster_toolkit import peak_height as peaks
from os.path import dirname, join
import numpy as np
import numpy.testing as npt

#Halo properties that are inputs
rhomconst = 2.77533742639e+11 #units are SM h^2/Mpc^3
Mass = 1e14 #Msun/h
Omega_m = 0.3 #arbitrary
datapath = "./data_for_testing/"
klin = np.loadtxt(join(dirname(__file__),datapath+"klin.txt"))#h/Mpc; wavenumber
plin = np.loadtxt(join(dirname(__file__),datapath+"plin.txt"))#[Mpc/h]^3 linear power spectrum
Ma = np.array([1e13, 1e14, 1e15]) #Msun/h
Ra = (Ma/(4./3.*np.pi*rhomconst*Omega_m))**(1./3.)


def test_s2_and_nu_functions():
    #Test the mass calls
    s2 = peaks.sigma2_at_M(Mass, klin, plin, Omega_m)
    nu = peaks.nu_at_M(Mass, klin, plin, Omega_m)
    npt.assert_equal(1.686/np.sqrt(s2), nu)
    s2 = peaks.sigma2_at_M(Ma, klin, plin, Omega_m)
    nu = peaks.nu_at_M(Ma, klin, plin, Omega_m)
    npt.assert_array_equal(1.686/np.sqrt(s2), nu)
    #Now test the R calls
    R = 1.0 #Mpc/h; arbitrary
    s2 = peaks.sigma2_at_R(R, klin, plin)
    nu = peaks.nu_at_R(R, klin, plin)
    npt.assert_equal(1.686/np.sqrt(s2), nu)

def test_single_vs_array():
    #First sigma2
    a1 = peaks.sigma2_at_M(Ma, klin, plin, Omega_m)
    a2 = np.array([peaks.sigma2_at_M(Mi, klin, plin, Omega_m) for Mi in Ma])
    npt.assert_array_equal(a1, a2)
    a1 = peaks.sigma2_at_R(Ra, klin, plin)
    a2 = np.array([peaks.sigma2_at_R(Ri, klin, plin) for Ri in Ra])
    npt.assert_array_equal(a1, a2)
