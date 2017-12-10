import pytest
from cluster_toolkit import concentration
from os.path import dirname, join
import numpy as np
import numpy.testing as npt

#Halo properties that are inputs
rhomconst = 2.77533742639e+11 #units are SM h^2/Mpc^3
Mass = 1e14 #Msun/h
Omega_m = 0.3 #arbitrary
datapath = "./data_for_testing/"
k = np.loadtxt(join(dirname(__file__),datapath+"klin.txt")) #h/Mpc; wavenumber
p = np.loadtxt(join(dirname(__file__),datapath+"plin.txt")) #[Mpc/h]^3 linear power spectrum
Marr = np.array([1e13, 1e14, 1e15]) #Msun/h

def test_exceptions():
    with pytest.raises(Exception):
        concentration.concentration_at_M(Mass, k, p, Omega_m, Mass_type="blah")

        
if __name__ == "__main__":
    print concentration.concentration_at_M(Mass, k, p, Omega_m, Mass_type="crit")
    print concentration.concentration_at_M(Mass, k, p, Omega_m, Mass_type="mean")
