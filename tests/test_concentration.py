import pytest
from cluster_toolkit import concentration
from os.path import dirname, join
import numpy as np
import numpy.testing as npt

#Halo properties that are inputs
h = 0.7
Omega_m = 0.3
Omega_b = 0.05
sigma8 = 0.847533481324
ns = 0.96
As = 2.215e-9
rhomconst = 2.77533742639e+11 #units are SM h^2/Mpc^3
Mass = 1e14 #Msun/h
datapath = "./data_for_testing/"
k = np.loadtxt(join(dirname(__file__),datapath+"klin.txt")) #h/Mpc; wavenumber
p = np.loadtxt(join(dirname(__file__),datapath+"plin.txt")) #[Mpc/h]^3 linear power spectrum
Marr = np.array([1e13, 1e14, 1e15]) #Msun/h

def test_exceptions():
    with pytest.raises(Exception):
        concentration.concentration_at_M(Mass, k, p, ns, Omega_b, Omega_m, h, Mass_type="vir")
        concentration.concentration_at_M(Mass, k, p, ns, Omega_b, Omega_m, h, delta=300)

if __name__ == "__main__":
    print concentration.concentration_at_M(Mass, k, p, ns, Omega_b, Omega_m, h, Mass_type="crit")
    print concentration.concentration_at_M(Mass, k, p, ns, Omega_b, Omega_m, h, Mass_type="mean")
    test_exceptions()
