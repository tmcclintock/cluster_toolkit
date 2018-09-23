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

def test_colossus():
    h = 0.7
    Omega_m = 0.3
    Omega_b = 0.05
    sigma8 = 0.847533481324
    ns = 0.96
    As = 2.215e-9
    #Try importing colossus
    try:
        from colossus.halo import concentration as colc
        from colossus.cosmology import cosmology
        cos = {'flat':True,'H0':h*100.,'Om0':Omega_m,'Ob0':Omega_b,'sigma8':sigma8,'ns':ns}
        cosmology.addCosmology('fiducial', cos)
        colcos = cosmology.setCosmology('fiducial')
        k = np.logspace(-5, 2, num=1000)/h #Mpc^-1
        pcol = colcos.matterPowerSpectrum(k) #(Mpc/h)^3
    except ImportError:
        print("colossus not installed, skipping test")

        
if __name__ == "__main__":
    print(concentration.concentration_at_M(Mass, k, p, ns, Omega_b, Omega_m, h, Mass_type="crit"))
    print(concentration.concentration_at_M(Mass, k, p, ns, Omega_b, Omega_m, h, Mass_type="mean"))
    test_exceptions()
    test_colossus()
