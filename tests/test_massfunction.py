import pytest
import numpy as np
from classy import Class
from cluster_toolkit import massfunction as mf
import numpy.testing as npt

cos = {"Omega_m":0.3, "Omega_b":0.05, "h":0.7, "A_s":2.19e-9, "sigma8":0.8, "n_s":0.96, "w0":-1.0, "N_eff":3.0146}
h = cos['h']
Omega_m = cos['Omega_m']
params = {
    'output': 'mPk', #linear only
    'h': cos['h'],
    'A_s': cos['A_s'],
    #'sigma8': cos['sigma8'],
    'n_s': cos['n_s'],
    'w0_fld': cos['w0'],
    'wa_fld': 0.0,
    'Omega_b': cos['Omega_b'],
    'Omega_cdm': cos['Omega_m'] - cos['Omega_b'],
    'Omega_Lambda': 1.- cos['Omega_m'],
    'N_eff':cos['N_eff'],
    'P_k_max_1/Mpc':10.,
    'z_max_pk':10.
}
cosmo = Class()
cosmo.set(params)
cosmo.compute()
z = 0
k = np.logspace(-5, 1, num=1000) #Mpc^-1
p = np.array([cosmo.pk_lin(ki, z) for ki in k])*h**3 #Mpc^3/h^3

M = np.logspace(12, 16, num=20)

def test_dndM():
    n = mf.dndM_at_M(M, k/h, p, Omega_m)*M
    k2 = k/h
    n2 = np.array([mf.dndM_at_M(Mi, k2, p, Omega_m) for Mi in M])*M
    npt.assert_array_equal(n, n2)

def test_dndM_M():
    n = mf.dndM_at_M(M, k/h, p, Omega_m)*M
    npt.assert_array_less(n[1:], n[:-1])

def test_dndM_z():
    #In the high mass end high z is always steeper
    M = np.logspace(14, 16, num=10) 
    nz0 = mf.dndM_at_M(M, k/h, p, Omega_m)*M
    z=1.
    p1 = np.array([cosmo.pk_lin(ki, z) for ki in k])*h**3 #Mpc^3/h^3
    nz1 = mf.dndM_at_M(M, k/h, p1, Omega_m)*M
    npt.assert_array_less(nz1, nz0)

if __name__ == "__main__":
    test_dndM()
