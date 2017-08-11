import clusterwl
import numpy as np
import matplotlib.pyplot as plt
import pyccl as ccl
from classy import Class
import time

h = 0.7
om = 0.3
A_s = 2.1e-9 #Power spectrum amplitude
n_s = 0.96 #Power spectrum index
Omega_b = 0.05 #NOTE: use Omega, not omega, since they mean different things
Omega_m = 0.3
Omega_cdm = Omega_m - Omega_b
k = np.loadtxt("test_data/klin.txt")
P = np.loadtxt("test_data/plin.txt")
NM = 1000

params = {
    'output': 'mPk',
    'h':h,
    'A_s':A_s,
    'n_s':n_s,
    'Omega_b':Omega_b,
    'Omega_cdm':Omega_cdm,
    'P_k_max_h/Mpc':1000.,
    'z_max_pk':3.0
    }
cosmo = Class()
cosmo.set(params)
start = time.time()
cosmo.compute()
end = time.time()
print "class compute:", end-start
k = np.logspace(-5, 2, base=10, num=1000) #h/Mpc
z = 0.0
start = time.time()
P = np.array([cosmo.pk_lin(ki*h, z) for ki in k])*h**3 #Mpc^3/h^3
end = time.time()
print "class P(k):", end-start


def plot_massfunction():
    Marr = np.logspace(12, 16, NM, base=10)
    dndM = np.zeros_like(Marr)
    clusterwl.massfunction.calc_dndM_at_M(Marr, k, P, om, dndM)
    plt.loglog(Marr, dndM, label="me")
    plt.legend()
    plt.show()

def plot_N():
    Marr = np.logspace(12, 16, NM, base=10)
    dndM = np.zeros_like(Marr)

    start = time.time()
    clusterwl.massfunction.calc_dndM_at_M(Marr, k, P, om, dndM)
    end = time.time()
    print "dndM calc:", end-start

    volume = 1050.**3 #Mpc^3/h^
    edges = np.logspace(12, 16, 10, base=10)

    start = time.time()
    N = clusterwl.massfunction.N_in_bins(edges, volume, Marr, dndM)
    end = time.time()
    print "N_bins calc:", end-start

    M = (edges[:-1]+edges[1:])/2.
    plt.loglog(M, N)
    plt.show()

if __name__ == "__main__":
    #plot_massfunction()
    plot_N()
