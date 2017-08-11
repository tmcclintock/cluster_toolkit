import clusterwl
import numpy as np
import matplotlib.pyplot as plt
import pyccl as ccl
import time

h = 0.7
om = 0.3
A_s = 2.1e-9 #Power spectrum amplitude
n_s = 0.96 #Power spectrum index
Omega_b = 0.05 #NOTE: use Omega, not omega, since they mean different things
Omega_m = 0.3
Omega_cdm = Omega_m - Omega_b
NM = 1000

start = time.time()
p = ccl.Parameters(Omega_c=Omega_cdm, Omega_b=Omega_b, h=h, A_s=A_s, n_s=n_s)
cosmo = ccl.Cosmology(p)
end = time.time()
print "ccl compute:", end-start
k = np.logspace(-5, 2, base=10, num=1000) #h/Mpc
z = 0.0
a = 1./(1+z)
start = time.time()
P = ccl.linear_matter_power(cosmo, k, a)*h**3
end = time.time()
print "ccl P(k):", end-start


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
