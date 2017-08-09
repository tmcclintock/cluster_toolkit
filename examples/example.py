import clusterwl
import numpy as np
import matplotlib.pyplot as plt

k = np.loadtxt("test_data/knl.txt")
P = np.loadtxt("test_data/pnl.txt")

NR = 1000
R = np.logspace(-2, 2.1, NR, base=10)
M = 1e14 
c = 5
om = 0.3
delta = 200
xi_nfw = np.zeros_like(R)
xi_mm = np.zeros_like(R)
xi_2halo = np.zeros_like(R)
xi_hm = np.zeros_like(R)

clusterwl.xi.calc_xi_nfw(R, M, c, delta, om, xi_nfw)
clusterwl.xi.calc_xi_mm(R, k, P, xi_mm)
clusterwl.xi.calc_xi_2halo(2, xi_mm, xi_2halo)
clusterwl.xi.calc_xi_hm(xi_nfw, xi_2halo, xi_hm)

plt.loglog(R, xi_nfw)
plt.loglog(R, xi_mm)
plt.loglog(R, xi_2halo)
plt.loglog(R, xi_hm, ls='--')
plt.show()
