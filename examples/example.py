import clusterwl
import numpy as np
import matplotlib.pyplot as plt

k = np.loadtxt("test_data/knl.txt")
P = np.loadtxt("test_data/pnl.txt")
klin = np.loadtxt("test_data/klin.txt")
Plin = np.loadtxt("test_data/plin.txt")

NR = 1000
R = np.logspace(-2, 2.1, NR, base=10)
M = 1e14 
c = 5
om = 0.3
xi_nfw = np.zeros_like(R)
xi_mm = np.zeros_like(R)
xi_2halo = np.zeros_like(R)
xi_hm = np.zeros_like(R)

Marr = np.logspace(12, 16, NR, base=10)
biases = np.zeros_like(Marr)
clusterwl.bias.calc_bias_at_M(Marr, k, P, om, biases, 500)

print clusterwl.xi.xi_nfw_at_R(R[0], M, c, om)
print clusterwl.xi.xi_mm_at_R(R[0], k, P)
bias = clusterwl.bias.bias_at_M(M, klin, Plin, om)
print "b(%.2e) = %.3f"%(M, bias)
clusterwl.xi.calc_xi_nfw(R, M, c, om, xi_nfw)
clusterwl.xi.calc_xi_mm(R, k, P, xi_mm)
clusterwl.xi.calc_xi_2halo(bias, xi_mm, xi_2halo)
clusterwl.xi.calc_xi_hm(xi_nfw, xi_2halo, xi_hm)

plt.loglog(R, xi_nfw)
plt.loglog(R, xi_mm)
plt.loglog(R, xi_2halo)
plt.loglog(R, xi_hm, ls='--')
plt.show()
