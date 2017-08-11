import clusterwl
import numpy as np
import matplotlib.pyplot as plt

k = np.loadtxt("test_data/knl.txt")
P = np.loadtxt("test_data/pnl.txt")
klin = np.loadtxt("test_data/klin.txt")
Plin = np.loadtxt("test_data/plin.txt")

NR = 1000
R = np.logspace(-2, 3, NR, base=10) #Xi_hm MUST be evaluated to higher than BAO
M = 1e14 
c = 5
om = 0.3
xi_nfw = np.zeros_like(R)
xi_mm = np.zeros_like(R)
xi_2halo = np.zeros_like(R)
xi_hm = np.zeros_like(R)

Marr = np.logspace(12, 16, NR, base=10)
biases = np.zeros_like(Marr)
clusterwl.bias.calc_bias_at_M(Marr, klin, Plin, om, biases, 500)

bias = clusterwl.bias.bias_at_M(M, klin, Plin, om)
clusterwl.xi.calc_xi_nfw(R, M, c, om, xi_nfw)
clusterwl.xi.calc_xi_mm(R, k, P, xi_mm)
clusterwl.xi.calc_xi_2halo(bias, xi_mm, xi_2halo)
clusterwl.xi.calc_xi_hm(xi_nfw, xi_2halo, xi_hm)
"""
plt.loglog(R, xi_nfw)
plt.loglog(R, xi_mm)
plt.loglog(R, xi_2halo)
plt.loglog(R, xi_hm, ls='--')
plt.show()
plt.clf()
"""

Rp = np.logspace(-2, 2.4, NR, base=10)
Sigma  = np.zeros_like(Rp)
Sigma_nfw = np.zeros_like(Rp)
DeltaSigma = np.zeros_like(Rp)
DeltaSigma_nfw = np.zeros_like(Rp)

clusterwl.deltasigma.calc_Sigma_at_R(Rp, R, xi_hm, M, c, om, Sigma)
clusterwl.deltasigma.calc_Sigma_nfw_at_R(Rp, M, c, om, Sigma_nfw)
clusterwl.deltasigma.calc_DeltaSigma_at_R(Rp, Rp, Sigma, M, c, om, DeltaSigma)
clusterwl.deltasigma.calc_DeltaSigma_at_R(Rp, Rp, Sigma_nfw, M, c, om, DeltaSigma_nfw)

Rmis = 0.25 #Mpc/h
Sigma_mis  = np.zeros_like(Rp)
Sigma_single  = np.zeros_like(Rp)
DeltaSigma_mis = np.zeros_like(Rp)
DeltaSigma_single = np.zeros_like(Rp)

clusterwl.miscentering.calc_Sigma_mis_at_R(Rp, Rp, Sigma, M, c, om, Rmis, Sigma_mis)
clusterwl.miscentering.calc_Sigma_mis_single_at_R(Rp, Rp, Sigma, M, c, om, Rmis, Sigma_single)

clusterwl.miscentering.calc_DeltaSigma_mis_at_R(Rp, Rp, Sigma_mis, DeltaSigma_mis)
clusterwl.miscentering.calc_DeltaSigma_mis_at_R(Rp, Rp, Sigma_single, DeltaSigma_single)

#plt.loglog(Rp, Sigma, ls="--", label=r"$\Sigma$")
#plt.loglog(Rp, Sigma_mis, ls=":", label=r"$\Sigma_{mis}$")
#plt.loglog(Rp, Sigma_single, ls=":", label=r"$\Sigma(R_{mis})$")
#plt.loglog(Rp, Sigma_nfw, label=r"$\Sigma_{nfw}$")
plt.loglog(Rp, DeltaSigma, label=r"$\Delta\Sigma$")
plt.loglog(Rp, DeltaSigma_mis, label=r"$\Delta\Sigma_{mis}$")
plt.loglog(Rp, DeltaSigma_single, label=r"$\Delta\Sigma(R_{mis})$")
plt.loglog(Rp, DeltaSigma_nfw, label=r"$\Delta\Sigma_{nfw}$")
plt.legend()
plt.subplots_adjust(wspace=0.01)
plt.show()
