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

def plot_bias():
    NM = 1000
    Marr = np.logspace(12, 16, NM, base=10)
    biases = np.zeros_like(Marr)
    clusterwl.bias.calc_bias_at_M(Marr, klin, Plin, om, biases, 500)
    plt.loglog(Marr, biases)
    plt.show()

bias = clusterwl.bias.bias_at_M(M, klin, Plin, om)
clusterwl.xi.calc_xi_nfw(R, M, c, om, xi_nfw)
clusterwl.xi.calc_xi_mm(R, k, P, xi_mm)
clusterwl.xi.calc_xi_2halo(bias, xi_mm, xi_2halo)
clusterwl.xi.calc_xi_hm(xi_nfw, xi_2halo, xi_hm)

def plot_xi():
    plt.loglog(R, xi_nfw)
    plt.loglog(R, xi_mm)
    plt.loglog(R, xi_2halo)
    plt.loglog(R, xi_hm, ls='--')
    plt.show()


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
Sigma_g2d  = np.zeros_like(Rp)
Sigma_exp  = np.zeros_like(Rp)
Sigma_single  = np.zeros_like(Rp)
DeltaSigma_g2d = np.zeros_like(Rp)
DeltaSigma_exp = np.zeros_like(Rp)
DeltaSigma_single = np.zeros_like(Rp)

clusterwl.miscentering.calc_Sigma_mis_at_R(Rp, Rp, Sigma, M, c, om, Rmis, Sigma_g2d)
clusterwl.miscentering.calc_Sigma_mis_at_R(Rp, Rp, Sigma, M, c, om, Rmis, Sigma_exp, kernel="exponential")
clusterwl.miscentering.calc_Sigma_mis_single_at_R(Rp, Rp, Sigma, M, c, om, Rmis, Sigma_single)
clusterwl.miscentering.calc_DeltaSigma_mis_at_R(Rp, Rp, Sigma_g2d, DeltaSigma_g2d)
clusterwl.miscentering.calc_DeltaSigma_mis_at_R(Rp, Rp, Sigma_exp, DeltaSigma_exp)
clusterwl.miscentering.calc_DeltaSigma_mis_at_R(Rp, Rp, Sigma_single, DeltaSigma_single)

def plot_Sigma():
    plt.loglog(Rp, Sigma, ls="-", label=r"$\Sigma$")
    plt.loglog(Rp, Sigma_g2d, ls=":", label=r"$\Sigma_{mis}^{g2d}$")
    plt.loglog(Rp, Sigma_exp, ls=":", label=r"$\Sigma_{mis}^{exp}$")
    plt.loglog(Rp, Sigma_nfw, ls="--", label=r"$\Sigma_{nfw}$")
    #plt.loglog(Rp, Sigma_single, ls=":", label=r"$\Sigma(R_{mis})$")
    plt.legend()
    plt.show()

def plot_DeltaSigma():
    plt.loglog(Rp, DeltaSigma, ls="-", label=r"$\Delta\Sigma$")
    plt.loglog(Rp, DeltaSigma_g2d, ls=":", label=r"$\Delta\Sigma_{mis}^{g2d}$")
    plt.loglog(Rp, DeltaSigma_exp, ls=":", label=r"$\Delta\Sigma_{mis}^{exp}$")
    plt.loglog(Rp, DeltaSigma_nfw, ls="--", label=r"$\Delta\Sigma_{nfw}$")
    #plt.loglog(Rp, DeltaSigma_single, ls=":", label=r"$\Delta\Sigma(R_{mis})$")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    #plot_Sigma()
    plot_DeltaSigma()
