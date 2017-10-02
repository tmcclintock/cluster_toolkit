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
rs = 0.2 #Mpc/h comoving
alpha = 0.18 #default value in Baxter
om = 0.3
xi_nfw   = clusterwl.xi.xi_nfw_at_R(R, M, c, om)
xi_ein   = clusterwl.xi.xi_einasto_at_R(R, M, rs, alpha, om)
#print xi_nfw[:5]/xi_ein[:5]

xi_mm    = clusterwl.xi.xi_mm_at_R(R, k, P)
bias = clusterwl.bias.bias_at_M(M, klin, Plin, om)
xi_2halo = clusterwl.xi.xi_2halo(bias, xi_mm)
xi_hm    = clusterwl.xi.xi_hm(xi_nfw, xi_2halo)


def plot_xi():
    plt.loglog(R, xi_nfw, label="nfw")
    plt.loglog(R, xi_ein, label="ein")
    #plt.loglog(R, xi_mm)
    plt.loglog(R, xi_2halo, label="2h")
    #plt.loglog(R, xi_hm, ls='--')
    plt.legend(loc=0)
    plt.show()

if __name__ == "__main__":
    #plot_bias()
    plot_xi()
    #plot_Sigma()
    #plot_DeltaSigma()
