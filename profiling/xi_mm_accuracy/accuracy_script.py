import pytest
from cluster_toolkit import xi
import numpy as np
import matplotlib.pyplot as plt

R = np.logspace(-1, 3, 100)
k = np.loadtxt("k.txt")
p = np.loadtxt("p.txt")

for i in range(-5,5,1):
    ximm = xi.xi_mm_at_R(R, k, p, step=0.005+i*0.001)
    plt.loglog(R, ximm, label='%d'%i)
plt.legend()
plt.show()


for i in range(195, 205, 1):
    ximm = xi.xi_mm_at_R(R, k, p, N=i)
    plt.loglog(R, ximm, label='%d'%i)
plt.legend()
plt.show()
