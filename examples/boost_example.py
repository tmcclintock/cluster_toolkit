import numpy as np
import matplotlib.pyplot as plt
import clusterwl

R = np.logspace(-2, 2, num=100)
B0 = 0.15
Rs = 1.0

boost = clusterwl.boostfactors.boost_nfw_at_R(R, B0, Rs)

plt.plot(R, boost)
plt.axhline(1, c='k')
plt.xscale('log')
plt.show()
