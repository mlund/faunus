#!/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp

r_mc, g_mc = np.loadtxt('rdf_mc.dat', unpack=True, usecols=(0,1))
r_ld, g_ld = np.loadtxt('rdf_ld.dat', unpack=True, usecols=(0,1))

stats = ks_2samp(r_mc, r_ld)
if stats[1] > 0.99:
   print('p-value = {} (no significant difference)'.format(stats[1]))
else:
   print('p-value = {} (significant difference)'.format(stats[1]))

plt.plot(r_mc, g_mc, label='Monte Carlo', linewidth=10, alpha=0.5)
plt.plot(r_ld, g_ld, label='Langevin', linewidth=2)
plt.legend(loc=0, frameon=False, fontsize=10)
plt.xlabel('Distance r (Å)', fontsize=14)
plt.ylabel('g(r)', fontsize=14)
plt.show()
