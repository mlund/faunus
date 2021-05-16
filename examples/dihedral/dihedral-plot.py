#!/bin/env python
import matplotlib.pyplot as plt
import numpy as np

r_mc, g_mc = np.loadtxt('dihedral-mc-rdf.dat', unpack=True, usecols=(0,1))
r_ld, g_ld = np.loadtxt('dihedral-md-rdf.dat', unpack=True, usecols=(0,1))
r_test, g_test = np.loadtxt('rdf.dat', unpack=True, usecols=(0,1))

plt.plot(r_mc, g_mc, label='Monte Carlo', linewidth=10, alpha=0.5)
plt.plot(r_ld, g_ld, label='Langevin', linewidth=2)
plt.plot(r_test, g_test, label='Test', linewidth=2)
plt.legend(loc=0, frameon=False, fontsize=10)
plt.xlabel('Distance r (Å)', fontsize=14)
plt.ylabel('g(r)', fontsize=14)
plt.savefig('dihedral-plot.png')
