#!/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import os

mc = np.loadtxt('improper_dihedral_mc_rdf.dat', unpack=True, usecols=(0,1))
ld = np.loadtxt('improper_dihedral_md_rdf.dat', unpack=True, usecols=(0,1))

plt.plot(mc[0], mc[1], label='Monte Carlo', linewidth=10, alpha=0.3)
plt.plot(ld[0], ld[1], label='Langevin', linewidth=2)

if os.path.isfile('rdf.dat'):
    test = np.loadtxt('rdf.dat', unpack=True, usecols=(0,1))
    plt.plot(test[0], test[1], label='Test', linewidth=2)

plt.xlim(0, 7)
plt.legend(loc=0, frameon=False, fontsize=10)
plt.xlabel('Distance, $r$ (Å)', fontsize=14)
plt.ylabel('$g(r)$', fontsize=14)
plt.savefig('improper_dihedral_plot.png')
