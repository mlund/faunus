#
# virial.py
#
# This script loads g(r) from disk, converts it to
# free energy, w(r), then shifts it so that the tail
# best fits a salt screened Debye-Huckel potential.
# The PMF is integrated to give the 2nd virial
# coefficient, B2.
#

import numpy as np
import os.path
from sys import exit
from math import sqrt, log, pi, exp, sinh
from scipy.optimize import curve_fit

infile      = "rdf.dat"           # input g(r) arbitrarily normalized
outfile     = "wofr.dat"          # output pmf
rmin        = 40.                 # fitting interval (min)
rmax        = 90.                 # - / / - (max)
lB          = 7.                  # bjerrum length
Q           = [  7.1,  7.1 ]      # charges of the g(r) particles
Mw          = [ 6.52e3, 6.52e3 ]  # molecular weight, g/mol
fit_debye   = True                # true if Debye length should be used to fit data
fit_radius  = True                # true if the radius should be used to fit data
guess_D     = 13.                 # debye length starting point
guess_shift = -2.0                # PMF shift starting point
guess_radius= 9.                  # distance of closest approach (angstrom)
guess       = [ guess_D, guess_shift, guess_radius ]

def f(r,D,shift,a):
  global Q,lB

  if (fit_debye==False):
    D=guess_D

  if (fit_radius==False):
    a=guess_radius

  s=sinh( a/D ) / ( a/D  ) 
  return lB * s**2 * Q[0] * Q[1] / r * np.exp(-r/D) - shift

if (os.path.isfile(infile)==False):
  print "Error: Input rdf file", infile, "is empty or could not be found."
  exit(1)

#
# LOAD AND FIT DATA
#
r, g = np.loadtxt(infile, usecols=(0,1), unpack=True)
w = -np.log(g)                               # convert to free energy
m = (r>=rmin) & (r<=rmax)                    # slice out x range (mask index array)
popt, pcov = curve_fit(f, r[m], w[m], guess) # fit data
debye = popt[0]                              # get fitted values
shift = popt[1]                              # -//-
rfit = r[m]                                  # r array in fitting interval
wfit = f(rfit, *popt) + shift                # w array -//-
w    = w + shift                             # shift PMF
print "Loaded g(r) file      = ", infile
print "Saved w(r) file       = ", outfile
print "Particle charges      = ", Q
print "Particle weights      = ", Mw, "g/mol"
print "Fit range [rmin,rmax] = ", rmin, rmax
print "Fitted Debye length   = ", debye, "A"
print "Fitted ionic strength = ", (3.04/debye)**2*1000., "mM"
print "Fitted shift          = ", shift, "kT"
print "Fitted radius         = ", popt[2], "A"

#
# 2ND VIRIAL COEFFICIENT
#
m      = (r<rmin)                    # slice out all data points below rmin 
r_dat  = r[m]
w_dat  = w[m]
inte   = -2*pi*( np.exp( -w_dat ) - 1 ) * r_dat**2
b2_dat = np.trapz(inte, r_dat)       # integrate data from contact -> rmin

infty=500.                           # assume pmf is ZERO after this point
r_dh  = np.arange( rmin, infty, 0.5 )# generate r array for rmin -> "infinity"
w_dh  = f(r_dh,debye,0,popt[2])              # generate w array -//-
inte  = -2*pi*(np.exp(-w_dh)- 1)*r_dh**2
b2_dh = np.trapz(inte, r_dh )        # integrate using Debye-Huckel
b2_hc = 2 * pi / 3 *r[0]**3          # zero -> contact assuming hard spheres
b2_tot = b2_hc+b2_dat+b2_dh          # total B2 (A**3)

print "Virial coefficient (cubic angstrom):"
print "  Hard sphere  [%5d:%5d] = " % ( 0, r[0] ), b2_hc
print "  Loaded data  [%5d:%5d] = " % ( r[0], rmin ), b2_dat 
print "  Debye-Huckel [%5d:%5d] = " % ( rmin, infty ), b2_dh
print "  TOTAL        [%5d:%5d] = " % ( 0, infty ), b2_tot, "A^3"
print "                             = ", b2_tot*6.022e23*1e-23/Mw[0]/Mw[1], "ml*mol/g^2"
print "  Reduced, B2/B2_HS          = ", b2_tot / b2_hc

#
# SAVE FINAL PMF TO DISK
#
r_final = np.concatenate( (r_dat, r_dh) )  # assemble data and DH r arrays
w_final = np.concatenate( (w_dat, w_dh) )  #                - //- w arrays
np.savetxt(outfile, np.transpose( (r_final,w_final) ) )

#
# PLOT POTENTIAL OF MEAN FORCE ?
# (requires matplotlib)
#
ans = raw_input("\nPlot the pmf (y/N) ? ")
if (ans!="y"):
  exit(0)

import matplotlib.pyplot as plt

plt.plot( r, w, "ro" )
plt.plot( rfit, wfit, linewidth=2.0 )
plt.xlabel("$r$", fontsize=18)
plt.ylabel("$w(r)/k_BT$", fontsize=18)
plt.show()

