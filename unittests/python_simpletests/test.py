#!/usr/bin/env python
import sys
from math import sqrt
import pyfaunus  as pan

inp  = pan.inputfile("test.conf")
test = pan.checkValue(inp)

# --- AVERAGE ---
a = pan.average_dbl()
b = pan.average_dbl()
a+=3.0
a+=1.1
b+=1.0
test.check("average_avg()", a.avg())
test.check("average_multiply", a*2.0)
c=a+b
test.check("average_add", c.avg())
test.check("average_add_control", (3+1.1+1)/3)

# --- POINT ---
a = pan.point()
b = pan.point()
a.x = -1.1
b.x = 10.0
test.check("point_len()", a.len())
test.check("point_dist()", a.dist(b))

# --- CELL ---
a = pan.particle()
b = pan.particle()
a.charge=+1
a.x=10.
b.charge=-1
b.z=50.
con = pan.cell(100.)
pot = pan.pot_coulomb(inp)
test.check("cell_getvolume()", con.getvolume())
test.check("cell_dist()", con.dist(a,b) )
test.check("cell_dist()_control1", a.dist(b) )
test.check("cell_dist()_control2", sqrt(pot.sqdist(a,b)) )
con.push_back(a)
test.check("cell_charge()", con.charge() )

print test.report()
sys.exit( test.returnCode() )

