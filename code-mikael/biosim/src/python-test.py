#
# swig -python -c++ point.h (see header modification)
# gcc ...something
gcc -c point.c example_wrap.c -I/usr/include/python2.2 -I/usr/lib/python2.2
# ld ...a_lot_of_mac_specific_flags
#
#

import point
import copy

p=[]

p.append( point.point() )
p.append( point.point() )

p[0].x=2.
p[1].y=5.

print p[0].x, p[1].x, p[0].dist(p[1])
