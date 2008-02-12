#
# VMD script to calculate charge distributions
# outside non-spherical objects.
#
#  M. Lund, Prague 2008
#
# METHODS:
#  sphere:  Cummulative charge sum in spheres
#           around the geometric center
#  surface: Cummulative charge sum in the non-spherical
#           volume outside the surface of the selection
#  shell:   Charge sum in non-spherical slices outside the
#           surface of central selection

from atomsel import *
from Molecule import *

# Select first molecule
m = moleculeList()[0]
n = m.numFrames()
n = 11000
method = "surface"
center = "protein"
volume = 0

beg=0     # Starting radius
end=30    # Ending radius
dr=0.5    # Resolution [AA]
frame = 0 # Starting frame

qavg=[0.]*int(end/dr+2)
vavg=[0.]*int(end/dr+2)

# Loop over frames
while frame<n:
  print "Frame: ",frame+1,"/",n
  # Loop over distances
  for b in range(beg, int(end/dr)):

    if method=="sphere":
      sel = atomsel(center, 0, frame)
      cm = sel.center()
      selstr = "sqr(x-" + str(cm[0]) + ")" \
          + " + sqr(y-" + str(cm[1]) + ")" \
          + " + sqr(z-" + str(cm[2]) + ")" \
          + " < sqr(" + str(b*dr) + ")"
      sel = atomsel(selstr,0,frame)

    elif method=="surface":
      selstr = "exwithin "+str(b*dr)+" of "+center
      sel = atomsel(selstr,0,frame)

    elif method=="shell":
      selstr = "(exwithin "+str(b*dr+dr)+" of "+center+") "\
          + "and (not exwithin "+str(b*dr)+" of "+center+")"
      sel = atomsel(selstr,0,frame)

    if volume==1:
      selstr = "name OW and ("+selstr+")"
      sel = atomsel(selstr,0,frame)
      vol = len(sel)/(55*6.022e23)*1e30 # Aangstrom^3
      vavg+=vol

    # Loop over charges
    qvec = sel.get("charge")
    qtot = 0
    for q in qvec:
      qtot+=q
    qavg[b]+=qtot
  frame+=1

# Save average charge as a function of r
f=open("radialcharge-"+method+".out", "w")
for b in range(beg, end/dr):
  f.write( str(b*dr)+" "+str(qavg[b]/frame)+"\n" )
f.close()

