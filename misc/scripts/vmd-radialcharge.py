#
# VMD script to calculate charge distributions
# outside non-spherical objects.
#
# M. Lund, Prague 2008
#

from atomsel import *
from Molecule import *

# Select first molecule
m = moleculeList()[0]
n = m.numFrames()
method = "shell" # [ sphere,surface,shell ]
center = "protein"

beg=0     # Starting radius
end=40    # Ending radius
dr=1.0    # Resolution [AA]
frame = 0 # Starting frame

qavg=[0]*int(end/dr)

# Loop over frames
while frame<n:
  print "Processing frame: ", frame+1, "/", n
  # Loop over distances
  for b in range(beg, int(end/dr)):

    # Sum spherical shells around central selection
    if method=="sphere":
      sel = atomsel(center, 0, frame)
      cm = sel.center()
      sel = atomsel( "sqr(x-" + str(cm[0]) + ")" \
          + " + sqr(y-" + str(cm[1]) + ")" \
          + " + sqr(z-" + str(cm[2]) + ")" \
          + " < sqr(" + str(b*dr) + ")", 0, frame)

    # Sum in non-spherical volume outside central selection
    elif method=="surface":
      sel = atomsel("exwithin "+str(b*dr)+" of "+center,0,frame)

    # Sum in non-spherical SLICES outside central selection
    elif method=="shell":
      sel = atomsel("(exwithin "+str(b*dr+dr)+" of "+center+") "\
          + "and (not exwithin "+str(b*dr)+" of "+center+")"\
          ,0,frame)

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

