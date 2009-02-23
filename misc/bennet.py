#!/sw/bin/python2.3
from Scientific import *
import sys, math, string

try:
  bennetin = sys.argv[1]; bennetout = sys.argv[2]
except:
  print "Usage: ", sys.argv[0], "bennetin bennetout (alpha dalpha, optional)"; sys.exit(1)

# Open files
g = open( bennetin , 'r'); 
f = open( bennetout, 'r'); 

glines = g.readlines()
flines = f.readlines()

gu = []; gp = []
fu = []; fp = []
for line in glines:
    guv, gpv = line.split()
    gu.append(float(guv)); gp.append(float(gpv))
for line in flines:
    fuv, fpv = line.split()
    fu.append(float(fuv)); fp.append(float(fpv))

g.close()
f.close()

# Initiate Bennet
diffgfok= float(1.e-7)
sumg    = float(0.)
sumf    = float(0.)
gint    = float(0.)
gintold = float(0.)
fint    = float(0.)
gintold = float(0.)
itr     = int(1)
mcnt    = int(1)

try:
  alpha=float(sys.argv[3])
  print 'Alpha set to ', alpha, 'initialy.'
except:
  alpha=1
  print 'No initial guess, alpha set to 1.'
try:
  dalpha  = float(sys.argv[4])
except:
  dalpha  = 0.1*alpha
  print 'dalpha set to ',dalpha

for i in range(len(gp)):
    sumg += gp[i]
    gint += gp[i]/(1. + alpha*math.exp(gu[i]))
for i in range(len(fp)):
    sumf += fp[i]
    fint += fp[i]/(1. + math.exp(fu[i])/alpha)

print 'Sumg  ', sumg ,' Sumf = ', sumf 
print 'Initial gint and fint = ', gint,' ',fint
# Iterate Bennet
print '##############'
while itr%100!=0:
  oldalpha=alpha; oldgint=gint; oldfint=fint 
  alpha+=dalpha
  gint=float(0.); fint=float(0.)
  for i in range(len(gp)):
      gint += gp[i]/(1. + alpha*math.exp(gu[i]))
  for i in range(len(fp)):
      fint += fp[i]/(1. + math.exp(fu[i])/alpha)
  diffgf=abs(fint-gint)
  if itr%10==0:
    print '#Macro ', mcnt
    print 'diffgf = ',diffgf,'  and alpha = ',alpha
    mcnt+=1
  
  itr+=1
  if diffgf<diffgfok:
     print '#############'
     print 'Good convergence!!'
     print 'alpha is equal to = ', alpha, ' or free energy difference (f->g) = ', math.log(alpha), sys.exit(1)
  if fint > gint:
     if oldfint < oldgint:
        dalpha*=-0.5
  elif oldfint > oldgint:
       dalpha*=-0.5

#Iteration not complete
print '############'
print 'Slow convergence, alpha = ', alpha, ', dalpha = ',dalpha, ' free energy difference = ', math.log(alpha)
print 'Last value of gint and fint =', gint, ' ', fint

gfile = open('gben.dat', 'w')
ffile = open('fben.dat', 'w')

for i in range(len(fu)):
    ffile.write('%g %12.5e\n' % (fu[i], fp[i]))
for i in range(len(gu)):
    gfile.write('%g %12.5e\n' % (gu[i], gp[i]))
ffile.close()
gfile.close()

#Energy distribution files should be formated as 

#umin p(umin)
#  ..
#  ..
#umax p(umax
