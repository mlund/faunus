#!/usr/bin/env python

#This program converts trajectory file with spherocylinders(cigars) from fuanus to pdb, psf and vmd.script files.
#Each spherocylinder is represented by two atoms at begining and end and bond between. VMD draw this as spherocylinder
#in cpk draw representation. (if radiuses are set right) Patch on particles is represented by a shifted spherocylinder
#in direction of patch using different color

import os
import sys
import math
import optparse
import commands
import string
import random

def write_pdb(outfilename,box,newdata):
    f = open(outfilename, 'w')
    outstring=""
    for frame in range(len(newdata)):
	newline="CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00 P 1           1\n" % (box[frame][0],box[frame][1],box[frame][2])
	outstring=outstring+newline
	for line in range(len(newdata[frame])):
	    newline="ATOM  %5d  N%1d  PSC F%4d    % 8.3f% 8.3f% 8.3f\n" % (line+1,(line)%4+1 ,(4+line-(line%4))/4,newdata[frame][line][0],newdata[frame][line][1],newdata[frame][line][2])
	    outstring=outstring+newline
	outstring=outstring+"END\n"
    f.write(outstring)
    f.close()

    return 0

def write_psf(outfilename,box,newdata):
    f = open(outfilename, 'w')
    outstring="PSF\n"
    frame0 = newdata[0]
    newline="%8d !NATOM\n" % (len(frame0))
    outstring=outstring+newline
    for line in range(len(frame0)):
	newline="%8d N%03d %4d PSC  N%1d   N%1d                   \n" % ( line+1,line%4+1, (4+line-(line%4))/4,line%4+1,line%4+1 )
	outstring=outstring+newline
    newline="%8d !NBOND\n" % (len(frame0)/2)
    outstring=outstring+newline
    newline=""
    for i in range(len(frame0)/2):
	newline+=" %7d %7d" % (2*i+1,2*i+2)
	if (i%4 ==3):
	    outstring=outstring+newline+"\n"
	    newline=""
    outstring=outstring+"%-64s"%(newline)
    f.write(outstring)
    f.close()

    return 0


def read_input(infilename):
    data=[]
    box=[]
    frame=[]
    inp=open(infilename)
    i=0
    for line in inp:
	linesplit=line.split()
	if (len(linesplit)==1):
	    atomnum = int(linesplit[0])
	    #print atomnum
	    i=0
	else:
	    if (len(linesplit)==6):
		[sweep,num,boxstr,bx,by,bz]=linesplit[:]
		box.append([float(bx),float(by),float(bz)])
	    else:
		[x,y,z,vx,vy,vz,px,py,pz]=linesplit[:]
		frame.append([float(x),float(y),float(z),float(vx),float(vy),float(vz),float(px),float(py),float(pz)])
		i=i+1
	if (i==atomnum):
	    #print frame
	    data.append(frame)
	    frame=[]

    return [box,data]

def vec_normalize(vec):
    # Find the magnitude
    mag2= 0.0
    for val in vec:
        mag2= mag2 + val*val
        mag= math.sqrt(mag2)
    nvec=[]
    # Divide by the magnitude to produce normalized vector
    for val in vec:
        nvec.append(val/mag)
    return nvec

def datatransform(data,leng,patch):
    newdata=[]
    for frame in data:
	newframe=[]
	for line in frame:
	    [x,y,z,vx,vy,vz,px,py,pz]=line[:]
	    vec=vec_normalize([vx,vy,vz])
	    newframe.append([x+leng/2*vec[0],y+leng/2*vec[1],z+leng/2*vec[2]])
	    newframe.append([x-leng/2*vec[0],y-leng/2*vec[1],z-leng/2*vec[2]])
	    patchmove=10*0.5*math.cos((patch+10)/2/180*math.pi)
	    newx=x+leng/2*vec[0]+patchmove*px
	    newy=y+leng/2*vec[1]+patchmove*py
	    newz=z+leng/2*vec[2]+patchmove*pz
	    newframe.append([newx,newy,newz])
	    newx=x-leng/2*vec[0]+patchmove*px
	    newy=y-leng/2*vec[1]+patchmove*py
	    newz=z-leng/2*vec[2]+patchmove*pz
	    newframe.append([newx,newy,newz])
	newdata.append(newframe)

    return newdata

def write_vmd(outfilename,outfilename2,patch):
    f = open("vmd.script", 'w')
    outstring=""
    outstring+="proc setlook {} {\n"
    outstring+="rotate stop\n"
    outstring+="color Display Background white\n"
    outstring+="display projection orthographic\n"
    outstring+="mol delrep 0 0\n"
    outstring+="mol selection \"name N1 N2\"\n"
    outstring+="mol addrep 0\n"
    outstring+="mol selection \"name N3 N4\"\n"
    outstring+="mol addrep 0\n"
    outstring+="mol modstyle 0 0 CPK 10.0 14 20 20\n"
    outstring+="mol modcolor 0 0 ColorID 0\n"
    outstring+="mol modmaterial 0 0 Edgy\n"
    rpatch=10.0*math.sin((patch+10)/2/180*math.pi)
    cylpatch=rpatch*1.4
    outstring+="mol modstyle 1 0 CPK %f %f 20 20\n" %(rpatch,cylpatch)
    outstring+="mol modcolor 1 0 ColorID 1\n"
    outstring+="mol modmaterial 0 0 Edgy\n"
    outstring+="axes location off\n"
    outstring+="}\n"
    outstring+="mol load psf %s \n" %(outfilename2)
    outstring+="mol addfile %s 0\n"%(outfilename)
    outstring+="setlook\n"

    f.write(outstring)
    f.close()

    return 0

def make(infilename,outfilename,outfilename2,leng,patch):

    [box,data]=read_input(infilename)
    newdata=datatransform(data,leng,patch)
    write_pdb(outfilename,box,newdata)
    write_psf(outfilename2,box,newdata)
    write_vmd(outfilename,outfilename2,patch)
    return 0

parser=optparse.OptionParser()
help="""Usage:
%prog [options] 
"""
parser.set_usage(help)
parser.add_option(
    "-i",
    "--input",
    help="Set file from which you want to load data",
    dest="infilename",
    default="movie"
    )
parser.add_option(
    "-o",
    "--output",
    help="Set to which file you want to save data",
    dest="outfilename",
    default="movie.pdb"
    )
parser.add_option(
    "--psf",
    help="Set to which file you want to save connectivity - psf",
    dest="outfilename2",
    default="movie.psf"
    )
parser.add_option(
    "-l",
    "--length",
    help="Set length of spherocylinder",
    dest="leng",
    default="30"
    )
parser.add_option(
    "-p",
    "--patch",
    help="Set size of patch in degrees",
    dest="patch",
    default="90"
    )

(options,arguments)=parser.parse_args()
make(options.infilename,options.outfilename,options.outfilename2,float(options.leng),float(options.patch))
