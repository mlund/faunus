#
# Simple script to convert molecular structures into Faunus AAM
# format. Uses OpenBabel's python interface to read a multitude
# of different input formats. Openbabel can be installed in
# anacona using:
# conda install --channel https://conda.anaconda.org/openbabel openbabel

# python 2/3 compatibility
from __future__ import print_function, division

import openbabel as ob
from math import pi
import sys, os, datetime


# see http://openbabel.org/docs/2.3.0/UseTheLibrary/PythonDoc.html

def MolecularWeight(residue):
    Mw = 0
    for atom in ob.OBResidueAtomIter(residue):
        Mw += atom.GetAtomicMass()
    return Mw


def Radius(residue):
    rho = 1.0
    Mw = MolecularWeight(residue)
    return (3. / (4 * pi) * Mw / rho) ** (1 / 3.)


def MassCenter(residue):
    wsum = 0
    v = [0, 0, 0]
    for atom in ob.OBResidueAtomIter(residue):
        w = atom.GetAtomicMass()
        wsum += w
        v[0] += w * atom.x()
        v[1] += w * atom.y()
        v[2] += w * atom.z()
    v[0] /= wsum
    v[1] /= wsum
    v[2] /= wsum
    return v


if len(sys.argv) == 1:
    print("First argument must be a structure file. Supported formats:")
    for s in ob.OBConversion().GetSupportedInputFormat():
        print(s)
    sys.exit()

mol = ob.OBMol()
obconv = ob.OBConversion()
infile = sys.argv[1]
informat = obconv.FormatFromExt(infile)
obconv.SetInFormat(informat)
obconv.ReadFile(mol, infile)

assert mol.NumResidues() > 0, infile + " not found or it is empty."

print("# Infile:", infile, "on", datetime.datetime.now(), os.uname()[1])
print(mol.NumResidues())
for res in ob.OBResidueIter(mol):
    cm = MassCenter(res)
    resname = res.GetName()
    resnum = res.GetNum()
    charge = 0
    radius = Radius(res)
    weight = MolecularWeight(res)
    print('{0:4} {1:5} {2:8.3f} {3:8.3f} {4:8.3f} {5:6.3f} {6:6.2f} {7:6.2f}'.format(
        resname, resnum, cm[0], cm[1], cm[2], charge, weight, radius))
