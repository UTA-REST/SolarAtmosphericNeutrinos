#!/bin/env python

import argparse
import numpy as np
import os 

BaseDir=os.environ['SOLARATMHOME']
MCEqPath=BaseDir+"/MCEq"

parser = argparse.ArgumentParser(description='Argument parsing for python')
parser.add_argument('--bvalues',dest='bvalues', type=float, default=[0.0],nargs="+",
                   help='what b values to run')
parser.add_argument('--dens_divs',dest='dens_divs', type=int,default=200,
                   help='how many slices of density to use')
parser.add_argument('--dens_index',dest='dens_index', type=int,default=2,
                   help='bin density power index')
parser.add_argument('--atmosphere',dest='atmosphere',default=BaseDir+"/AtmProfiles/SolarAtmosphereCut.csv",
                    help="atmosphere to use")
parser.add_argument('--hadronicmodel',dest='hadronicmodel',default='SIBYLL2.3_rc1',
                   help='hadronic model to use')
parser.add_argument('--primarymodel',dest='primarymodel',default='HillasGaisser_H3a',
                   help='primary model to use')
parser.add_argument('--outputpath',dest='outputpath',default=BaseDir+"/Outputs",
                    help='path of output directory')
parser.add_argument('--runname',dest='runname',default='testrun',
                    help='name of this run')

args=parser.parse_args()

print "Configuration of this run:"
print args

AtmospherePath=args.atmosphere
bvalues=       args.bvalues
runname=       args.runname
outputpath=    args.outputpath
DensDivs=      args.dens_divs
DensIndex=     args.dens_index

Settings={}
Settings['MakeFiles']=True
Settings['HadronicModel']=args.hadronicmodel
Settings['PrimaryModel']=args.primarymodel


import os
import sys
from os.path import join
os.chdir(MCEqPath)
sys.path.append(MCEqPath+"/MCEqRunner")

from scipy import interpolate as interp

from MCEq.core import MCEqRun
from mceq_config import config, mceq_config_without
import MCEq.geometry as geom
import MCEqRunF


#Load atmospheric model and set constants
kmtocm=1000*100
SolarRadius=696*pow(10,3)*kmtocm
Atmosph=np.loadtxt(AtmospherePath,usecols=[1,2,3],delimiter=',')
SolarHeight=(Atmosph[:,0]-Atmosph[-1,0])*kmtocm
Densities=Atmosph[:,2]
ChromosHeight=SolarHeight[0]
MaxB=(ChromosHeight+SolarRadius)/SolarRadius
LogDensFunction=interp.interp1d(SolarHeight, np.log10(Densities),bounds_error=False)
DensFunction=lambda x: pow(10.,LogDensFunction(x))

TestPoints=MCEqRunF.SetupTestPoints(bvalues,DensFunction,SolarRadius,ChromosHeight,DensDivs,BinSchemeIndex=DensIndex)



mceqs=[]
for i in range(0,len(TestPoints)):
    Settings['FileName']=outputpath+"/"+runname+"_b{:1.4f}".format(TestPoints[i][0])
#    print "Running test point " +str(i) + " output to " + Name 
    MCEqRunF.RunMCEq(TestPoints[i][1],TestPoints[i][2],Settings)

