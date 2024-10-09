#!/usr/bin/env python
#=======================================================================================
# get_yields_call.py
# Calculate sputtering and reflection yield(s) obtained with TRIDYN
# for a single or multiple angle(s)
#=======================================================================================
 
import numpy as np
import math
import numpy.polynomial.polynomial as poly
import sys
import pickle
import os

#import matplotlib.pyplot as plt
#from   pylab import spicy # *
#from scipy.optimize import curve_fit
#from   scipy import stats
#from scipy import interpolate

#only implemented for a range of angles as the script is expected to be used for:
#      a) one energy, one angle --> one run (no loop)
#      b) a distribution of angles; a distribution of energies for each angle
#            --> run Ftridyn for each angle, with the given Ein

def sputtering_and_reflection(ftridynOneOutOutput='He_WOUT.dat',
                              ftridynFolder="angle",
                              angle=[0.0],
                              weightAngle=[1.0],
                              fNImpacts=1.0e5,
                              logFile=None,
                              print_test=False
                              ):


    cwd = os.getcwd()

    if logFile  is not None:
        print('\t redirect getYields (call) output of to:')
        print('\t ' , logFile)
        outF = open(logFile, "a")
        sys.stdout = outF

    else:
        print ('\t no log file defined in getYields (call)')
        print ('\t print output to default sys.stdout')

    print(' ')
    sys.stdout.flush()

    print(' ')
    print('\t getYields:')

    if print_test:
        print('\t run with inputs:')
        print('\t \t ftridynOneOutOutput = ', ftridynOneOutOutput)
        print('\t \t ftridynFolder = ', ftridynFolder)
        print('\t \t angle = ', angle)
        print('\t \t weightAngle = ', weightAngle)
        print('\t \t fNImpacts = ', fNImpacts)
        print('\t \t logFile = ', logFile)
        print(' ')
        
    totalSpYield=0.0;
    totalRYield=0.0
    yieldOutput=[]


    if len(angle)>1:
        print('\t \t reading the impact energy distribution for ', (len(angle)), ' angles') 
    else:
        print('\t \t single, fixed angle used')

    totalWeight=np.sum(weightAngle)
    print('\t \t the sum of weights is ', totalWeight)
    for a in np.arange(0,len(angle),1):

        if weightAngle[a]==0.0:
            print('\t \t for angle ', angle[a], 'with weight = ', weightAngle[a], 'skipping all yields analysis')

        elif weightAngle[a] >0.0:
            angleFolder=ftridynFolder+str(angle[a])
            
            #calculate the sputtering yield for each run and take the average
            #if this does not work, use method of He_WSPYIELD.OUT (in xolotl's component)
            ftridynCurrentOutOutput=angleFolder+"/"+ftridynOneOutOutput
            if not os.path.isfile(ftridynCurrentOutOutput):
                continue
            ftridynCurrentOutFile=open(ftridynCurrentOutOutput,"r")
            ftridynCurrentOutData=ftridynCurrentOutFile.read().split('\n')
            searchStringSputter='SPUTTERED PARTICLES(2)'
            for line in ftridynCurrentOutData:
                if searchStringSputter in line:
                    break
            sputterStringWithEmptyFields=line.strip().split(" ")
            sputteringNparticlesString=[x for x in sputterStringWithEmptyFields if x]
            sputteringNparticles=sputteringNparticlesString[2]
            spYield=float(sputteringNparticles)/float(fNImpacts)
            weightedSpYield=spYield*weightAngle[a]/totalWeight
            totalSpYield += weightedSpYield

            #idem for reflection:
            searchStringReflect='BACKSCATTERED PROJECTILES(1)'
            for line in ftridynCurrentOutData:
                if searchStringReflect in line:
                    break
            reflectStringWithEmptyFields=line.strip().split(" ")
            reflectNparticlesString=[x for x in reflectStringWithEmptyFields if x]
            reflectNparticles=reflectNparticlesString[2]
            reflectYield=float(reflectNparticles)/float(fNImpacts)
            weightedRYield=reflectYield*weightAngle[a]/totalWeight
            totalRYield += weightedRYield

            print('\t \t for angle ',angle[a],', weight = ',weightAngle[a],': spY = ',spYield,' ; weighted spY = ',weightedSpYield,' ; rY = ',reflectYield, ' ; weighted rY = ', weightedRYield)
            sys.stdout.flush()

    print(' ')
    print("\t \t average sputtering yield is ", totalSpYield)
    yieldOutput.append(totalSpYield)
    print("\t \t average reflection yield is ", totalRYield)
    yieldOutput.append(totalRYield)

    if print_test:
        print('\t \t yieldOutput = ', yieldOutput)
    
    sys.stdout.flush()

    return yieldOutput
    

################# END OF NEW PYTHON SCRIPT  ####################

if __name__ == '__main__':

   import shutil

   sputtering_and_reflection()
