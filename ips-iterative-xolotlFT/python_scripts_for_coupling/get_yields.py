#!/usr/bin/env python
#=======================================================================================
# get_yields.py
# Calculate sputtering and reflection yield(s) obtained with TRIDYN
# for a single or multiple angle(s)
#=======================================================================================
 
import numpy as np
import math
import numpy.polynomial.polynomial as poly
import sys
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
                      fNImpacts=1.0e5,
                      angle=[0.0],
                      weightAngle=[1.0],
                      logFile=None
                  ):

    if logFile  is not None:
        print(('redirect getYields output of to:', logFile))
        outF = open(logFile, "a")
        sys.stdout = outF

    else:
        print ('no log file defined in getYields')
        print ('print output to default sys.stdout')
        
    print(' ')
    print('getYields:')

    totalSpYield=0.0;
    totalRYield=0.0
    yields=[]


    if len(angle)>1:
        print('\t reading the impact energy distribution for ', (len(angle)), ' angles') 
    else:
        print('\t single, fixed angle used')

    totalWeight=np.sum(weightAngle)
    print('\t the sum of weights is ', totalWeight)
    for a in np.arange(0,len(angle),1):

        if weightAngle[a]==0.0:
            print('\t for angle ', angle[a], 'with weight = ', weightAngle[a], 'skipping all yields analysis')

        elif weightAngle[a] >0.0:
            angleFolder=ftridynFolder+str(angle[a])
            
            #calculate the sputtering yield for each run and take the average
            #if this does not work, use method of He_WSPYIELD.OUT (in xolotl's component)
            ftridynCurrentOutOutput=angleFolder+"/"+ftridynOneOutOutput
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

            print('\t for angle ',angle[a],', weight = ',weightAngle[a],': spY = ',spYield,' ; weighted spY = ',weightedSpYield,' ; rY = ',reflectYield, ' ; weighted rY = ', weightedRYield)
            sys.stdout.flush()

    print(' ')
    print("\t average sputtering yield is ", totalSpYield)
    yields.append(totalSpYield)
    print("\t average reflection yield is ", totalRYield)
    yields.append(totalRYield)

    sys.stdout.flush()

    return yields

################# END OF NEW PYTHON SCRIPT  ####################

if __name__ == '__main__':

   import shutil

   sputtering_and_reflection()
