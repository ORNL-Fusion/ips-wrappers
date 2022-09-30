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

def sputtering_and_reflection_launch(ftridynOneOutOutput='He_WOUT.dat',
                                     ftridynFolder="angle",
                                     angle=[0.0],
                                     weightAngle=[1.0],
                                     fNImpacts=1.0e5,
                                     logFile=None,
                                     print_test=False
):

    #If pikle file exists, read from pkl file:
    cwd = os.getcwd()
    pkl_file=cwd+'/ft_getYields.pkl'
    if os.path.exists(pkl_file):
        #dic = pickle.load( open( pkl_file, "rb" ) )
        with open(pkl_file, "rb") as pf:
            dic = pickle.load(pf)
        #first check the log file, to print everything there
        if 'logFile' in dic:
            logFile=dic['logFile']

    if logFile  is not None:
        print('\t redirect getYields output of to:')
        print('\t ' , logFile)
        outF = open(logFile, "a")
        sys.stdout = outF

    else:
        print ('\t no log file defined in getYields')
        print ('\t print output to default sys.stdout')

    print(' ')
    sys.stdout.flush()

    #if pikle file exists, read from pkl file:
    if os.path.exists(pkl_file):
        #dic = pickle.load( open( pkl_file, "rb" ) ) #already above
        print('\t In get_yields,  dictionary defined in pkl file: ')
        print('\t', pkl_file)
        print(' ')
        sys.stdout.flush()
        if 'print_test' in dic:
            print_test=dic['print_test']

        if print_test:
            print('\t TEST: dictionary in pkl file contains:')
            print('\t \t', dic)
        
        #for each of the possible inputs to the script, check if given in pkl file

        #print lines if in test mode:
        if 'ftridynOneOutOutput' in dic:
            ftridynOneOutOutput=dic['ftridynOneOutOutput']
            if print_test:
                print('\t from pkl file, set dict value to ftridynOneOutOutput=',dic['ftridynOneOutOutput'])
        else:
            print('\t no value defined in pkl; use default for ftridynOneOutOutput=', ftridynOneOutOutput)

        if 'ftridynFolder' in dic:
            ftridynFolder=dic['ftridynFolder']
            if print_test:
                print('\t from pkl file, set dict value to ftridynFolder=', dic['ftridynFolder'])
        else:
            print('\t no value defined in pkl; use default for ftridynFolder=', ftridynFolder)

        if 'angle' in dic:
            angle=dic['angle']
            if print_test:
                print('\t from pkl file, set dict value to angle=', dic['angle'])
        else:
            print('\t no value defined in pkl; use default for angle=', angle)

        if 'weightAngle' in dic:
            weightAngle=dic['weightAngle']
            if print_test:
                print('\t from pkl file, set dict value to weightAngle=', dic['weightAngle'])
        else:
            print('\t no value defined in pkl; use default for weightAngle=', weightAngle)

        if 'fNImpacts' in dic:
            fNImpacts=dic['fNImpacts']
            if print_test:
                print('\t from pkl file, set dict value to fNImpacts=', dic['fNImpacts'])
        else:
            print('\t no value defined in pkl; use default for fNImpacts=', fNImpacts)

        if 'logFile' in dic:
            if print_test:
                print('\t from pkl file, set dict value to logFile=', dic['logFile'])
        else:
            print('\t no value defined in pkl; use default for logFile=', logFile)

    else:
        print('get_yields did not find pkl file, run with default values')
    sys.stdout.flush()

    #close file after reading values so that we can use it to print yields
    pf.close()
    
    print(' ')
    print('\t getYields:')

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
        print('\t \t TEST: yieldOutput = ', yieldOutput)
    
    sys.stdout.flush()

    #if pkl file exists, use it to return yields
    if os.path.exists(pkl_file):
        yields_dic={}
        yields_dic['yields']=yieldOutput
        with open(pkl_file, "wb") as pf:
            pickle.dump(yields_dic, pf)
        pf.close()
        return
    else: #if pkl file / dic not present, try returning yields directly
        return yieldOutput
    

################# END OF NEW PYTHON SCRIPT  ####################

if __name__ == '__main__':

   import shutil

   sputtering_and_reflection_launch()
