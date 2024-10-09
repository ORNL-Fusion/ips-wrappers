#!/usr/bin/env python
#=======================================================================================
# translate_ftridyn_to_xolotl_call.py
# fits the profile obtained with TRIDYN
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

def ftridyn_to_xolotl(ftridynOnePrjOutput='He_WDUMPPRJ.dat',                  
                      ftridynFolder="angle",       
                      angle=[0.0],
                      weightAngle=[1.0],
                      nBins=200,
                      fitOrder=15,
                      prjRange=50.0, #in [A]
                      logFile=None,
                      print_test=False
                      ):

    cwd = os.getcwd()

    if logFile  is not None:
        outF = open(logFile, "a")
        sys.stdout = outF
        print('\t redirect tridynPlotting output of to:')
        print('\t ', logFile)
        
    else:
        print ('\t no log file defined in tridynPlotting')
        print ('\t print output to default sys.stdout')
    print(' ')
    sys.stdout.flush()

    print('\t tridynPlotting:')

    if print_test:
        print('\t run with inputs:')
        print('\t \t ftridynOnePrjOutput = ', ftridynOnePrjOutput)
        print('\t \t ftridynFolder = ', ftridynFolder)
        print('\t \t angle = ', angle)
        print('\t \t weightAngle = ', weightAngle)
        print('\t \t nBins = ', nBins)
        print('\t \t prjRange = ', prjRange)
        print('\t \t logFile = ', logFile)
        print(' ')
        
    totalSpYield=0.0;
    totalRYield=0.0
    yields=[]

    if len(angle)>1:
        print('\t \t reading the impact energy distribution for ', (len(angle)), ' angles') 
    else:
        print('\t \t single, fixed angle used')

    totalWeight=np.sum(weightAngle)
    print('\t \t the sum of weights is ', totalWeight, ' and projectile range', prjRange , ' [A]')

    nonZeroAngle=0
    for a in np.arange(0,len(angle),1):
        if weightAngle[a]==0.0:
            print('\t \t \t angle',a,' of weight ', weightAngle[a], ': skip all analysis with no contribution to implantation profile')

        elif weightAngle[a] >0.0:
            angleFolder=ftridynFolder+str(angle[a])
    
            ## Open files ("bla1" is not used but I could not figure out how to easily open a file without using two columns)
            ftridynCurrentPrjOutput=angleFolder+"/"+ftridynOnePrjOutput
            #print "reading file: %s" %(ftridynCurrentPrjOutput)
            num_lines_prj = sum(1 for line in open(ftridynCurrentPrjOutput))

            if num_lines_prj==0:
                print("\t \t \t WARNING: prj file is empty for angle ", angle[a])
            else:
                if print_test:
                    print("\t \t calculate sputtering yield for angle ", angle[a])
                nonZeroAngle+=1 #identify the first time that an angle contributes, to initialize n[]
                depth1, bla1 = np.loadtxt(ftridynCurrentPrjOutput, usecols = (2,3) , unpack=True)

                ## Put first data into the plot; '/10.0' is to convert A -> nm                       
                #print 'in translate script: using range', prjRange , ' [A]'
                m, bins = np.histogram(depth1/10.0, bins=nBins, range=(0.0,prjRange/10.0), density=True)

                ## "bins" is actually the edges on the histogram bins so it needs to be formated to give the center of each bin in "b"
                b = np.delete(bins, len(bins) - 1)
                offset = (b[1] - b[0]) / 2.0
                for i in range(len(b)):
                    b[i] = b[i] + offset

                #weight each binned distribution and add contribution from each angle / energy
                #m is the values of each run; n is the weighted sum of m's
                #initialize n[] the first time that an angle contributes (when nonZeroAngle=1)
                if nonZeroAngle==1:
                    n=[]
                    if print_test:
                        print("\t \t initialized n[] for angle ", angle[a])
                    for i in range(len(m)):
                        n.append(0)

                for i in range(len(m)):
                    n[i]+=m[i]*weightAngle[a] #/numLines[a] ; no need to normalize by number of lines, as that's caused by reflection

                sys.stdout.flush()

    ## Fit with polynomials
    fit = poly.polyfit(b, n, fitOrder)

    ## Get the fit function
    fitFunc = poly.Polynomial(fit)

    ## Open 'tridyn.dat' where results will be printed
    outputFile = open('tridyn.dat', 'w')
    
    ## Write in the output file
    for i in range(fitOrder+1):
        #outputFile.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \n" %(fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6], fit[7], fit[8], fit[9], fit[10], fit[11], fit[12], fit[13], fit[14], fit[15]))
        outputFile.write("%s "%(fit[i]))
    for i in range(15-fitOrder):
        outputFile.write("0.0 ")
    outputFile.write("\n")
        
    ## Close the output file
    outputFile.close()
    
    sys.stdout.flush()
    return

################# END OF NEW PYTHON SCRIPT  ####################

if __name__ == '__main__':

   import shutil

   ftridyn_to_xolotl()
