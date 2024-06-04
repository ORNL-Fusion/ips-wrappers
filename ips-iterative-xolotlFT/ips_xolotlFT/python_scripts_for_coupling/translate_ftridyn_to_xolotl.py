#!/usr/bin/env python
#=======================================================================================
# tridynPlotting.py
# Plots and fit the profile obtained with TRIDYN
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

def ftridyn_to_xolotl_launch(ftridynOnePrjOutput='He_WDUMPPRJ.dat',                  
                             ftridynFolder="angle",       
                             angle=[0.0],
                             weightAngle=[1.0],
                             nBins=200,
                             fitOrder=15,
                             prjRange=50.0, #in [A]
                             logFile=None,
                             print_test=False
                  ):

    #if pikle file exists, read from pkl file: 
    cwd = os.getcwd()
    pkl_file=cwd+'/ft_implProfiles.pkl'
    if os.path.exists(pkl_file):
        dic = pickle.load( open( pkl_file, "rb" ) )
        #first check the log file, to print everything there
        if 'logFile' in dic:
            logFile=dic['logFile']

    if logFile  is not None:
        outF = open(logFile, "a")
        sys.stdout = outF
        print('\t redirect tridynPlotting output of to:')
        print('\t ', logFile)
        
    else:
        print ('\t no log file defined in tridynPlotting')
        print ('\t print output to default sys.stdout')
    print(' ')
    #cwd = os.getcwd() #already above
    sys.stdout.flush()

    if 'print_test' in dic:
        print_test=dic['print_test']

    if print_test:
        print('\t dictionary in pkl file contains:')
        print('\t \t', dic)

    #if pikle file exists, read from pkl file:   
    #pkl_file=cwd+'/ft_implProfiles.pkl' #already above 
    if os.path.exists(pkl_file):
        #dic = pickle.load( open( pkl_file, "rb" ) ) #already above
        print('\t In translate_ft_to_xol, dictionary defined in pkl file: ')
        print('\t', pkl_file)
        sys.stdout.flush()

        print(' ')
        #print lines if in test mode:
        #for each of the possible inputs to the script, check if given in pkl file
        if 'ftridynOnePrjOutput' in dic:
            ftridynOnePrjOutput=dic['ftridynOnePrjOutput']
            if print_test:
                print('\t from pkl file, set dict value to ftridynOnePrjOutput=',dic['ftridynOnePrjOutput'])
        else:
            print('\t no value defined in pkl; use default for ftridynOnePrjOutput=', ftridynOnePrjOutput)

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
            
        if 'nBins' in dic:
            nBins=dic['nBins']
            if print_test:
                print('\t from pkl file, set dict value to nBins=', dic['nBins'])
        else:
            print('\t no value defined in pkl; use default for nBins=', nBins)

        if 'fitOrder' in dic:
            fitOrder=dic['fitOrder']
            if print_test:
                print('\t from pkl file, set dict value to fitOrder=', dic['fitOrder'])
        else:
            print('\t no value defined in pkl; use default for fitOrder=', fitOrder)
            
        if 'prjRange' in dic:
            prjRange=dic['prjRange']
            if print_test:
                print('\t from pkl file, set dict value to prjRange=', dic['prjRange'])
        else:
            print('\t no value defined in pkl; use default for prjRange=', prjRange)
            
        if 'logFile' in dic:
            if print_test:
                print('\t from pkl file, set dict value to logFile=', dic['logFile'])
        else:
            print('\t no value defined in pkl; use default for logFile=', logFile)

    else:
        print('In translate_ft_to_xol, did not find pkl file, run with default values')
    sys.stdout.flush()

    print(' ')
    print('\t tridynPlotting:')

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
    
    ### Plot the fit
    #y = []
    #for i in range(len(b)):
    #    y.append(fitFunc(b[i]))
    #plot1.plot(b, y, lw=5, color=plt.cm.jet(200), label='pure W')
    # 
    ## Plot the legend
    #l = plot1.legend(loc='best')
    #setp(l.get_texts(), fontsize=40)
    # 
    ### Some shaping
    #plot1.set_xlabel("depth (nm)",fontsize=35)
    ##plot1.set_ylabel("A.U.",fontsize=25)
    #plot1.set_xlim([0.0, 12.0])
    #plot1.set_ylim([0.0, 0.25])
    #plt.xticks(np.arange(0, 12.0, 5.0))
    ## plot1.set_yscale('log')
    ##plot1.grid()
    #plot1.tick_params(axis='both', which='major', labelsize=30)
    #plot1.tick_params(axis='both', which='minor', labelsize=30)
    # 
    ### Show the plots
    #plt.show()

    sys.stdout.flush()
    return

################# END OF NEW PYTHON SCRIPT  ####################

if __name__ == '__main__':

   import shutil

   ftridyn_to_xolotl_launch()
