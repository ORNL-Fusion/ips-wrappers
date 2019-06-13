#!/usr/bin/env python
#=======================================================================================
# tridynPlotting.py
# Plots and fit the profile obtained with TRIDYN
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

def ftridyn_to_xolotl(ftridynOnePrjOutput='He_WDUMPPRJ.dat',                  
                      ftridynFolder="angle",       
                      angle=[0.0],
                      weightAngle=[1.0],
                      nBins=200,
                      prjRange=50.0, #in [A]
                      logFile=None
                  ):

    if logFile  is not None:
        print(('redirect tridynPlotting output of to:', logFile))
        outF = open(logFile, "a")
        sys.stdout = outF

    else:
        print ('no log file defined in tridynPlotting')
        print ('print output to default sys.stdout')
        
    print(' ')
    print('tridynPlotting:')

    totalSpYield=0.0;
    totalRYield=0.0
    yields=[]


    if len(angle)>1:
        print('\t reading the impact energy distribution for ', (len(angle)), ' angles') 
    else:
        print('\t single, fixed angle used')

    totalWeight=np.sum(weightAngle)
    print('\t the sum of weights is ', totalWeight, ' and projectile range', prjRange , ' [A]')
    
    nonZeroAngle=0
    for a in np.arange(0,len(angle),1):
        if weightAngle[a]==0.0:
            print('\t \t angle',a,' of weight ', weightAngle[a], ': skip all analysis with no contribution to implantation profile')

        elif weightAngle[a] >0.0:
            angleFolder=ftridynFolder+str(angle[a])
            
            #calculate the sputtering yield for each run and take the average
            ## Open files ("bla1" is not used but I could not figure out how to easily open a file without using two columns)
            ftridynCurrentPrjOutput=angleFolder+"/"+ftridynOnePrjOutput

            #        nLines = sum(1 for line in open(ftridynCurrentPrjOutput,"r"))
            #        if (a==1):
            #            numLines=[]
            
            #        numLines.append = nLines
            #        print "FTridyn's file has " , nLines, " lines"
            #        if nLines==0:
            #            continue
            
            #print "reading file: %s" %(ftridynCurrentPrjOutput)
            num_lines_prj = sum(1 for line in open(ftridynCurrentPrjOutput))
            if num_lines_prj==0:
                print("\t \t WARNING: prj file is empty for angle ", angle[a])
            else:
                #print "TEST: calculate sputtering yield for angle ", angle[a]
                nonZeroAngle+=1 #identify the first time that an angle contributes, to initialize n[]
                depth1, bla1 = np.loadtxt(ftridynCurrentPrjOutput, usecols = (2,3) , unpack=True)
                                
                ## Put first data into the plot; '/10.0' is to convert A -> nm
                #print 'in translate script: using range', prjRange , ' [A]'
                m, bins = np.histogram(depth1/10.0, nBins,(0.0,prjRange/10.0),normed=True)
            
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
                    #print "TEST: initialized n[] for angle ", angle[a]
                    for i in range(len(m)):
                        n.append(0)

                for i in range(len(m)):
                    n[i]+=m[i]*weightAngle[a] #/numLines[a] ; no need to normalize by number of lines, as that's caused by reflection
            
                sys.stdout.flush()

    ## Fit with polynomials
    fit = poly.polyfit(b, n, 15)
    ## Get the fit function
    fitFunc = poly.Polynomial(fit)
    
    ## Open 'tridyn.dat' where results will be printed
    outputFile = open('tridyn.dat', 'w')
    
    ## Write in the output file
    outputFile.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \n" %(fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6], fit[7], fit[8], fit[9], fit[10], fit[11], fit[12], fit[13], fit[14], fit[15]))
    
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

   ftridyn_to_xolotl()
