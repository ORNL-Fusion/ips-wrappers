#!/usr/bin/env python
#=======================================================================================
# tridynPlotting.py
# Plots and fit the profile obtained with TRIDYN
#=======================================================================================
 
import numpy as np
import math
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
#from   pylab import spicy # *
#from scipy.optimize import curve_fit
#from   scipy import stats
#from scipy import interpolate

#only implemented for a range of angles as the script is expected to be used for:
#      a) one energy, one angle --> one run (no loop)
#      b) a distribution of angles; a distribution of energies for each angle
#            --> run Ftridyn for each angle, with the given Ein

def ftridyn_to_xolotl(ftridynOneOutput='He_WDUMPPRJ.dat',
                      ftridynFolder="angle",                      
                      GITRoutput='impactDistribution.dat'
                  ):
    
    
    #read GITR output to get weight -- MODIFY WHEN I'VE GOT AN EXAMPLE
    angleIndex, weightAngle = np.loadtxt(GITRoutput, usecols = (0,1) , unpack=True)
    print "read GITRs output, %s" %(GITRoutput)
    print "number of angles is ", (len(angleIndex))
    
    for a in np.arange(1,len(angleIndex)+1,1):
        
        numLines== sum(1 for line in open('ftridynOneOutput'))
        num_lines[]
        num_lines.append = sum(1 for line in open('myfile.txt'))

        print "---------------------"
        print "folder number=", a
        print "value of the angle= ", angleIndex[a-1]
        print "weight of angle ", a , " is ", weightAngle[a-1]
        print "and the file is of length",
    
        ## Open files ("bla1" is not used but I could not figure out how to easily open a file without using two columns)
        angleFolder=ftridynFolder+"_%d" %(a)
        ftridynCurrentOutput=angleFolder+"/"+ftridynOneOutput

        nLines = sum(1 for line in open('ftridynCurrentOutput'))
        numLines[]
        numLines.append = nLines
        print "FTridyn's file has " , nLines, " lines"
        if nLines==0:
            continue
            
        print "reading file: %s" %(ftridynCurrentOutput)
        depth1, bla1 = np.loadtxt(ftridynCurrentOutput, usecols = (2,3) , unpack=True)
        
        ## Put first data into the plot
        m, bins = np.histogram(depth1/10.0, 200,normed=True)

        ## "bins" is actually the edges on the histogram bins so it needs to be formated to give the center of each bin in "b"
        b = np.delete(bins, len(bins) - 1)
        offset = (b[1] - b[0]) / 2.0
        for i in range(len(b)):
            b[i] = b[i] + offset
 
        #weight each binned distribution and add contribution from each angle / energy
        #m is the values of each run; n is the weighted sum of m's
        #initialize array of total distribution, n
        if (a==1):
            n=[]
            for i in range(len(m)):
                n.append(0)

        for i in range(len(m)):
            n[i]+=m[i]*weightAngle[a-1]/numLines[a-1]
            
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


    return


################# END OF NEW PYTHON SCRIPT  ####################

if __name__ == '__main__':

   import shutil

   ftridyn_to_xolotl()
