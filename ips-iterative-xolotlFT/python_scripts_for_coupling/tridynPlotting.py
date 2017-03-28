#!/usr/bin/env python
#=======================================================================================
# tridynPlotting.py
# Plots and fit the profile obtained with TRIDYN
#=======================================================================================
 
import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *
from scipy.optimize import curve_fit
from   scipy import stats
from scipy import interpolate
import numpy.polynomial.polynomial as poly
 
## Create plots
#fig1 = plt.figure()
#plot1 = plt.subplot(111)
 
## Open files ("bla1" is not used but I could not figure out how to easily open a file without using two columns)
depth1, bla1 = loadtxt('He_WDUMPPRJ.dat', usecols = (2,3) , unpack=True)
 
## Put first data into the plot
#n, bins, patches = plt.hist(depth1/10.0, 200, normed=True, facecolor=plt.cm.jet(200), alpha=0.4)
n, bins = np.histogram(depth1/10.0, 200,normed=True)
## "bins" is actually the edges on the histogram bins so it needs to be formated to give the center of each bin in "b"
b = np.delete(bins, len(bins) - 1)
offset = (b[1] - b[0]) / 2.0
for i in range(len(b)):
    b[i] = b[i] + offset
 
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
