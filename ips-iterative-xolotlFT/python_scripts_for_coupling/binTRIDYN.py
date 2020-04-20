#!/usr/bin/env python

## Bin the output from Xolotl

import numpy as np
import math
import h5py
import sys
import os.path
from os import path
#from   pylab import *

def binTridyn(inFile='last_TRIDYN_toBin.h5', outFile='last_TRIDYN.dat'):

    print(' ')
    print('from binTRIDYN, called with input')
    print(inFile)
    sys.stdout.flush()

    print ("File exists:"+str(path.exists(inFile)))
    
## Open the files
    f = h5py.File(inFile, 'r')    
    print('from binTRIDYN, opened inFile to read')
    concDset = f['concs']
    print('read concDset')
    sys.stdout.flush()

## Bin the concentrations
    depthBin = []
    heBin = []
    dBin = []
    tBin = []
    vBin = []
    iBin = []
    TBin = []

    nBins=int(math.floor(concDset[len(concDset)-1][0] * 2.0)+1)

    print('from binTRIDYN, initialized bins')
    print('calculated bins')
    sys.stdout.flush()
    
    for k in range (0, nBins):#200):
        depthBin.append(k/2.0)
        heBin.append(0.0)
        dBin.append(0.0)
        tBin.append(0.0)
        vBin.append(0.0)
        iBin.append(0.0)
        TBin.append(0.0)

    print('from binTRIDYN, zeros in bins')
    sys.stdout.flush()
        
    oldIndice = -10
    i = 0

    for k in range(0, len(concDset)):
        indice = int(math.floor(concDset[k][0] * 2.0))
        if (indice != oldIndice):
            if (oldIndice >= 0):
                heBin[oldIndice] = heBin[oldIndice] / float(i)
                dBin[oldIndice] = dBin[oldIndice] / float(i)
                tBin[oldIndice] = tBin[oldIndice] / float(i)
                vBin[oldIndice] = vBin[oldIndice] / float(i)
                iBin[oldIndice] = iBin[oldIndice] / float(i)
                TBin[oldIndice] = TBin[oldIndice] / float(i)
            i = 0
            oldIndice = indice
    
        i = i + 1

        heBin[indice] = heBin[indice] + concDset[k][1]
        dBin[indice] = dBin[indice] + concDset[k][2]
        tBin[indice] = tBin[indice] + concDset[k][3]
        vBin[indice] = vBin[indice] + concDset[k][4]
        iBin[indice] = iBin[indice] + concDset[k][5]
        TBin[indice] = TBin[indice] + concDset[k][6]

    print('from binTRIDYN, rebinned values')
    print('end of loop')
    sys.stdout.flush()
        
## Open 'outputFile.dat' where results will be printed
    outputFile = open(outFile, 'w')

    print('from binTRIDYN, write output to output file')
    sys.stdout.flush()
    
## Loop on all the elements
    for i in range(0, len(depthBin)):
    
    ## Write in the output file
        if (TBin[i] > 0.0):
            outputFile.write("%s %s %s %s %s %s %s\n" %(depthBin[i], heBin[i], dBin[i], tBin[i], vBin[i], iBin[i], TBin[i]))

## Close the output file
    outputFile.close()

    print(' ... ')
    print('successfully produced output file')
    print(outFile)
    print(' ')
    sys.stdout.flush()
    
    return
