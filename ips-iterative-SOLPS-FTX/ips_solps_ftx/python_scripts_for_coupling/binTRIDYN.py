#!/usr/bin/env python

## binTridyn version for tempGrip executable of Xolotl
## Bin the output from Xolotl

import numpy as np
import math
import h5py
import sys
#from   pylab import *

def binTridyn(inFile='TRIDYN_0.h5', outFile='last_TRIDYN.dat'):

    print(' ')
    print('TEST binTRIDYN tempGrid TEST binTRIDYN tempGrid TEST binTRIDYN tempGrid TEST')
    print('from binTRIDYN tempGrid, called with input')
    print(inFile)
    sys.stdout.flush()
    
    ## Open the files
    f = h5py.File(inFile, 'r')
    concDset = f['concs']
    print('from binTRIDYN, opened inFile to read')
    sys.stdout.flush()
    
    ## Bin the concentrations
    print('read concDset')
    depthBin = []
    heBin = []
    vBin = []
    iBin = []
    dBin = []
    tBin = []
    TBin = []
    
    nBins=int(math.floor(concDset[len(concDset)-1][0] * 2.0)+1)

    print('from binTRIDYN, initialized bins')
    print('calculated bins')
    
    for k in range (0, nBins):#200):
        depthBin.append(k/2.0)
        heBin.append(0.0)
        vBin.append(0.0)
        iBin.append(0.0)
        TBin.append(0.0)
        dBin.append(0.0)
        tBin.append(0.0)

    print('from binTRIDYN, zeros in bins')
        
    oldIndice = -10
    i = 0

    for k in range(0, len(concDset)):
        indice = int(math.floor(concDset[k][0] * 2.0))
        if (indice != oldIndice):
            if (oldIndice >= 0):
                heBin[oldIndice] = heBin[oldIndice] / float(i)
                vBin[oldIndice] = vBin[oldIndice] / float(i)
                iBin[oldIndice] = iBin[oldIndice] / float(i)
                TBin[oldIndice] = TBin[oldIndice] / float(i)
                dBin[oldIndice] = dBin[oldIndice] / float(i)
                tBin[oldIndice] = tBin[oldIndice] / float(i)
            i = 0
            oldIndice = indice
    
        i = i + 1
        heBin[indice] = heBin[indice] + concDset[k][1]
        vBin[indice] = vBin[indice] + concDset[k][2]
        iBin[indice] = iBin[indice] + concDset[k][3]
        TBin[indice] = TBin[indice] + concDset[k][4]
        dBin[indice] = dBin[indice] + concDset[k][5]
        tBin[indice] = tBin[indice] + concDset[k][6]


    print('from binTRIDYN, rebinned values')
    print('end of loop')
    sys.stdout.flush()
    
    ## Open 'outputFile.dat' where results will be printed
    outputFile = open(outFile, 'w')

    ## Loop on all the elements
    for i in range(0, len(depthBin)):
    
        ## Write in the output file
        if (TBin[i] > 0.0):
            outputFile.write("%s %s %s %s %s %s %s\n" %(depthBin[i], heBin[i], vBin[i], iBin[i], TBin[i], dBin[i], tBin[i]))

    ## Close the output file
    outputFile.close()

    print(' ... ')
    print('successfully produced output file')
    print(outFile)
    sys.stdout.flush()
    return
