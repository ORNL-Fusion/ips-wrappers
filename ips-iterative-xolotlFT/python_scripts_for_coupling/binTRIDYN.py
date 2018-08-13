#!/usr/bin/env python

## Bin the output from Xolotl

import numpy as np
import math
import h5py
#from   pylab import *

def binTridyn(inFile='last_TRIDYN_toBin.h5', outFile='last_TRIDYN.dat'):

## Open the files
    f = h5py.File(inFile, 'r')    
    concDset = f['concs']

## Look for indiceCut
    indiceCut = -1
    for k in range(1, len(concDset)):
        step = concDset[k][0] - concDset[k-1][0]
        if (step > 0.95):
            indiceCut = int(concDset[k-1][0] * 2) + 1
            break

## Bin the concentrations
    depthBin = []
    heBin = []
    dBin = []
    tBin = []
    vBin = []
    iBin = []

    nBins=int(math.floor(concDset[len(concDset)-1][0] * 2.0)+1)

    for k in range (0, nBins):#200):
        depthBin.append(k/2.0)
        heBin.append(0.0)
        dBin.append(0.0)
        tBin.append(0.0)
        vBin.append(0.0)
        iBin.append(0.0)

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
            i = 0
            oldIndice = indice
    
        i = i + 1
        if (indice < indiceCut):
            heBin[indice] = heBin[indice] + concDset[k][1]
            dBin[indice] = dBin[indice] + concDset[k][2]
            tBin[indice] = tBin[indice] + concDset[k][3]
            vBin[indice] = vBin[indice] + concDset[k][4]
            iBin[indice] = iBin[indice] + concDset[k][5]
        else:
            heBin[indice] = heBin[indice] + concDset[k][1]
            heBin[indice-1] = heBin[indice-1] + concDset[k][1]
            dBin[indice] = dBin[indice] + concDset[k][2]
            dBin[indice-1] = dBin[indice-1] + concDset[k][2]
            tBin[indice] = tBin[indice] + concDset[k][3]
            tBin[indice-1] = tBin[indice-1] + concDset[k][3]
            vBin[indice] = vBin[indice] + concDset[k][4]
            vBin[indice-1] = vBin[indice-1] + concDset[k][4]
            iBin[indice] = iBin[indice] + concDset[k][5]
            iBin[indice-1] = iBin[indice-1] + concDset[k][5]
            
## Open 'outputFile.dat' where results will be printed
    outputFile = open(outFile, 'w')

## Loop on all the elements
    for i in range(0, len(depthBin)):
    
    ## Write in the output file
        if (heBin[i] > 0.0 or dBin[i] > 0.0 or tBin[i] > 0.0):
            outputFile.write("%s %s %s %s %s %s\n" %(depthBin[i], heBin[i], dBin[i], tBin[i], vBin[i], iBin[i]))

## Close the output file
    outputFile.close()

    return
