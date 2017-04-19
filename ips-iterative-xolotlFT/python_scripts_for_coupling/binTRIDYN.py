#!/usr/bin/env python

## Bin the output from Xolotl

import numpy as np
import math
from   pylab import *

def binTridyn(inFile='last_TRIDYN_toBin.dat', outFile='last_TRIDYN.dat'):

## Open the files
    depth, He, V, I = loadtxt(inFile, usecols = (0,1,2,3) , unpack=True)

## Look for indiceCut
    indiceCut = -1
    for k in range(1, len(depth)):
        step = depth[k] - depth[k-1]
        if (step > 0.95):
            indiceCut = int(depth[k-1] * 2) + 1
            break

## Bin the concentrations
    depthBin = []
    heBin = []
    vBin = []
    iBin = []
    for k in range (0, 200):
        depthBin.append(k/2.0)
        heBin.append(0.0)
        vBin.append(0.0)
        iBin.append(0.0)

    oldIndice = -10
    i = 0
    for k in range(0, len(depth)):
        indice = int(floor(depth[k] * 2.0))
        if (indice != oldIndice):
            if (oldIndice >= 0):
                heBin[oldIndice] = heBin[oldIndice] / double(i)
                vBin[oldIndice] = vBin[oldIndice] / double(i)
                iBin[oldIndice] = iBin[oldIndice] / double(i)
            i = 0
            oldIndice = indice
    
        i = i + 1
        if (indice < indiceCut):
            heBin[indice] = heBin[indice] + He[k]
            vBin[indice] = vBin[indice] + V[k]
            iBin[indice] = iBin[indice] + I[k]
        else:
            heBin[indice] = heBin[indice] + He[k]
            heBin[indice-1] = heBin[indice-1] + He[k]
            vBin[indice] = vBin[indice] + V[k]
            vBin[indice-1] = vBin[indice-1] + V[k]
            iBin[indice] = iBin[indice] + I[k]
            iBin[indice-1] = iBin[indice-1] + I[k]
            
## Open 'outputFile.dat' where results will be printed
    outputFile = open(outFile, 'w')

## Loop on all the elements
    for i in range(0, len(depthBin)):
    
    ## Write in the output file
        if (heBin[i] > 0.0):
            outputFile.write("%s %s %s %s\n" %(depthBin[i], heBin[i], vBin[i], iBin[i]))

## Close the output file
    outputFile.close()

    return
