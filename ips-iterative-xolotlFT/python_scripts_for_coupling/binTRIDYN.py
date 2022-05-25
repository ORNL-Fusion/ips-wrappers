#!/usr/bin/env python

import numpy as np
import math
import h5py
import sys
import os.path
from os import path

def v1(inFile='last_TRIDYN_toBin.h5', outFile='last_TRIDYN.dat', print_test=False):

## binTridyn version for v1 (e.g., master) executable of Xolotl
## Bin the output from Xolotl v1
    
    #print(' ')
    print('\t called binTRIDYN_master')

    #if TEST
    if print_test:
        print('\t \t TEST: with input')
        print('\t \t \t', inFile)
        print ("\t \t File exists:", str(path.exists(inFile)), '\n')
    sys.stdout.flush()
    
## Open the files
    f = h5py.File(inFile, 'r')    
    #if TEST
    if print_test:
        print('\t \t TEST: from binTRIDYN, opened inFile to read')
    concDset = f['concs']
    #if TEST
    if print_test:
        print('\t \t read concDset \n')
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

    #if TEST
    if print_test:
        print('\t \t TEST: from binTRIDYN, initialized bins')
        sys.stdout.flush()
    
    for k in range (0, nBins):#200):
        depthBin.append(k/2.0)
        heBin.append(0.0)
        dBin.append(0.0)
        tBin.append(0.0)
        vBin.append(0.0)
        iBin.append(0.0)
        TBin.append(0.0)

    #if TEST
    if print_test:
        print('\t \t from binTRIDYN, zeros in bins')
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

    #if TEST
    if print_test:
        print('\t \t from binTRIDYN, rebinned values')
        print('\t \t end of loop \n')
        sys.stdout.flush()
        
## Open 'outputFile.dat' where results will be printed
    outputFile = open(outFile, 'w')
    
    #if TEST
    if print_test:
        print('\t \t TEST: from binTRIDYN, write output to output file \n')
        sys.stdout.flush()
    
## Loop on all the elements
    for i in range(0, len(depthBin)):
    
    ## Write in the output file
        #if (TBin[i] > 0.0):
        if (heBin[i] + vBin[i] + iBin[i] + dBin[i] + tBin[i] + TBin[i] > 0.0):
            outputFile.write("%s %s %s %s %s %s %s\n" %(depthBin[i], heBin[i], dBin[i], tBin[i], vBin[i], iBin[i], TBin[i]))

## Close the output file
    outputFile.close()

    print('\t ... binTRIDYN_master successfully produced output file')
    print('\t \t ', outFile)
    #print(' ')
    sys.stdout.flush()
    
    return


def v2(inFile='TRIDYN_0.h5', outFile='last_TRIDYN.dat', print_test=False):

## binTridyn version for v2 (e.g. tempGrip) executable of Xolotl
## Bin the output from Xolotl v2
    
    #print(' ')
    print('\t called binTRIDYN_tempGrid')
    
    if print_test:
        print('\t \t TEST: with input')
        print('\t \t ', inFile)
        sys.stdout.flush()
    
    ## Open the files
    f = h5py.File(inFile, 'r')
    concDset = f['concs']
    if print_test:
        print('\t \t from binTRIDYN, opened inFile to read \n')
        sys.stdout.flush()
    
    ## Bin the concentrations
    if print_test:
        print('\t \t TEST: read concDset \n')
    depthBin = []
    heBin = []
    vBin = []
    iBin = []
    dBin = []
    tBin = []
    TBin = []
    
    nBins=int(math.floor(concDset[len(concDset)-1][0] * 2.0)+1)

    if print_test:
        print('\t \t TEST: from binTRIDYN, initialized bins')
        print('\t \t calculated bins')
    
    for k in range (0, nBins):#200):
        depthBin.append(k/2.0)
        heBin.append(0.0)
        vBin.append(0.0)
        iBin.append(0.0)
        TBin.append(0.0)
        dBin.append(0.0)
        tBin.append(0.0)

    if print_test:
        print('\t \t from binTRIDYN, zeros in bins')
        
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

    if print_test:
        print('\t \t from binTRIDYN, rebinned values')
        print('\t \t end of loop \n')
    sys.stdout.flush()
    
    ## Open 'outputFile.dat' where results will be printed
    outputFile = open(outFile, 'w')

    if print_test:
        print('\t \t TEST: write output to output file \n')
    
    ## Loop on all the elements
    for i in range(0, len(depthBin)):
    
        ## Write in the output file
        #if (TBin[i] > 0.0):
        if (heBin[i] + vBin[i] + iBin[i] + dBin[i] + tBin[i] + TBin[i] > 0.0):
            outputFile.write("%s %s %s %s %s %s %s\n" %(depthBin[i], heBin[i], vBin[i], iBin[i], TBin[i], dBin[i], tBin[i]))

    ## Close the output file
    outputFile.close()

    print(' \t ... binTRIDYN_tempGrid successfully produced output file')
    print('\t \t', outFile)
    #print(' ')
    sys.stdout.flush()
    return
