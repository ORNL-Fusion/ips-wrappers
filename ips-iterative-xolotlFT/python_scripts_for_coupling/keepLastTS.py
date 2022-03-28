#!/usr/bin/env python
#keepLastTS version with function def to call from IPS-FTX

import numpy as np
import math
from   pylab import *
import h5py
import sys

def keepLastTS(inFile='xolotlStop.h5', outFile='netFile', print_test=False):

    print('\t called keepLastTS')
    #if TEST
    if print_test:
        print('\t \t TEST: with input:')
        print('\t \t', inFile)
    sys.stdout.flush()
    
    ## Open the file we want to copy from
    f = h5py.File(inFile, 'r')
    #if TEST
    if print_test:
        print('\t \t succesfully opened ')
        print('\t \t', inFile)
        print(' ')
        print('\t \t TEST: read concentrations of last TS')
    
    ## Get the last time step saved in the file
    concGroup = f['concentrationsGroup']
    timestep = concGroup.attrs['lastTimeStep']

    #if TEST
    if print_test:
        print('\t \t ... read concs succesfully')
        print(' ')
    
    #if TEST
    if print_test:
        print('\t \t TEST: write into outFile: ')
        print('\t \t', outFile)
        print(' ')
        sys.stdout.flush()
    
    ## Create the file to copy to
    fNew = h5py.File(outFile, 'a')

    ## Create the concentration group
    concGroupNew = fNew.create_group('concentrationsGroup')

    ## Set the last time step
    concGroupNew.attrs['lastTimeStep'] = timestep

    ## Copy the last timestep group
    groupName ='concentration_' + str(timestep)
    concGroup.copy(groupName, concGroupNew)

    ## Copy the other groups
    f.copy('headerGroup', fNew)
    f.copy('networkGroup', fNew)
    #if TEST
    if print_test:
        print('\t \t TEST: all information succesfully written into ')
        print('\t \t', outFile)

    print('\t ...keepLastTS done!')
    sys.stdout.flush()

    return
