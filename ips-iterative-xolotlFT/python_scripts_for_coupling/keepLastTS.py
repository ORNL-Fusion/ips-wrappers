#!/usr/bin/env python
#keepLastTS version with function def to call from IPS-FTX

import numpy as np
import math
from   pylab import *
import h5py
import sys

def keepLastTS(inFile='xolotlStop.h5', outFile='netFile'):

    print(' ')
    print('keepLastTS, called with input:')
    print(inFile)
    sys.stdout.flush()
    
    ## Open the file we want to copy from
    f = h5py.File(inFile, 'r')
    print('\t succesfully opened ', inFile)
    print('\t read concentrations of last TS')
    
    ## Get the last time step saved in the file
    concGroup = f['concentrationsGroup']
    timestep = concGroup.attrs['lastTimeStep']

    print('\t read concs succesfully')
    print('\t write them into outFile: ', outFile)
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
    print('\t all information succesfully written into ', outFile)
    print('keepLastTS done!')
    sys.stdout.flush()

    return
