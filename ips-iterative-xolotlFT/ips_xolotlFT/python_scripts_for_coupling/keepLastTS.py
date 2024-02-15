#!/usr/bin/env python
#keepLastTS version with function def to call from IPS-FTX

import numpy as np
import math
from   pylab import *
import h5py
import sys

def keepLastTS(inFile='xolotlStop.h5', outFile='netFile', print_test=False):

    #if TEST
    if print_test:
        print('\t called keepLastTS')
        print('\t \t with input:')
        print('\t \t', inFile)
    sys.stdout.flush()
    
    ## Open the file we want to copy from
    f = h5py.File(inFile, 'r')
    #if TEST
    if print_test:
        print('\t \t succesfully opened ')
        print('\t \t', inFile)
        print(' ')
        print('\t \t read concentrations of last TS')
    
    ## Get the last time step saved in the file
    concGroup = f['concentrationsGroup']
    timestep = concGroup.attrs['lastTimeStep']
    if print_test:
        print(f"\t \t => Foud attribute 'lastTimeStep' with value {timestep}")

    ## Get the last loop saved in the file, if any
    if "lastLoop" in concGroup.attrs.keys():
        loop = concGroup.attrs['lastLoop']
        if print_test:
            print(f"\t \t => Foud attribute 'lastLoop' with value {loop}")

    #if TEST
    if print_test:
        print('\t \t ... read concs succesfully')
        print(' ')
    
    #if TEST
    if print_test:
        print('\t \t write into outFile: ')
        print('\t \t', outFile)
        print(' ')
        sys.stdout.flush()
    
    ## Create the file to copy to
    fNew = h5py.File(outFile, 'a')

    ## Create the concentration group
    concGroupNew = fNew.create_group('concentrationsGroup')

    ## Set the last time step
    concGroupNew.attrs['lastTimeStep'] = timestep
    if print_test:
        print(f"\t \t => Set attribute 'lastTimeStep' to {timestep}")

    ## Set the last loop if needed
    if "lastLoop" in concGroup.attrs.keys():
        concGroupNew.attrs['lastLoop'] = loop
        if print_test:
            print(f"\t \t => Set attribute 'lastLoop' to {loop}")

    ## Copy the last timestep group
    groupName = 'concentration_' + str(timestep)
    if not groupName in concGroup.keys():
        groupName = 'concentration_' + str(loop) + '_' + str(timestep)
    if not groupName in concGroup.keys():
        print(f"ERROR: group with name {groupName} not found in concentration group!")
    if print_test:
        print(f"\t \t => Group name of last time step is {groupName}")
    concGroup.copy(groupName, concGroupNew)

    ## Copy the other groups
    if "headerGroup" in f.keys():
        f.copy('headerGroup', fNew)
    f.copy('networkGroup', fNew)
    #if TEST
    if print_test:
        print('\t \t all information succesfully written into ')
        print('\t \t', outFile)

        print('\t ...keepLastTS done!')
        sys.stdout.flush()

    return
