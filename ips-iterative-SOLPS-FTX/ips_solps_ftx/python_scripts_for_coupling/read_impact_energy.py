#!/usr/bin/env python                                                                                                                                                                                                 
import os
import sys

#def scalar(time=0.0, inputEnergyFile='hpicEnergies.txt', outputEnergyFile='energyForFTX.txt', print_test=False, logFile=None): #with output file
def scalar(time=0.0, inputEnergyFile='hpicEnergy.txt', print_test=False, logFile=None): #without output file

    if logFile  is not None:
        print('\t redirect read_impact_energy output to:')
        print('\t' , logFile)
        orig_stdout = sys.stdout
        logF = open(logFile, "a")
        sys.stdout = logF
    else: 
        print('\t no logFile defined for read_impact_energy')

    print('From read_impact_energy.scalar')
    if print_test:
        print('\t launched script with inputs:')
        print('\t time = ', time)
        print('\t inputEnergyFile = ', inputEnergyFile)
        #print('\t outputEnergyFile = ', outputEnergyFile)
        print('\t print_test = ', print_test)
        print('\t logFile =', logFile)

    print(' ')
    print('RUN read_impact_energy.scalar:')
    print(' ')
    sys.stdout.flush()

    #search for line with: init_t <= time <= end_t
    inF=open(inputEnergyFile, "r")
    inLines=inF.readlines()
    for x in inLines:
        l=x.strip()
        init_t=float(l.split('\t')[0])
        end_t=float(l.split('\t')[1])
        if (time>=init_t) and (time<=end_t):
            energy_D=float(l.split('\t')[2])
            energy_C=float(l.split('\t')[3])
            print('\t for t=', time, ' init_t = ', init_t, ' end_t=', end_t)
            print('\t found energy_D = ', energy_D, ' and energy_C = ', energy_C, ' eV')
            break;
    inF.close()

    #outF=open(outputEnergyFile, "w") #overwrite previous file
    #outF.write(str(energy_C, energy_D))
    #outF.close()
    #sys.stdout.flush()

    sys.stdout.flush()
    if logFile  is not None:
        if print_test:
            print('\t close log file ', logF.name)
            print('\t set sys.stdout back to: ', orig_stdout)
        sys.stdout.flush()
        sys.stdout = orig_stdout
        logF.close()
        
    return [energy_D,energy_C]


def distrib(time=0.0, fileName=None, print_test=False, logFile=None):
    #NOT IMPLEMENTED YET
    print('this function is defined, but not implemeted yet')
    print('returns nothing!')
    sys.stdout.flush()
    return 



### END OF PYTHON SCRIPT
if __name__ == '__main__':

    import shutil

    scalar()
