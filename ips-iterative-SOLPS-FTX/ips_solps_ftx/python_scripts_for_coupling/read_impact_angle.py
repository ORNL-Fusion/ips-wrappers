#!/usr/bin/env python                                                                                                                                                                                                 
import os
import sys

#def scalar(time=0.0, inputAngleFile='hpicAngles.txt', outputAngleFile='angleForFTX.txt', print_test=False, logFile=None): #with output file
def scalar(time=0.0, inputAngleFile='hpicAngles.txt', print_test=False, logFile=None): #without output file

    if logFile  is not None:
        print('\t redirect read_impact_angle output to:')
        print('\t' , logFile)
        orig_stdout = sys.stdout
        logF = open(logFile, "a")
        sys.stdout = logF
    else: 
        print('\t no logFile defined for read_impact_angle')

    print('From read_impact_angle.scalar')
    if print_test:
        print('\t launched script with inputs:')
        print('\t time = ', time)
        print('\t inputAngleFile = ', inputAngleFile)
        #print('\t outputAngleFile = ', outputAngleFile)
        print('\t print_test = ', print_test)
        print('\t logFile =', logFile)

    print(' ')
    print('RUN read_impact_angle.scalar:')
    print(' ')
    sys.stdout.flush()

    #search for line with: init_t <= time <= end_t
    inF=open(inputAngleFile, "r")
    inLines=inF.readlines()
    for x in inLines:
        l=x.strip()
        init_t=float(l.split('\t')[0])
        end_t=float(l.split('\t')[1])
        if (time>=init_t) and (time<=end_t):
            angle_D=float(l.split('\t')[2])
            angle_C=float(l.split('\t')[3])
            print('\t for t=', time, ' init_t = ', init_t, ' end_t=', end_t)
            print('\t found angle_D = ', angle_D , ' and angle_C = ', angle_C, ' deg')
            break;
    inF.close()

    #outF=open(outputAngleFile, "w") #overwrite previous file
    #outF.write(str(angle))
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
        
    return [angle_D,angle_C]


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
