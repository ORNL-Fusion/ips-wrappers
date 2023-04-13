#!/usr/bin/env python                                                                                                                                                                                                 
import os
import sys

def average_SOLPS_input(ftxOut_file='ftxOut.txt', print_test=False, logFile=None):

    # SOLPS can only take average the FTX recycling factors
    # move to script when looks like it works                                                                                                                                                                 

    if logFile  is not None:
        print('\t redirect write_ftxOut output of to:')
        print('\t' , logFile)
        outF = open(logFile, "a")
        sys.stdout = outF


    print('From average_SOLPS_input')
    if print_test:
        print('\t launched script with inputs:')
        print('\t ftxOut_file = ', ftxOut_file)
        print('\t print_test = ', print_test)
        print('\t logFile =', logFile)

    print(' ')
    print('RUN average_SOLPS_input:')
    print(' ')
    sys.stdout.flush()
    
    t_wall=[]
    recycling_FT=[]
    recycling_xol=[]
    recycling_total=[]
    f=open(ftxOut_file, "r")
    lines=f.readlines()

    
    for x in lines:
        l=x.strip()
        t_wall.append(float(l.split('\t')[1]))
        recycling_FT.append(float(l.split('\t')[2]))
        recycling_xol.append(float(l.split('\t')[3]))
        recycling_total.append(float(l.split('\t')[4]))

    f.close()

    try:
        ave_twall=sum(t_wall)/len(t_wall)
        ave_Rft=sum(recycling_FT)/len(recycling_FT)
        ave_Rxol=sum(recycling_xol)/len(recycling_xol)
        ave_Rtot=sum(recycling_total)/len(recycling_total)
        if print_test:
            print('succesfully calculated average FTX outputs are:')
            print('\t average T wall = ', ave_twall)
            print('\t average R_FT = ', ave_Rft)
            print('\t average R_Xol = ', ave_Rxol)
            print('\t average R_tot = ', ave_Rtot)
            print('Done with average_SOLPS_input')
            sys.stdout.flush()

            return  ave_twall, ave_Rft, ave_Rxol, ave_Rtot
    except:
        print("some average calculation went wrong:")
        if print_test:
            print("\t Twall sum = ", sum(t_wall), " and len = ", len(t_wall))
            print("\t Twall sum = ", sum(recycling_FT), " and len = ", len(recycling_FT))
            print("\t Twall sum = ", sum(recycling_xol), " and len = ", len(recycling_xol))
            print("\t Twall sum = ", sum(recycling_total), " and len = ", len(recycling_total))
        print("do not return any values")
        print('Done with average_SOLPS_input')
        sys.stdout.flush()
        return

### END OF PYTHON SCRIPT
if __name__ == '__main__':

    import shutil

    average_SOLPS_input()
