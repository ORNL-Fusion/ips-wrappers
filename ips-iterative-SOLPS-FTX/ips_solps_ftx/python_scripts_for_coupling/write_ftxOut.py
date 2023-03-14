#!/usr/bin/env python
import os
from ips_xolotlFT.python_scripts_for_coupling import param_handler
import sys
import pickle


def write_ftxOut(grid=20,
                 last_tridyn='last_TRIDYN.dat',
                 log_ftx='log.ftx',
                 ftxIn='ftxInput.txt',
                 retentionFile='retentionOut.txt',
                 print_test=False,
                 logFile=None):


    cwd = os.getcwd()
    print(' ')
    print('\t Called write_ftxOut')
    print('\t from directory: ')
    print('\t', cwd)
    sys.stdout.flush()

    #if pikle file exists, read from pkl file:
    pkl_file=cwd+'/write_ftxOut.pkl'
    if os.path.exists(pkl_file):
        dic = pickle.load( open( pkl_file, "rb" ) )
        #first check the log file, to print everything there
        if 'logFile' in dic:
            logFile=dic['logFile']
        grid=dic['grid']
        last_tridyn=dic['last_tridyn']
        log_ftx=dic['log_ftx']
        ftxIn=dic['ftxIn']
        retentionFile=dic['retentionFile']
        print_test=dic['print_test']
    else:
        print('no pkl file found, continue with function-call-inputs or defaults')

    if logFile  is not None:
        print('\t redirect write_ftxOut output of to:')
        print('\t' , logFile)
        outF = open(logFile, "a")
        sys.stdout = outF

    sys.stdout.flush()

    if print_test:
        print('\t launched script with inputs:')
        print('\t grid = ', grid)
        print('\t last_tridyn = ', last_tridyn)
        print('\t log_ftx = ', log_ftx)
        print('\t ftxIn = ', ftxIn)
        print('\t retentionFile = ', retentionFile)
        print('\t print_test = ', print_test)
        print('\t logFile =', logFile)

    print(' ')
    print('RUN write_ftxOut:')
    print(' ')


    if (print_test):
        print('\t this is just a test')
        print('\t this script is emptry for now')
        print('\t will create a ftxOut.txt file with dummy values')

    
    grid=20    
    Twall=400
    RFT=0.9
    RXol=0.9
    Rtot=RFT+(1-RFT)*RXol
    
    outputFTFile=open('ftxOut.txt', "w")
    outputFTFile.write("#grid  T_wall  R_FT R_Xol  R_tot \n")
    outputFTFile.write(str(grid)+'\t'+str(Twall)+'\t'+str(RFT)+'\t'+str(RXol)+'\t'+str(Rtot)+'\n' )
    outputFTFile.close()

    return

################# END OF NEW PYTHON SCRIPT  ####################

if __name__ == '__main__':

   import shutil

   write_ftxOut()
