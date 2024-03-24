#!/usr/bin/env python
import os
from ips_xolotlFT.python_scripts_for_coupling import param_handler
import sys
import pickle


def write_ftxOut(grid=20,
                 last_tridyn='last_TRIDYN.dat',
                 log_ftx='log.ftx',
                 tridyn='tridyn.dat',
                 retentionFile='retentionOut.txt',
                 H_plasma='no',
                 outFile='ftxOut.txt',
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
        tridyn=dic['tridyn']
        retentionFile=dic['retentionFile']
        H_plasma=dic['H_plasma']
        outFile=dic['outFile']
        print_test=dic['print_test']
    else:
        print('no pkl file found, continue with function-call-inputs or defaults')

    if logFile  is not None:
        print('\t redirect write_ftxOut output of to:')
        print('\t' , logFile)
        outF = open(logFile, "a")
        sys.stdout = outF

    sys.stdout.flush()

    print('From write_ftxOut')
    if print_test:
        print('\t launched script with inputs:')
        print('\t grid = ', grid)
        print('\t last_tridyn = ', last_tridyn)
        print('\t log_ftx = ', log_ftx)
        print('\t tridyn = ', tridyn)
        print('\t retentionFile = ', retentionFile)
        print('\t H_plasma = ', H_plasma)
        print('\t print_test = ', print_test)
        print('\t logFile =', logFile)

    print(' ')
    print('RUN write_ftxOut:')
    print(' ')


    #if (print_test):
    #    print('\t this is just a test')
    #    print('\t this script is emptry for now')
    #    print('\t will create a ftxOut.txt file with dummy values')

    
    #grid=20    
    #Twall=400
    #RFT=0.9
    #RXol=0.9
    #Rtot=RFT+(1-RFT)*RXol
    
    #T wall : 1st line, column 6 of last_TRIDYN.dat:
    Twall_line=open(last_tridyn).readline().rstrip()
    if (print_test):
        print('\t Twall_line = ', Twall_line)
    Twall=float(Twall_line.split(' ')[5])
    print('\t Twall = ', Twall)
    print(' ')

    
    #R_FT : "D :  spY = spY and rY = "  $R_FT
    #we only want the last occurrence of the string (latest value of RFT)
    #--> open file in reverse & the 1st time we find the string
    R_FT_lines=reversed(open(log_ftx).readlines())
    if H_plasma=='yes':
        R_FT_string="H :  spY"
        if (print_test):
            print('\t looking for R_FT_string for H : ', R_FT_string)
    else:
        R_FT_string="D :  spY"        
        if (print_test):
            print('\t looking for R_FT_string for D : ', R_FT_string)    
    for row in R_FT_lines:
        #if R_FT_string in R_FT_lines:
        if row.find(R_FT_string) != -1:
            if (print_test):
                print('\t found R_FT_string in row ', row)
                print("\t row.split(' ') = ", row.split(' '))
            RFT = float(row.split(' ')[9])
            break;
    try:
        print('\t R_FT = ', RFT)
    except:
        print('\t WARNING! did not find string, assume RFT = 0.5')
        RFT=0.5
    print(' ')

    #D fraction
    #needed to calculate Rxol
    #in tridyn.dat, line with "D" column 3 
    D_frac_lines=open(tridyn).readlines()
    if (print_test):
        print('\t looking for string (D) in ', tridyn)
    for row in D_frac_lines:
        if row.find("D ") != -1:
            if (print_test):
                print('\t found D line', row)
                print('\t row.split(' ') =', row.split(' '))
            D_frac=float(row.split(' ')[2])
            break;
    try:
        print('\t D_frac = ', D_frac)
    except:
        print('\t WARNING! did not find string, assume D_frac = 1')
        D_frac=1
    print(' ')
    
    #R_Xol:
    #old method: 1 - (D_content / fluence*D_fraction)
    #            fluence, column 2: 1st line  - last line
    #            D_content, column 4: 1st line  - last line
    #new method: 1 - (D_content / D_fluence)
    #            D fluence, column 3: 1st line  - last line
    #            D_content, column 5: 1st line  - last line
    
    retLines=open(retentionFile).readlines()
    R_Xol_firstLine=retLines[1]
    if (print_test):
        print('\t first line of retention file is : ', R_Xol_firstLine)
        print('\t R_Xol_firstLine.split(' ') = ', R_Xol_firstLine.split(' '))    
    #init_fluence=float(R_Xol_firstLine.split(' ')[1]) #this is total fluence ; old method
    init_fluence=float(R_Xol_firstLine.split(' ')[2]) #new method
    init_Dconc=float(R_Xol_firstLine.split(' ')[4])
    print('\t initial D fluence = ', init_fluence, ' and D_conc = ', init_Dconc)
    print(' ')
    R_Xol_lastLine=retLines[-1]
    if (print_test):
        print('\t last line of retention file is : ', R_Xol_lastLine)
        print('\t R_Xol_lastLine.split(' ') = ', R_Xol_lastLine.split(' '))
    #last_fluence=float(R_Xol_lastLine.split(' ')[1]) #this is total fluence ; old method 
    last_fluence=float(R_Xol_lastLine.split(' ')[2]) #new method
    last_Dconc=float(R_Xol_lastLine.split(' ')[4])
    print('\t last D fluence = ', last_fluence, ' and D_conc = ', last_Dconc)
    print(' ')
    fluence=last_fluence-init_fluence
    D_conc=last_Dconc-init_Dconc
    print('\t total fluence =', fluence, ' and D_conc = ', D_conc)
    #RXol=1-(D_conc/(fluence*D_frac)) # old method
    if fluence>0:
        RXol=1-(D_conc/fluence)
        print('\t RXol = ', RXol)
    else:
        print('\t WARNING: fluence <= 0 ; RXol = 1')
        RXol=1
    print(' ') 
    
    print('\t calculate Rtot=RFT+(1-RFT)*RXol, using RFT = ', RFT, ', RXol = ', RXol)
    Rtot=RFT+(1-RFT)*RXol
    print('\t Rtot = ', Rtot)
    print(' ')
    
    print('\t and write to ', outFile)
    outputFTFile=open(outFile, "a+")
    #outputFTFile.write("#grid \t  T_wall \t R_FT \t R_Xol \t R_tot \n") skip the header for now
    outputFTFile.write(str(grid)+'\t'+str(Twall)+'\t'+str(RFT)+'\t'+str(RXol)+'\t'+str(Rtot)+'\n' )
    outputFTFile.close()
    print('Done with write_ftxOut')
    print(' ')

    return

################# END OF NEW PYTHON SCRIPT  ####################

if __name__ == '__main__':

   import shutil

   write_ftxOut()
