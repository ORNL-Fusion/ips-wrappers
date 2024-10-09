#!/usr/bin/env python

import os
from ips_xolotlFT.python_scripts_for_coupling import param_handler
import sys
import pickle
import numpy as np

def plasmaOut2ftxIn(plasmaOutFile='plasmaOut.pkl', ftxInFile='ftxIn.txt',print_test=True,logFile=None):

    cwd = os.getcwd()
    print(' ')
    print('\t Called plasmaOut2ftxIn')
    print('\t from directory: ')
    print('\t', cwd)
    sys.stdout.flush()
    
    if logFile  is not None:
        print('\t redirect plasmaOut2ftxIn output of to:')
        print('\t' , logFile)
        orig_stdout = sys.stdout
        outF = open(logFile, "a")
        sys.stdout = outF
    sys.stdout.flush()
    
    if print_test:
        print('Called script with inputs:')
        print('\t plasmaOutFile = ', plasmaOutFile)
        print('\t ftxInFile = ', ftxInFile)
        print('\t print_test = ', print_test)
        print('\t logFile = ', logFile)

        
    print(' ')
    print('RUN plasmaOut2ftxIn:')
    print(' ')

    solpsParams=pickle.load(open(plasmaOutFile, "rb" ))
    if print_test:
        print('\t Reading output of SOLPS:')
        print('\t', solpsParams)
        print(' ')
    sys.stdout.flush()
    
    #solpsParams=param_handler.read(plasmaOutFile) #(INPUT_DIR+'/'+solps_outFile)
    #print('\t reading output of SOLPS:')
    #for k,v, in solpsParams.items():
    #    print(('\t {0} : {1}'.format(k, v)))
    #print(' ')
    #sys.stdout.flush()
    
    #format float to list if needed
    #currently implemented for inputEnergy, Ti, Te, inputAngle, bfieldAngle
    #might need to add Z (all same values or add 1.0?)
    
    #idea is: if solpsParams[k] is a single value (float or int),
    #turn into list with one Ein value per species (all same value)
    #otherwise (if not float/int, no plasmaSpecies given, etc), leave as is
    print('\t Format Ein, Ti, Te, Ain and/or Bin as lists of len(plasmaSpecies) :')
    if ('inputEnergy' in solpsParams):
        Ein=[]
        #TO-DO!! can be re-written as isinstance(param, (int, float)) , together
        if (isinstance(solpsParams['inputEnergy'],float) or isinstance(solpsParams['inputEnergy'],int) or (isinstance(solpsParams['inputEnergy'],np.ndarray) and (solpsParams['inputEnergy'].size < len(solpsParams['plasmaSpecies'])))):
            if ('plasmaSpecies' in solpsParams):
                for n in range(len(solpsParams['plasmaSpecies'])):
                    Ein.append(float(solpsParams['inputEnergy']))
            else: #without plasmaSpecies, just leave as is
                Ein=solpsParams['inputEnergy']
        else: #if it's not float or int, leave as is (whether it's a list or not)
            Ein=solpsParams['inputEnergy']
        print('\t Ein =', Ein)

    if ('Ti' in solpsParams):
        Ti=[]
        if (isinstance(solpsParams['Ti'],float) or isinstance(solpsParams['Ti'],int) or (isinstance(solpsParams['Ti'],np.ndarray) and (solpsParams['Ti'].size < len(solpsParams['plasmaSpecies'])))):
            if ('plasmaSpecies' in solpsParams):
                for n in range(len(solpsParams['plasmaSpecies'])):
                    #add numpy.squeeze to FIX:
                    #Conversion of an array with ndim > 0 to a scalar is deprecated, and will error in future.
                    #Ensure you extract a single element from your array before performing this operation. (Deprecated NumPy 1.25.)                    
                    solpsParams['Ti']=np.squeeze(solpsParams['Ti'])
                    print('\t converted array Ti into scalar ', solpsParams['Ti'])
                    Ti.append(float(solpsParams['Ti']))
            else: #without plasmaSpecies, just leave as is
                Ti=solpsParams['Ti']
        else: #if it's not float or int, leave as is (whether it's a list or not)
            Ti=solpsParams['Ti']
        print('\t Ti =', Ti)

    if ('Te' in solpsParams):
        Te=[]
        if (isinstance(solpsParams['Te'],float) or isinstance(solpsParams['Te'],int) or (isinstance(solpsParams['Te'],np.ndarray) and (solpsParams['Te'].size < len(solpsParams['plasmaSpecies'])))):
            if ('plasmaSpecies' in solpsParams):
                for n in range(len(solpsParams['plasmaSpecies'])):
                    #add numpy.squeeze to FIX:
                    #Conversion of an array with ndim > 0 to a scalar is deprecated, and will error in future.
                    #Ensure you extract a single element from your array before performing this operation. (Deprecated NumPy 1.25.)
                    solpsParams['Te']=np.squeeze(solpsParams['Te'])
                    print('\t converted array Te into scalar ', solpsParams['Te'])
                    Te.append(float(solpsParams['Te']))
            else: #without plasmaSpecies, just leave as is
                Te=solpsParams['Te']
        else: #if it's not float or int, leave as is (whether it's a list or not)
            Te=solpsParams['Te']
        print('\t Te =', Te)

    if ('inputAngle' in solpsParams):
        Ain=[]
        if (isinstance(solpsParams['inputAngle'],float) or isinstance(solpsParams['inputAngle'],int) or (isinstance(solpsParams['inputAngle'],np.ndarray) and (solpsParams['inputAngle'].size < len(solpsParams['plasmaSpecies'])))):
            if ('plasmaSpecies' in solpsParams):
                for n in range(len(solpsParams['plasmaSpecies'])):
                    Ain.append(solpsParams['inputAngle'])
            else: #without plasmaSpecies, just leave as is
                Ain=solpsParams['inputAngle']
        else: #if it's not float or int, leave as is (whether it's a list or not)
            Ain=solpsParams['inputAngle']
        print('\t Ain =', Ain)


    if ('bfieldAngle' in solpsParams):
        Bin=[]
        if (isinstance(solpsParams['bfieldAngle'],float) or isinstance(solpsParams['bfieldAngle'],int) or (isinstance(solpsParams['bfieldAngle'],np.ndarray) and (solpsParams['bfieldAngle'].size < len(solpsParams['plasmaSpecies'])))):
            if ('plasmaSpecies' in solpsParams):
                for n in range(len(solpsParams['plasmaSpecies'])):
                    Bin.append(solpsParams['bfieldAngle'])
            else: #without plasmaSpecies, just leave as is
                Bin=solpsParams['bfieldAngle']
        else: #if it's not float or int, leave as is (whether it's a list or not)
            Bin=solpsParams['bfieldAngle']
            print('\t Bin =', Bin)
        print('\n')


    sys.stdout.flush()
    ## 2 - do ops to calculate/format inputs for FTX:
    print ('\t Estimates for SOLPS output --> FTX input')
    ftxInputs={}

    #if print_test:
    #    print('\t solps output dict = ', solpsParams)
    #sys.stdout.flush()
    
    print('\n')
    ## 2.a : check for parameters expected from SOLPS & 'translate' accordingly
    
    ### 2.a.i: parameters with a single value
    
    ### particle flux: /m2s --> /nm2s done in ftx driver
    if ('flux' in solpsParams):
        ftxInputs['flux'] = solpsParams['flux'] #/1.0e18
        if print_test:
            print('\t ftx[flux] = ', ftxInputs['flux'])
    sys.stdout.flush()

    ### heat flux: W/m2 --> W/nm2
    ### in some (older) version of Xolotl, heat is one line; in others is two lines
    ### for now implement single version suited for xolotl-tempGrid-build (print warning)
    ### might need to make changes to Xolotl param template too (tbd)
    if 'heatFlux' in solpsParams:
        print ('\t \t WARNING: heat flux in Xolotl handled in the newer way, using 2 lines: temperature model & values')
        ftxInputs['tempHandler'] = 'heat'
        #if solps heatflux contains 2 fields, assume it's bulkT given, then assume room T
        if (len(solpsParams['heatFlux']) == 2):
            ftxInputs['tempParam'] = (solpsParams['heatFlux'][0], solpsParams['heatFlux'][1]) #/m2s --> /nm2s done in ftx driver
        #if a single field given, assume bulkT = room T
        elif (len(solpsParams['heatFlux']) == 1):
            ftxInputs['tempParam'] = (solpsParams['heatFlux'], 300.0)  #/m2s --> /nm2s done in ftx driver
            print ('\t \t WARNING: no bulkT given; assume room T = 300 K')
        if print_test:
            print('\t \t heat flux tempHandler =', ftxInputs['tempHandler'], 'and temparam =', ftxInputs['tempParam'])
        print('\n')
    sys.stdout.flush()
            
    ### 2.a.ii: inputs that are given for each species:
    # flux fraction
    # input energy
    # input angle
    
    ### for now, include all
    if ('plasmaSpecies' in  solpsParams):
        #ftxInputs['plasmaSpecies'] = 'He W D T C'
        #C_f : different options for C in F-TRIDYN: _a, _b, _d... ; _f noted as "C for fusion"
        #ftxInputs['plasmaSpecies'] = ['He','W', 'D', 'T', 'C_f']
        std_plasmaSpecies = ['He','W', 'D', 'T'] #output from SOLPS will contain C, not C_f ; handle that in FTRIDYN
        ftxInputs['plasmaSpecies']=std_plasmaSpecies
        print('\t hard coded (std) plasma species in ftx = ', std_plasmaSpecies) #, ftxInputs['plasmaSpecies'] )
        print('\t Append any other that might be in SOLPS to ftxInputs[plasmaSpecies]')
        for n in range(len(solpsParams['plasmaSpecies'])):
            if solpsParams['plasmaSpecies'][n] in ftxInputs['plasmaSpecies']:
                print('\t TEST: ', solpsParams['plasmaSpecies'][n] , ' already exists in ftxInputs[plasmaSpecies]')
                print('\t TEST: skip adding')
            else:
                print('\t TEST: ', solpsParams['plasmaSpecies'][n] , ' does not exists in ftxInputs[plasmaSpecies]')
                print('\t TEST: Add to list so Ein, Ain... are also computed for it')
                ftxInputs['plasmaSpecies'].append(solpsParams['plasmaSpecies'][n])
        print('\t TEST: final ftxInputs[plasmaSpecies] = ', ftxInputs['plasmaSpecies'])
        print('\n')
        
        #create a variable that's 1 if species exists ; 0 if it doesn't:
        #probably won't need it, but leave it commented out for now:
        #for n in range(len(ftxInputs['plasmaSpecies'])):
        #    s = ftxInputs['plasmaSpecies'][n]
        #    if (s in solpsParams['plasmaSpecies']):
        #        speciesExists[n]=1.0
        #    else:
        #        speciesExists[n]=0.0
        
    sys.stdout.flush()
    ### reformat to have a vaolue for all species, even if 0.0
    ### might not be the most efficient method, but it should work for all scenarios
    ### example of what I'm trying to do in loop below
       #if ('He' in solpsParams['plasmaSpecies']):
       #    ftxInputs['fluxFraction'][0] = solpsParams['fluxFraction'][nCount-1]
       #    nCount+=1
       #else:
       #    ftxInputs['fluxFraction'][0] = 0.0
       #
       #if ('W' in solpsParams['plasmaSpecies']):
       #    ftxInputs['fluxFraction'][1] = solpsParams['fluxFraction'][nCount-1]
       #else:
       #    ftxInputs['fluxFraction'][1] = 0.0
       
    print('\t reformat to have Ein, Ain for each species:')
    nCount=0
    for n in range(len(std_plasmaSpecies)): #use std list to avoid C vs C_f conflict #range(len(ftxInputs['plasmaSpecies'])):
        s = std_plasmaSpecies[n] #ftxInputs['plasmaSpecies'][n]
        print('\t for n  = ', n , ' ; s = ', s, ' :')
        
        #initialize dictionary entries:
        #create dictionary entries for the first species
        if (n==0):
            if print_test:
                print('\t \t initialize Ein, Ti, Te, Ain, Bin lists')
            #always write flux fraction, even if zero:
            ftxInputs['fluxFraction']=[] #[0.0] DEL
            if print_test:
                print('\t \t as of loop for ', s, 'ftxInputs[fluxFraction] = ', ftxInputs['fluxFraction'])
            #write inputEnergy if input provided by SOLPS:
            if  (('inputEnergy' in solpsParams) or ('Ti' in solpsParams) or ('Te' in solpsParams)):
                ftxInputs['inputEnergy']=[] #[0.0] DEL
                if print_test:
                    print('\t \t as of loop for ', s, 'ftxInputs[inputEnergy] = ', ftxInputs['inputEnergy'])
            #write inputAngle if input provided by SOLPS:
            if (('inputAngle' in solpsParams) or ('bfieldAngle' in solpsParams)):
                ftxInputs['inputAngle']=[] #[0.0] DEL
                if print_test:
                    print('\t \t as of loop for ', s, 'ftxInputs[inputAngle] = ', ftxInputs['inputAngle'])
            sys.stdout.flush()
            print('\n')

        #append values to dictionary entry in subsequent species
        if (s in solpsParams['plasmaSpecies']):
            nCount+=1
            m=nCount-1
            #flux fraction
            print('\t check for fluxFraction:')
            if ('fluxFraction' in solpsParams):
                ftxInputs['fluxFraction'].append(solpsParams['fluxFraction'][m]) #[n] = solpsParams['fluxFraction'][m] DEL
            else:
                print('\t \t WARNING: no fluxFraction provided by input file. set all to zero')
                ftxInputs['fluxFraction'].append(0.0) #[n] = 0.0 DEL
            print('\t \t ftx[fluxFraction] = ', ftxInputs['fluxFraction'])
            print('\n')
            sys.stdout.flush()

	    #inputEnergy
            #check whether passed as energy, or need to estimate based on Te & Ti
            print('\t check for inputEnergy:')
            if  ('inputEnergy' in solpsParams):
                if print_test:
                    print('\t \t solpsParams[inputEnergy][m] = ', Ein[m]) # solpsParams['inputEnergy'][m]) DEL
                    sys.stdout.flush()
                ftxInputs['inputEnergy'].append(Ein[m]) #[n] = solpsParams['inputEnergy'][m] DEL
                print('\t \t ftx[inputEnergy] = ', ftxInputs['inputEnergy'])
                print('\n')
            elif ('Ti' in solpsParams):
                if ('Te' in solpsParams):
                    if ('Z' in solpsParams):
                        print('\t \t WARNING: no inputEnegy given ; estimate as 2Ti + 3ZTe')
                        if print_test:
                            print('\t \t solpsParams[Ti] = ', Ti[m]) #solpsParams['Ti'][m]) DEL
                            print('\t \t solpsParams[Z] = ', solpsParams['Z'][m]) #may need to change to Z[m]
                            print('\t \t solpsParams[Te] = ', Te[m]) #solpsParams['Te'][m]) DEL
                            sys.stdout.flush() 
                        ftxInputs['inputEnergy'].append(2*Ti[m]+3*solpsParams['Z'][m]*Te[m])
                        #[n] = 2*solpsParams['Ti'][m]+3*solpsParams['Z'][m]*solpsParams['Te'][m]
                    else:
                        print('\t \t WARNING: no inputEnegy given or Z given ; estimate as 2Ti + 3Te (Z=1)')
                        if print_test:
                            print('\t \t solpsParams[Ti] = ', Ti[m]) #solpsParams['Ti'][m]) DEL
                            print('\t \t solpsParams[Te] = ', Te[m]) #solpsParams['Te'][m]) DEL
                            sys.stdout.flush() 
                        ftxInputs['inputEnergy'].append(2*Ti[m]+3*Te[m])
                        #[n] = 2*solpsParams['Ti'][m]+3*solpsParams['Te'][m]
                elif ('Te' not in solpsParams):
                    if ('Z' in solpsParams):
                        print('\t \t WARNING: no inputEnegy or Te given ; estimate as Ti * (2+3Z)')
                        if print_test:
                            print('\t \t solpsParams[Ti] = ', Ti[m]) #solpsParams['Ti'][m]) DEL
                            print('\t \t solpsParams[Z] = ', solpsParams['Z'][m]) #may need to change to Z[m]
                            sys.stdout.flush() 
                        ftxInputs['inputEnergy'].append(Ti*(2+3*solpsParams['Z'][m]))
                        #[n] = solpsParams['Ti'][m]*(2+3*solpsParams['Z'][m])
                    else:
                        print('\t \t WARNING: no inputEnegy, Te or Z given ; estimate as 5*Ti (Z=1)')
                        if print_test:
                            print('\t \t solpsParams[Ti] = ', Ti[m]) #solpsParams['Ti'][m]) DEL
                            sys.stdout.flush() 
                        ftxInputs['inputEnergy'].append(5*Ti[m])
                        #[n] = 5*solpsParams['Ti'][m]
                print('\t \t using Ti-Te, ftx[inputEnergy] = ', ftxInputs['inputEnergy'])
                print('\n')
            elif ('Te' in solpsParams):
                if ('Z' in solpsParams):
                    print('\t \t WARNING: no inputEnegy or Ti given ; estimate as Te * (2+3Z)')
                    if print_test:
                        print('\t \t solpsParams[Z] = ', solpsParams['Z'][m]) #may need to change to Z[m]
                        print('\t \t solpsParams[Te] = ', Te[m]) #solpsParams['Te'][m]) DEL
                        sys.stdout.flush() 
                    ftxInputs['inputEnergy'].append(Te*(2+3*solpsParams['Z'][m]))
                    #[n] = solpsParams['Te'][m]*(2+3*solpsParams['Z'][m])
                else:
                    print('\t \t WARNING: no inputEnegy, Ti or Z given ; estimate as 5*Te (Z=1)')
                    if print_test:
                        print('\t \t solpsParams[Te] = ', solpsParams['Te'][m])
                        sys.stdout.flush() 
                    ftxInputs['inputEnergy'].append(5*Te[m])
                    #[n] = 5*solpsParams['Te'][m]
                print('\t \t using Te, ftx[inputEnergy] = ', ftxInputs['inputEnergy'])
                print('\n')
            else:
                print ('\t \t WARNING: no inputEnergy, Te or Ti given ; continue with default values in config file')
                print('\n')
            sys.stdout.flush()

            #inputAngle
            print('\t check for inputAngle:')
            if ('inputAngle' in solpsParams):
                ftxInputs['inputAngle'].append(Ain[m])
                #[n] = solpsParams['inputAngle'][m]
                print('\t \t ftx[inputAngle] = ', ftxInputs['inputAngle'])
                print('\n')
            elif ('bfieldAngle' in solpsParams):
                print ('\t \t WARNING: no inputAngle given; use B field angle ')
                ftxInputs['inputAngle'].append(Bin[m])
                print('\t \t ftx[inputAngle] = ', ftxInputs['inputAngle'])
                print('\n')
                #[n] = solpsParams['bfieldAngle'][m]
            else:
                print('\t \t WARNING: no inputAngle or bfieldAngle given; continue with default values in config file (likely normal incidence)')
                print('\n')
            sys.stdout.flush()

        else: # if s not in solpsParams['plasmaSpecies']
            print('\t ', s, ' is not in solps[plasmaSpecies] ; set fluxFraction, Ein, Ain to zero')
            #fluxFraction
            if ('fluxFraction' in solpsParams):
                ftxInputs['fluxFraction'].append(0.0) #[n] = 0.0
                print('\t \t fluxFraction = 0.0')
            else:
                print('\t \t WARNING: no fluxFraction provided by input file. set all to zero')
                ftxInputs['fluxFraction'].append(0.0) #[n] = 0.0
            #inputEnergy
            if  (('inputEnergy' in solpsParams) or ('Ti' in solpsParams) or ('Te' in solpsParams)):
                ftxInputs['inputEnergy'].append(0.0) #[n] = 0.0
                print('\t \t inputEnergy = 0.0')
            else:
                print ('\t \t WARNING: no inputEnergy, Te or Ti given ; continue with default values in config file')
            #inputAngle
            if (('inputAngle' in solpsParams) or ('bfieldAngle' in solpsParams)):
                ftxInputs['inputAngle'].append(0.0) #[n] = 0.0
                print('\t \t inputAngle = 0.0')
            else:
                print('\t \t WARNING: no inputAngle or bfieldAngle given; continue with default values in config file (likely normal incidence)')
            print('\n')
            sys.stdout.flush()
        
    ## 2.b : write down all other parameters 'as is':
                        
    print('\t write down all other parameters "as is":')
    for k,v, in solpsParams.items():
        
        if k in ftxInputs:
            if print_test:
                print('\t \t k = ', k , ' ; v = ', v)
                sys.stdout.flush()
            print('\t \t ftx[',k,'] = ', ftxInputs[k])
        else:
            print('\t \t param ', k , 'not in ftx dictionary yet ; add as is')
            if print_test:
                print('\t \t v = ', v)
            ftxInputs[k] = v
    print('\n')
    sys.stdout.flush()

    # 3 - print ftxIn.txt - DONE
    print('\t write FTX input file:')
    print('\t \t to file name: ', ftxInFile)
    ftxIn=open(ftxInFile, "w")
    for k,v, in ftxInputs.items():
        #first check if its numpy array ; if so, turn into list (for Te and Ti, bc sqeeze above)
        if isinstance(v,np.ndarray):
            print('\t \t ', k, ' is a numpy array. Convert to list/ float / int before writing to file')
            v_list=v.tolist()
            v=v_list
        
        if (isinstance(v, int) or isinstance(v, float) or isinstance(v, dict) or isinstance(v, str)):
            print(('\t \t {0} : {1}'.format(k, v)))
            ftxIn.write(('{0}={1}\n'.format(k, v)))
        else: #not float, int or string; assume it's list
            print('\t \t type of', k, ' is ', type(v))
            sys.stdout.flush()
            v_string=''
            for i in range(len(v)):
                v_string+=str(v[i])+' '
            v_string = v_string[:-1]
            print(('\t \t {0} : {1}'.format(k, v_string)))
            ftxIn.write(('{0}={1}\n'.format(k, v_string)))
            del v_string
    ftxIn.close()
    print(' ')
    
    sys.stdout.flush()

    if logFile  is not None:
        if print_test:
            print('\t close log file ', outF.name)
            print('\t set sys.stdout back to: ', orig_stdout)
        sys.stdout.flush()
        sys.stdout = orig_stdout
        outF.close()



    
    return

if __name__ == '__main__':

    import shutil
    plasmaOut2ftxIn()
   

