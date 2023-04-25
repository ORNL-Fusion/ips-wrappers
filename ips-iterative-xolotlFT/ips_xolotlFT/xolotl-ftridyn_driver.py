#! /usr/bin/env python

from ipsframework import Component
import sys
import os
import os.path
import subprocess
import numpy
import shutil
from ips_xolotlFT.python_scripts_for_coupling import translate_xolotl_to_ftridyn
#import translate_ftridyn_to_xolotl_launch
#import get_yields_launch
from ips_xolotlFT.python_scripts_for_coupling import binTRIDYN
from ips_xolotlFT.python_scripts_for_coupling import param_handler
import traceback
from ips_xolotlFT.python_scripts_for_coupling import transferGrid
import pickle
from ips_xolotlFT.python_scripts_for_coupling import keepLastTS
from ips_xolotlFT.python_scripts_for_coupling import write_tridynDat
from ips_xolotlFT.python_scripts_for_coupling import handle_tempModel
from ips_xolotlFT.python_scripts_for_coupling import handle_gridModel
import inspect

class xolotlFtridynDriver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeStamp=0.0, **keywords):

        cwd = self.services.get_working_dir()

        print('  FT-X driver:init called ')
        print('\t with keywords: ',keywords)
        print(' ')
        
        print('\t output file of the FT-X workflow:')
        if 'LOG_FILE' in keywords:
            logFile=keywords['LOG_FILE']
            outFile=cwd+'/'+logFile
            print('\t log file defined in keywords: ')
            print('\t', outFile)
            outF=open(outFile , 'a')
            sys.stdout = outF
        else:
            try:
                self.LOG_FILE
                logFile = self.LOG_FILE
                outFile=logFile
                print('\t log file defined in config file', outFile)
                outF=open(outFile , 'a')
                sys.stdout = outF            
            except:
                print('\t No log file defined; using default sys.stdout')
                outfile=None

        print(' ')
        print('---------------------------')
        print('---------------------------')
        print('  FT-X driver:init  ')
        print('---------------------------')
        print('---------------------------')
        print(' ')
        
        try:
            if self.PRINT_TEST=='on':
                self.print_test=True
                print('Print TEST lines (devel / debug)')
            elif self.PRINT_TEST=='off':
                self.print_test=False
                print('Print only std output lines (no TEST lines)')
            else:
                self.print_test=False
                print('print_TEST values not valid')
                print('Print only std output lines (no TEST lines)')
        except Exception as e:
            print(e)
            print('print_TEST values not found')
            print('Print only std output lines (no TEST lines)')
            self.print_test=False
        print(' ')

        try:
            if self.PRINT_UQ=='on':
                self.print_uq=True
                print('Print UQ = on ; print ">>PR" lines')
            elif self.PRINT_UQ=='off':
                self.print_uq=False
                print('print UQ = off ; Do not print UQ lines')
            else:
                self.print_uq=False
                print('invalid print UQ value ; Do not print UQ lines')
        except Exception as e:
            print(e)
            print('print UQ not found ; Do not print UQ lines')
            self.print_uq=False
        print(' ')

        #test giving explicit wrapper path in modernized FTX workflow
        if self.print_test:
            try:
                self.SCRIPT
                if self.SCRIPT == "":
                    print('no explicit script path provided. use module loaded in environment')
                else:
                    print('using explicit path to wrapper')
                    print(self.SCRIPT)
            except Exception as e:
                print(e)
                print('no script variable defined. use module loaded in environment')
            print(' ')
            
        #stage input files
        print('staging input files {} '.format(self.INPUT_FILES))
        self.services.stage_input_files(self.INPUT_FILES)
        print('\t ...input files staged succesfully')
        print(('input directory for this simulation is {} \n'.format( self.INPUT_DIR)))


        #round the time to the decimal point
        if 'time_decimal' in keywords:
            self.time_decimal=int(keywords['time_decimal'])
        else:
            self.time_decimal=5
        print('round time to the ', self.time_decimal, 'th digit')
        
        plasma_state_file = self.STATE_FILES #to only stage/update what the driver needs ; self.services.get_config_param('PLASMA_STATE_FILES')
        plasma_state_list = plasma_state_file.split()
        for index in range(len(plasma_state_list)):
            open(plasma_state_list[index], 'a').close()
            if self.print_test:
                print('\t created: ', plasma_state_list[index])
        print(' ')
        #A MORE ELEGANT WAY --  FOR THE FUTURE
            #for file in plasma_state_list:
            #    open(file, 'a').close()
            
        self.services.update_state()
        self.services.stage_state()

        #### DRIVER PARAMETERS #####
        
        print('reading DRIVER parameters from FTX config file: \n')
        
        self.driver={}
        for k,v in self.DRIVER_INPUT_PARAMETERS.items():
            if param_handler.is_int(v):
                print(('\t integer input parameter {0} = {1}'.format(k, v)))
                self.driver[k]=int(v)
            elif param_handler.is_float(v):
                print(('\t float input parameter {0} = {1}'.format(k , v)))
                self.driver[k]=float(v)
            else:
                print(('\t other {0} input parameter {1} = {2}'.format(type(v), k , v)))
                self.driver[k]=v
        print('\n')
        sys.stdout.flush()

        ## TIME CHARACTERISTICS: START, STEP, END
        ## migth be given in config file or in separate (input) file
        self.time = {}

        #check what/if given in config file in its own section:
        try:
            self.TIME_PARAMETERS
            print('reading TIME parameters from FTX config file: \n')
            for k,v in self.TIME_PARAMETERS.items():
                if param_handler.is_int(v):
                    print(( '\t reading integer input parameter {0} = {1}'.format( k, v)))
                    self.time[k]=int(v)
                elif param_handler.is_float(v):
                    print(( '\t reading float input parameter {0} = {1}'.format( k, v)))
                    self.time[k]=float(v)
                elif param_handler.is_list(v):
                    values = []
                    for val in v.split(' '):
                        if param_handler.is_int(val):
                            values.append(int(val))
                        elif param_handler.is_float(val):
                            values.append(float(val))
                        else:
                            values.append(val)
                    print(('\t reading list input parameter {0} = {1}'.format( k, values)))
                    self.time[k]=values
                else:
                    print(('\t reading string input parameter {0} = {1} '.format( k, v)))
                    self.time[k]=v
            print('\n')
            if self.print_test:
                print('\t TEST: after reading parameters from config file time section, the time dictionary is:')
                print('\t', self.time)
            print(' ')
            sys.stdout.flush()
                
        except Exception as e:
            print(e)
            print('no [TIME_PARAMETERS] defined in config file')
            sys.stdout.flush()
            
            print('check if INIT, STEP, END TIME and other paramters (LPPS_TS, START_MODE) are defined in driver section of config file')
            if ('INIT_TIME' in self.driver):
                self.time['INIT_TIME']=self.driver['INIT_TIME']
            if ('END_TIME' in self.driver):
                self.time['END_TIME']=self.driver['END_TIME']
            if ('LOOP_TIME_STEP' in self.driver):
                self.time['LOOP_TIME_STEP']=self.driver['LOOP_TIME_STEP']
            if ('LOOP_N' in self.driver):
                self.time['LOOP_N']=self.driver['LOOP_N']
            if ('LOOP_TS_FACTOR' in self.driver):
                self.time['LOOP_TS_FACTOR']=self.driver['LOOP_TS_FACTOR']
            if ('LOOP_TS_NLOOPS' in self.driver):
                self.time['LOOP_TS_NLOOPS']=self.driver['LOOP_TS_NLOOPS']
            if ('START_MODE' in self.driver):
                self.time['START_MODE']=self.driver['START_MODE']

        if self.print_test:
            print('\t TEST: after reading parameters from config file, the time dictionary (in driver) is:')
            print('\t', self.time)
            print(' ')
        #Check if/what's given in input file
        #overwrite config file values
        try:
            self.TIME_FILE
            print('\t use input from time file defined in config file (overwirte config file values)')
            print('\t', self.INPUT_DIR+'/'+self.TIME_FILE)
            sys.stdout.flush()
            self.time_temp=param_handler.read(self.INPUT_DIR+'/'+self.TIME_FILE)            
            
            if self.print_test:
                print('\t TEST: after reading parameters from time file, the time_temp dictionary is:')
                print('\t', self.time_temp)
                print(' ')
                
            ## for testing purposes, a more thorough print line:
            ## use these three afterwards
            if self.print_test:
                for k,v, in self.time.items():
                    print(('\t {0} : {1}'.format(k, v)))
                print(' ')
            sys.stdout.flush()
            
            for k,v in self.time_temp.items():
                self.time[k]=v
                if param_handler.is_int(v):
                    print(( '\t reading integer input parameter {0} = {1}'.format( k, v)))
                elif param_handler.is_float(v):
                    print(( '\t reading float input parameter {0} = {1}'.format( k, v)))
                elif param_handler.is_list(v):
                    values = []
                    for val in v.split(' '):
                        if param_handler.is_int(val):
                            values.append(int(val))
                        elif param_handler.is_float(val):
                            values.append(float(val))
                        else:
                            values.append(val)
                    print(('\t reading list input parameter {0} = {1}'.format( k, values)))
                else:
                    print(('\t reading string input parameter {0} = {1} '.format( k, v)))

            if self.print_test:
                print('\t TEST: after copyting values from time_temp, the time dictionary is:')
                print('\t', self.time_temp)
            sys.stdout.flush()
            
        except Exception as e:
            print(e)
            print('no TIME_FILE defined. Use time-parameters defined in config file')
        print(' ')
        sys.stdout.flush()
            
        print('after reading from config file and file, time characteristics are:')
        for k,v, in self.time.items():
            print(('\t {0} : {1}'.format(k, v)))
        print(' ')
        print(('running FTX from t = {0} , to t = {1} , in steps of dt = {2}'.format(self.time['INIT_TIME'], self.time['END_TIME'], self.time['LOOP_TIME_STEP'])))
        print('\t with START_MODE = ', self.time['START_MODE'])
        
        self.driverMode=self.time['START_MODE']
        print('\n')
        
        #### FTRIDYN PARAMETERS #####

        print(' ')
        print('reading FTRIDYN parameters from ips config file: \n')

        sys.stdout.flush()
        self.ftridyn={}
        for k,v in self.FTRIDYN_INPUT_PARAMETERS.items():
            if param_handler.is_int(v):
                print(('\t integer input parameter {0} = {1}'.format(k, v)))
                self.ftridyn[k]=int(v)
            elif param_handler.is_float(v):
                print(('\t float input parameter {0} = {1}'.format(k, v)))
                self.ftridyn[k]=float(v)
            else:
                print(('\t other {0} input parameter {1} = {2}'.format(type(v), k, v)))
                self.ftridyn[k]=v
        print('\n')
        sys.stdout.flush()
        
        #### XOLOTL PARAMETERS ##### 

        print('\n')
        print('XOLOTL paramters: \n')

        #get dimension to read the correct input parameter template
        dim=int(self.XOLOTL_INPUT_PARAMETERS['dimensions'])
        #xolotl_param_template=self.XOLOTL_PARAM_TEMPLATE_PATH+'/paramXolotl_'+str(dim)+'D.txt'
        xolotl_param_template=self.INPUT_DIR+'/paramXolotl_'+str(dim)+'D.txt' 
        print(('\t reading Xolotl default parameters from {}'.format(xolotl_param_template)))

        sys.stdout.flush()
        self.xp = param_handler.xolotl_params()
        self.xp.read(xolotl_param_template)
        print(('\t running Xolotl in {} D'.format(dim)))
        print(' ')
        sys.stdout.flush()

        if (self.print_uq):
            print(f"PR: >>> I have read these Xolotl parameters from the template file {xolotl_param_template}:")
            for k, v in self.xp.parameters.items():
                print(f"PR: >>> \t{k} : {v}")
            print(f"PR: >>> I have read these Xolotl parameters from ips.ftx.config:")
            for k, v in self.XOLOTL_INPUT_PARAMETERS.items():
                print(f"PR: >>> \t{k} : {v}")
        
        #overwrite default Xolotl parameters that are specified in ips.config
        print('modify XOLOTL paramters with parameters read from the simulations config file: \n')

        sys.stdout.flush()
        for k,v in self.XOLOTL_INPUT_PARAMETERS.items():
            if param_handler.is_int(v):
                if k in self.xp.parameters:
                    print(('\t integer input parameter {0} = {1} with {2}'.format(k, self.xp.parameters[k] , v)))
                else:
                    print(('\t add integer input parameter {0} = {1} to xolotl parameters'.format(k , v)))
                self.xp.parameters[k]=int(v)
            elif param_handler.is_float(v):                
                if k in self.xp.parameters:
                    print(('\t float input parameter {0} = {1} with {2}'.format(k, self.xp.parameters[k] , v)))
                else:
                    print(('\t add float input parameter {0} = {1} to xolotl parameters'.format(k , v)))
                self.xp.parameters[k]=float(v)
            elif param_handler.is_list(v):
                values = []
                for val in v.split(' '):
                    if param_handler.is_int(val):
                        values.append(int(val))
                    elif param_handler.is_float(val):
                        values.append(float(val))
                    else:
                        values.append(val)
                print(('\t reading list input parameter {0} = {1} '.format(k, values)))
                self.xp.parameters[k]=values
            elif v=='False':
                if k in self.xp.parameters:
                    print(('\t input parameter {0} = {1} and thus will delete it from the xolotl parameters'.format(k, v)))
                    del self.xp.parameters[k]
                else:
                    print(('\t WARNING: input parameter {0} = {1} not present in xolotl parameters. Do nothing. '.format(k , v)))
            else:                
                print(('\t reading string input parameter {0} = {1}'.format(k, v )))
                self.xp.parameters[k]=v
        sys.stdout.flush()
            
        print('\n')
        print('replacing PETSC arguments from ips config file:\n')
        for k,v in self.XOLOTL_INPUT_PETSCARGS.items():
            if (v == 'False'):
                if (k in self.xp.parameters['petscArgs']) :
                    del self.xp.parameters['petscArgs'][k]
                    print('removed ', k, 'from petsc arguments')
                else:
                    print ('WARNING: could not remove ', k, ' from petsc arguments. It does not exists in template ')
            elif param_handler.is_int(v):
                print(('\t integer {0} = {1} with {2} '.format(k , self.xp.parameters['petscArgs'][k] , v)))
                self.xp.parameters['petscArgs'][k]=int(v)
            elif param_handler.is_float(v):
                print(('\t float {0} = {1} with {2}'.format(k , self.xp.parameters['petscArgs'][k] , v))) 
                self.xp.parameters['petscArgs'][k]=float(v)
            else:
                print(('\t other {0} argument {1} = {2} with {3}'.format(type(v), k , self.xp.parameters['petscArgs'][k] , v)))
                self.xp.parameters['petscArgs'][k]=v
        sys.stdout.flush()

        #make sure grid is in the right format:
        print(' ')
        print('check grid format:')
        if self.driver['xolotl_v']==1:
            if self.print_test:
                print('calling handle grid model from:')
                print(os.path.abspath(inspect.getfile(handle_gridModel.v1)))
            gridVal,rm_gridType, rm_gridParam=handle_gridModel.v1(xp_parameters=self.xp.parameters, print_test=self.print_test)
            self.xp.parameters['grid']=gridVal
            if rm_gridType:
                del self.xp.parameters['gridType']
            if rm_gridParam:
                del self.xp.parameters['gridParam']
        elif self.driver['xolotl_v']==2:
            if self.print_test:
                print('calling handle grid model from:')
                print(os.path.abspath(inspect.getfile(handle_gridModel.v2)))
            gridType,gridVal,rm_grid,rm_regularGrid = handle_gridModel.v2(xp_parameters=self.xp.parameters, print_test=self.print_test)
            self.xp.parameters['gridType']=gridType
            self.xp.parameters['gridParam']=gridVal
            if rm_grid:
                del self.xp.parameters['grid']
            if rm_regularGrid:
                del self.xp.parameters['regularGrid']

        #CONTROL WHICH PROCESSES ARE ON:
        #delete Xolotl processes that are specified as false in ips.config
        print(' ')
        xp_processes={}
        for k,v in self.XOLOTL_INPUT_PROCESSES.items(): #specified in ips.config
            xp_processes[k]=v

        processList=self.xp.parameters['process'] #all processes, specified in param template
        processString=''

        for p in range(len(processList)):
            key=processList[p]
            if (key in xp_processes) and (xp_processes[key]=='false'):
                print(('delete {} from Xolotl params, as it was set to false in ips.config \n'.format(key))) 
            else:
                processString+=key+' '

        print('processes included in Xolotl are:')
        print('\t {}\n'.format(processString.strip()))
        self.xp.parameters['process']=processString.strip()

        #performance handler can only be "dummy" or "os" now:
        if not (self.xp.parameters['perfHandler'] == "os" or self.xp.parameters['perfHandler'] == "dummy"):
            print ('invalid performance handler in Xolotl.')
            print(' \t Change from ', self.xp.parameters['perfHandler'], 'to default, "os" ')
            self.xp.parameters['perfHandler'] = "os"
        sys.stdout.flush()
        
        ### PLASMA RELATED PARAMETERS ###
        ## code agnostic (gitr, solps... all handled together)
        
        #keep track if plasma input exists (file or config file) with string: F = false, T = true
        plasmaExists=False
        self.plasma = {}
        try:
            self.PLASMA_INPUT_PARAMETERS
            plasmaExists=True
            print('Read PLASMA inputs defined in config file.')
            sys.stdout.flush()
            
            for k,v in self.PLASMA_INPUT_PARAMETERS.items():
                if param_handler.is_int(v):
                    print(( '\t reading integer PLASMA input parameter {0} = {1}'.format( k, v)))
                    self.plasma[k]=int(v)
                elif param_handler.is_float(v):
                    print(( '\t reading float PLASMA input parameter {0} = {1}'.format( k, v)))
                    self.plasma[k]=float(v)
                elif param_handler.is_list(v):
                    values = []
                    for val in v.split(' '):
                        if param_handler.is_int(val):
                            values.append(int(val))
                        elif param_handler.is_float(val):
                            values.append(float(val))
                        else:
                            values.append(val)
                    print(('\t reading list PLASMA input parameter {0} = {1}'.format( k, values)))
                    self.plasma[k]=values                
                else:
                    print(('\t reading string PLASMA input parameter {0} = {1} '.format( k, v))) 
                    self.plasma[k]=v
            sys.stdout.flush()
                    
        except Exception as e:
            print(e)
            print('no PLASMA_INPUT_PARAMETERS defined in config file. ')
            sys.stdout.flush()
        print(' ')
        
        try:
            self.PLASMA_OUTPUT_FILE
            plasmaExists=True
            print('use input from PLASMA output file (as defined in config file)')
            print('\t read PLASMA parameters (from file) :')
            print('\t \t {}'.format(self.INPUT_DIR+'/'+self.PLASMA_OUTPUT_FILE))
            print('\t these inputs will overwrite values given in config file ')
            self.plasma_temp=param_handler.read(self.INPUT_DIR+'/'+self.PLASMA_OUTPUT_FILE)

            for k,v, in self.plasma_temp.items():
                print(('\t \t {0} : {1}'.format(k, v)))
                self.plasma[k]=v
            print(' ')
            sys.stdout.flush()

            #also take the opportunity to copy it here and update the plasma state:
            print('\t copy ', self.PLASMA_OUTPUT_FILE, ' to cwd, so it can be updated in plasma state')
            shutil.copyfile(self.INPUT_DIR+'/'+self.PLASMA_OUTPUT_FILE, self.PLASMA_OUTPUT_FILE)
            try:
                print('\t copy ', self.PLASMA_OUTPUT_ORIG, ' to cwd, so it can be updated in plasma state')
                shutil.copyfile(self.INPUT_DIR+'/'+self.PLASMA_OUTPUT_ORIG, self.PLASMA_OUTPUT_ORIG)
            except Exception as e:
                print('\t', e)
            sys.stdout.flush()
            self.services.update_state()
            
        except Exception as e:
            print(e)
            print('no PLASMA_OUTPUT_FILE defined in config file. ')
            sys.stdout.flush()

        print(' ')
            
        if self.print_test:
            print('TEST: final PLASMA dictionary is: ')
            for k,v in self.plasma.items():
                print('\t ', k ,'=', v)
            print(' ')
        
        if plasmaExists:
            if 'plasmaOutputDir' in self.plasma:
                print('\t PLASMA output will be read from {} \n'.format(self.plasma['plasmaOutputDir']+'_'+prj))

            ## OVERWRITE FTX PARAMETERS WITH INPUTS FROM PLASMA 

            print('overwrite FTX parameters with inputs from plasma:')
            print(' ')
            
            if 'flux' in self.plasma:
                self.xp.parameters['flux']=self.plasma['flux']*1.0e-18
                print(('\t replaced flux in Xolotl (in ion/nm2s = 1e-18 ion/m2s) by value given by PLASMA {} (in ion/m2s) \n'.format(self.plasma['flux'])))
            else:
                print(('\t no flux specified in PLASMA, so using values in Xolotl {} \n'.format(self.xp.parameters['flux'])))
            sys.stdout.flush()
            
            ## parameters for the temperature model depend on the Xolotl executable
            ## v1 (e.g., master) uses a single line: heat = [value], startTemp = [value], etc.
            ## v2 (e.g., tempGrid) uses two lines:  tempHandler = heat / ...
            ##                                      heat / tempParam = [value]
            ## this ends up being quite long, so better in its own file:
            ## Xolotl v1 needs one line, v2 needs 2 lines --> handle_tempModel contains 2 functions

            print('\t temperature model handled by handle_tempModel: ')
            if (self.driver['xolotl_v']==1):
                mod,val,rm_startTemp = handle_tempModel.v1(xp_parameters=self.xp.parameters,plasma=self.plasma, print_test=self.print_test) 
                sys.stdout.flush()
                self.xp.parameters[mod]=val
                print('\t ... write_tempModel succesfully returned: for Xolotl v1, temperature model: ', mod, val)
            elif(self.driver['xolotl_v']==2):
                mod,val,rm_startTemp = handle_tempModel.v2(xp_parameters=self.xp.parameters,plasma=self.plasma, print_test=self.print_test)
                sys.stdout.flush()
                self.xp.parameters['tempHandler']=mod
                self.xp.parameters['tempParam']=val
                print('\t ... handle_tempModel succesfully returned: for Xolotl v2, temperature model: ', mod, val)
            else:
                print('\t WARNING: Xolotl version not recognized: ')
                print( '\t \t xolotl_v = ', self.driver['xolotl_v'])
                mod='startTemp'
                val=300
                rm_startTemp=False
                print('\t \t assume standard values: ', mod, ' = ', val)
            if rm_startTemp:
                del self.xp.parameters['startTemp']
                print('\t \t and delete startTemp')
            print(' ')

        else:
            print('no PLASMA input file or paramters defined in config file')
        sys.stdout.flush()

        print(' ')
        
        #initialize lists that will contain values for each species
        #input parameters
        self.energyIn=[]
        self.inAngle=[] #values in input file
        self.angleIn=[] #values used in simulation: read from ips.config or from angular distrib file
        self.weightAngle=[]
        #self.fluxFraction=[]
        #other parameters
        self.spYield=[]
        self.rYield=[]
        self.spYieldMode=[]        
        self.rYieldMode=[]
        self.maxRangeXolotl=[]
        #files
        #given as string of 4 (name for W, He, D and T) -> split into list
        self.ft_input_file=self.FT_INPUT_FILE.split()        
        self.ft_energy_input_file=self.FT_ENERGY_INPUT_FILE.split()
        self.ftx_lay_file=self.FTX_LAY_FILE.split()
        self.ft_output_file=self.FT_OUTPUT_FILE.split()
        self.ft_output_prj_file=self.FT_OUTPUT_PRJ_FILE.split()
        #others
        self.FT_energy_file_name=[]
        self.angleFile=[]
        self.aWeightFile=[]
        self.eadist_output_path=[]
        self.eadist_output_file=[]

        
        try:
            self.PLASMA_OUTPUT_FILE
            print('read inputEnergy, inputAngle, plasmaSpecies and fluxFraction from PLASMA input file')
            
            if 'inputEnergy' in self.plasma:
                inputEnergy=self.plasma['inputEnergy']
            if 'inputAngle' in self.plasma:
                inputAngle=self.plasma['inputAngle']
            if 'plasmaSpecies' in self.plasma:
                self.plasmaSpecies=self.plasma['plasmaSpecies']
            if 'fluxFraction' in self.plasma:
                self.fluxFraction=self.plasma['fluxFraction']

        except Exception as e2:
                print(e2)
                print('no input from PLASMA')
                print('expect inputEnergy, inputAngle, plasmaSpecies and fluxFraction from config file (in FTRIDYN dictionary)')
                print('will crash if not provided')
                inputEnergy=self.ftridyn['inputEnergy']
                inputAngle=self.ftridyn['inputAngle']
                self.plasmaSpecies=self.ftridyn['plasmaSpecies']
                self.fluxFraction=self.ftridyn['fluxFraction']
        
        inputSpYield=self.ftridyn['inputSpYield'].split(' ')
        inputRYield=self.ftridyn['inputRYield'].split(' ')
        #inputFluxFraction=self.GITR_INPUT_PARAMETERS['inputFluxFraction'].split(' ')        

        for i in range(len(self.plasmaSpecies)): #self.plasmaSpecies.iteritems():
            prj=self.plasmaSpecies[i]

            self.energyIn.append(float(inputEnergy[i]))
            self.inAngle.append(float(inputAngle[i]))
            self.spYield.append(float(inputSpYield[i]))
            self.rYield.append(float(inputRYield[i]))
            #self.fluxFraction.append(float(inputFluxFraction[i]))
            print(('\t index {0}, species {1} '.format( i, self.plasmaSpecies[i])))
            print(('\t energy {0}, angle {1} '.format( self.energyIn[i] ,self.inAngle[i])))
            print(('\t spYield {0} , rYield {1} , fluxFraction {2} \n'.format(self.spYield[i], self.rYield[i], self.fluxFraction[i] )))

            if self.inAngle[i] < 0 :
                try:
                    self.PLASMA_OUTPUT_FILE
                    print('angular distributions given by PLASMA')
                    plasma_output_dir='plasmaOutputDir'+'_'+prj
                    print('\t PLASMA output of angular distributions will be read from')
                    print('\t \t {}'.format(self.plasma[plasma_output_dir]))
                    self.angleFile.append(self.plasma[plasma_output_dir].strip()+'/'+self.PLASMA_ANGLE_FILE.strip()) #self.gitr['gitrOutputDir'].strip()
                    self.aWeightFile.append(self.plasma[plasma_output_dir].strip()+'/'+self.PLASMA_AWEIGHT_FILE.strip()) #self.gitr['gitrOutputDir'].strip()

                except Exception as e2:
                    print(e2)
                    print('no input from PLASMA')
                
                print('\t reading angles and weights for {0} from: '.format(prj))
                print('\t \t {0}'.format(self.angleFile[i]))
                print('\t \t {0}'.format(self.aWeightFile[i]))
                a = numpy.loadtxt(self.angleFile[i], usecols = (0,) , unpack=True)
                w = numpy.loadtxt(self.aWeightFile[i], usecols = (0,) , unpack=True)
                self.angleIn.append(a)
                self.weightAngle.append(w)
            else:
                self.angleIn.append([self.inAngle[i]])
                self.weightAngle.append([1.0])
                self.angleFile.append('')
                self.aWeightFile.append('')
                print(('\t {} angle value as defined by user \n'.format(prj)))
            print(' ')

            #AND MAYBE SOMETHING SIMILAR WITH ENERGIES?

            if self.spYield[i]<0:
                self.spYieldMode.append('calculate')                
            else:
                self.spYieldMode.append('fixed')

            if self.rYield[i]<0:
                self.rYieldMode.append('calculate')
            else:
                self.rYieldMode.append('fixed')

            #FTRIDYN FILES
            #prepare input files; i.e., those transferred from FT init (generateInput) to FT step (run code)
            #leave 'others' empty for a pure FT run
            try:
                self.PLASMA_OUTPUT_FILE                
                if self.energyIn[i] < 0:
                    plasma_output_dir='plasmaOutputDir'+'_'+prj
                    print('\t PLASMA output of energy distributions will be read from:')
                    print ('\t \t {}'.format(self.plasma[plasma_output_dir]))
                    self.FT_energy_file_name.append(self.ft_energy_input_file[i]) #"He_W0001.ED1"
                    #where all the energy distribution files are located
                    self.eadist_output_path.append(self.plasma[plasma_output_dir].strip())
                    try:
                        self.PLASMA_EADIST_FILE 
                        self.eadist_output_file.append(self.PLASMA_EADIST_FILE)
                        print('\t using energy distribution file format given in PLASMA eadist file: ', self.PLASMA_EADIST_FILE)
                    except : #default file format ['dist','.dat']
                        self.eadist_output_file.append(['dist','.dat'])
                        print("\t using default energy distribution file format, ['dist','.dat']")
                else:       
                    self.FT_energy_file_name.append('')         
                    self.eadist_output_path.append('')
                    self.eadist_output_file.append([' ',' '])
            except Exception as e2:
                print(e2)
                print('no input from PLASMA')

            print(' ')
            #initialize maxRangeXolotl list
            self.maxRangeXolotl.append(0.0)


        #stage initial network File (INIT mode) OR restart files (RESTART mode)
        if (self.driverMode=='INIT'):
            #save copy to last_TRIDYN.dat (and networkFile, below), as they'll be overwritten later
            if (os.path.exists(self.INPUT_DIR+'/last_TRIDYN.dat')):
                shutil.copyfile(self.INPUT_DIR+'/last_TRIDYN.dat',self.INPUT_DIR+'/last_TRIDYN_init.dat')
            if (os.path.exists(self.INPUT_DIR+'/'+self.NETWORK_FILE)):
                #save copy of initial network file, because it'll be overwritten in restart
                shutil.copyfile(self.INPUT_DIR+'/'+self.NETWORK_FILE, self.INPUT_DIR+'/networkFile_init.h5') 
                if ('networkFile' in self.xp.parameters):
                    print(('\t INIT mode: stage initial network file {}\n'.format(self.NETWORK_FILE)))
                    self.services.stage_input_files(self.NETWORK_FILE)
                    print('\t \t ...initial network file staged succesfully {}\n')
            else:
                print('\t WARNING: INIT mode network file:')
                print('\t \t either networkFile not defined in config, ')
                print('\t \t or could not find initial network file in input file directory {}\n'.format(self.INPUT_DIR))
                print('\t \t SKIP staging  {}\n'.format(self.NETWORK_FILE))
                

        elif (self.driverMode=='RESTART'):
            restart_files = self.services.get_config_param('RESTART_FILES') 
            print(('\t RESTART mode: stage restart files {} \n'.format(restart_files)))
            restart_list = restart_files.split()
            for index in range(len(restart_list)): 
                print('\t staging input file ', restart_list[index]) 
                self.services.stage_input_files(restart_list[index])
            print('\t \t ...restart files staged succesfully')

        sys.stdout.flush()
        self.services.update_state()


    def step(self, timeStamp=0.0,**keywords):

        cwd = self.services.get_working_dir()

        print('  FT-X driver:step called')
        print('\t with keywords: ',keywords)
        
        print('\t output file of the FT-X workflow:')
        if 'LOG_FILE' in keywords:
            logFile=keywords['LOG_FILE']
            outFile=cwd+'/'+logFile 
            print('\t \t log file defined in keywords: ')
            print('\t \t ', outFile)
            outF=open(outFile , 'a')
            sys.stdout = outF
        else:
            try:
                self.LOG_FILE
                logFile = self.LOG_FILE
                outFile=logFile
                print('\t \t log file defined in config file', outFile)
                outF=open(outFile , 'a')
                sys.stdout = outF
            except:
                print('\t \t No log file defined; using default sys.stdout')
                outFile=None

        print('\n')
        print('---------------------------')
        print('---------------------------')
        print('  FT-X driver:step  ')
        print('---------------------------')
        print('---------------------------')
        print(' ')


        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')
        
        self.services.stage_state() 

        #round the time to the decimal point
        #if 'time_decimal' in keywords:
        #    self.time_decimal=int(keywords['time_decimal'])
        #else:
        #    self.time_decimal=5
        #print('round time to the ', self.time_decimal, 'th digit')
        #print(' ')
        #sys.stdout.flush()
        
        #check that loop doesnt go over the end time
        ##if (self.print_uq):
        #print("\t Rounding time step here:")
        #time=self.time['INIT_TIME']
        #rounded_time = round(time, self.time_decimal)
        ##if (self.print_uq):
        #print(f"\t Given time step is {time}, rounding to {rounded_time}")
        #time = rounded_time
        ##if (self.print_uq):
        #print("\t Rounding end time here:")
        #end_time=self.time['END_TIME']
        #rounded_end_time = round(end_time, self.time_decimal)
        ##if (self.print_uq):
        #print(f"\t end time is {end_time}, rounding to {rounded_end_time}")
        #end_time = rounded_end_time

        #if time >= end_time:
        #    print('init time ', time , ' >= end time', end_time)
        #    print("EXIT SIMULATION")
        #    sys.stdout.flush()
        #    return

        #estimated_end_time=round(time+self.time['LOOP_TIME_STEP'], self.time_decimal)
        #print('TEST: time+self.time[LOOP_TIME_STEP] = ', time+self.time['LOOP_TIME_STEP'] , ' ; estimated_end_time (rounded) = ', estimated_end_time)
        #print('\t before starting time-loop, checked that time step given in config file is not longer than needed to reach the end of the simulation')
        #if estimated_end_time > end_time: #time+self.time['LOOP_TIME_STEP']>end_time:
        #    self.time['LOOP_TIME_STEP']=round(end_time-time,self.time_decimal)
        #    self.xp.parameters['petscArgs']['-start_stop']=end_time/10.0
        #    if (self.print_uq):
        #        print(f"PR: >>> Updated -start_stop to {self.xp.parameters['petscArgs']['-start_stop']:.4f}")
        #    print('\t WARNING: time step given in config file longer than needed for last loop ')
        #    print('\t \t before starting time-loop, adapt driver time step to {} to reach exactly endTime (rounded to {}th decimal) '.format(self.time['LOOP_TIME_STEP'], self.time_decimal))
        #    self.xp.parameters['petscArgs']['-start_stop']=self.time['LOOP_TIME_STEP']/10.0
        #    print(('\t \t accordingly, Xolotls data is saved every (start_stop) = {} '.format( self.xp.parameters['petscArgs']['-start_stop'])))
        #else:
        #    print('\t time-step is shorter than needed. Continue with it')
        #    self.xp.parameters['petscArgs']['-start_stop']=(time+self.time['LOOP_TIME_STEP'])/10.0
        #    if (self.print_uq):
        #        print(f"PR: >>> Updated -start_stop to {self.xp.parameters['petscArgs']['-start_stop']:.4f}")

        #print(' ')
        #sys.stdout.flush()

        time=self.time['INIT_TIME']
        
        print('before starting time-loop, checked that time-parameters:')

        #round the time to the decimal point
        if 'time_decimal' in keywords:
            self.time_decimal=int(keywords['time_decimal'])
        else:
            self.time_decimal=5
        print('round time to the ', self.time_decimal, 'th digit')
        print(' ')
        sys.stdout.flush()

        
        print("\t Round end time:")
        end_time=self.time['END_TIME']
        rounded_end_time = round(end_time, self.time_decimal)
        print(f"\t end time is {end_time}, rounding to {rounded_end_time}")
        end_time = rounded_end_time
        self.time['END_TIME']=rounded_end_time
        print('\t assign rounded value back : self.time[END_TIME] = ', self.time['END_TIME'])
        
        #for time in numpy.arange(self.initTime,self.endTime,self.timeStep):
        while time<end_time: #self.driver['END_TIME']:

            self.services.stage_state()
            print(('driver time (in loop) {}'.format(time)))
            self.services.update_state()

            #round the time to the decimal point
            print("\t Round time step:")
            rounded_time = round(time, self.time_decimal)
            print(f"\t Given time is {time}, rounding to {rounded_time}")
            time = rounded_time
            
            print('\t Round loop time step: ')
            rounded_loop_time_step = round(self.time['LOOP_TIME_STEP'], self.time_decimal)
            print(f"\t loop time step is {self.time['LOOP_TIME_STEP']}, rounding to {rounded_loop_time_step}")
            
            print('\t Round start_stop: ')
            rounded_start_stop = round(self.xp.parameters['petscArgs']['-start_stop'], self.time_decimal)
            print(f"\t loop time step is {self.xp.parameters['petscArgs']['-start_stop']}, rounding to {rounded_start_stop}")
            
            print('Check that step given in config file (rounded) is not longer than needed to reach the end of the simulation')
            if time >= end_time:
                print('time ', time , ' >= end time', end_time)
                print("EXIT SIMULATION")
                sys.stdout.flush()
                return

            #TO-DO: once it's working, add print_test
            estimated_end_time=round(rounded_time+rounded_loop_time_step, self.time_decimal)
            print('TEST: rounded_time+rounded_loop_time_step = ', time+rounded_loop_time_step , ' ; estimated_end_time (rounded) = ', estimated_end_time)
            print('\t before starting time-loop, checked that time step given in config file (rounded) is not longer than needed to reach the end of the simulation')
            if estimated_end_time > end_time: #time+self.time['LOOP_TIME_STEP']>end_time:
                rounded_loop_time_step=round(end_time-time,self.time_decimal)
                rounded_start_stop=round(rounded_loop_time_step/10, self.time_decimal) #end_time/10.0
                print('\t WARNING: time step given in config file longer than needed for last loop ')
                print('\t          adapt driver time step to {} to reach exactly endTime (rounded) '.format(rounded_loop_time_step))
                print('\t          and save Xolotls every (start_stop) = {} (rounded)'.format(rounded_start_stop))
            else:
                print('\t time-step is shorter than or exactly as needed to reach end time. Continue with:')
                print('\t rounded time step = ', rounded_loop_time_step)                
                print('\t rounded start_stop = ', rounded_start_stop)

            final_time = time+rounded_loop_time_step #self.time['LOOP_TIME_STEP']
            rounded_final_time = round(final_time, self.time_decimal)
            print('\t final time of this loop (rounded) is = ', rounded_final_time)
            
            #copy rounded values to time dictionary
            print('assign rounded values back to time dictionary:')
            self.time['LOOP_TIME_STEP']=rounded_loop_time_step
            self.xp.parameters['petscArgs']['-ts_final_time']=rounded_final_time
            self.xp.parameters['petscArgs']['-start_stop']=rounded_start_stop
            print('\t self.time[LOOP_TIME_STEP] = ', self.time['LOOP_TIME_STEP'])
            print('\t self.xp.parameters[petscArgs][-ts_final_time] = ', self.xp.parameters['petscArgs']['-ts_final_time'])
            print('\t self.xp.parameters[petscArgs][-start_stop] = ', self.xp.parameters['petscArgs']['-start_stop'])
            print(' ')
            sys.stdout.flush()
            
            #keep all files to be saved (not plasma state) in folder with time stamp
            timeFolder='t'+str(time)
            if not os.path.exists(timeFolder):
                os.makedirs(timeFolder)
            print(('\t output of this time-loop will be saved in {} \n'.format(timeFolder)))

            self.time['LOOP_N']+=1
            self.collapsedLoops=0 #reset 
            self.xolotlExitStatus='collapsed'

            print('set xolotl exit status back to collapsed \n')  
            sys.stdout.flush()

            ###################################### 
            ############## run FTridyn ############
            #### for each (Tg,Prj) combination ####
            ###################################### 

            print('\n')
            print('-----------')
            print('F-TRIDYN:  ')
            print('-----------\n')
            # A) GET INPUT THAT MIGHT CHANGE EVERT LOOP READY

            #determine parameters related to init/restart
            iW=self.plasmaSpecies.index('W')
            if (self.driverMode == 'INIT'):
                print('init mode yes')
                self.ftridyn['iQ0']=0

                targetList=[]
                targetList.append(self.plasmaSpecies[iW]) #only W in the first loop
                for i in range(1,4):
                    targetList.append('') #leave empty

                self.ftridyn['nDataPts'] = 100 #same as default value in generate_ftridyn_input
                if (self.ftridyn['totalDepth']==0.0):
                    self.ftridyn['nTT']=self.ftridyn['initialTotalDepth']
                else:
                    self.ftridyn['nTT']=self.ftridyn['totalDepth']
                
            else:
                print('init mode no')
                sys.stdout.flush()
                self.ftridyn['iQ0']=-1
                if self.print_test:
                    print('\t call xolotlToLay with args:')
                    print('\t \t totalDepth=',self.ftridyn['totalDepth'])
                    print('\t \t logFile=',outFile)
                    print('\t \t print_test=',self.print_test)
                    sys.stdout.flush()
                self.ftridyn['nDataPts'] = translate_xolotl_to_ftridyn.xolotlToLay(totalDepth=self.ftridyn['totalDepth'],logFile=outFile,print_test=self.print_test)
                #prepare target strings for F-Tridyn:
                #we always have W in the substrate  (tg1); others are optional; mixed material composition given by LAY file
                targetList=[]
                targetList.append(self.plasmaSpecies[iW])
                for prj in ['He', 'D', 'T']: #generate_ftridyn_input expects max 4 target species; declare all, even if empty
                    if prj in self.plasmaSpecies:
                        i=self.plasmaSpecies.index(prj)
                        if self.fluxFraction[i]>0.0: #species exists and fraction > 0
                            targetList.append(prj)
                        else:                        
                            targetList.append('') #leave empty
                    else:
                        targetList.append('') #leave empty
                        
                print(('\t passing to F-Tridyn the list of targets t{}'.format(targetList)))

                #Xolotl only outputs He_W0001.LAY; but it's same substrate composition for running all projectiles
                iHe=self.plasmaSpecies.index('He')
                for i in range(len(self.plasmaSpecies)): #self.plasmaSpecies.iteritems():
                    prj=self.plasmaSpecies[i]
                    if prj!='He': #do not self-copy
                        print(('\t copy {0} as {1} '.format( self.ftx_lay_file[iHe], self.ftx_lay_file[i]))) 
                        shutil.copyfile(self.ftx_lay_file[iHe],self.ftx_lay_file[i])
                
                if (self.ftridyn['totalDepth']==0.0):
                    print('\t Totaldepth from last_TRIDYN.dat') 
                    self.ftridyn['nTT']=10*numpy.max(numpy.loadtxt('last_TRIDYN.dat')[:,0])
                else:
                    print(('\t totalDepth fixed to {}'.format(self.ftridyn['totalDepth'])))
                    self.ftridyn['nTT']=self.ftridyn['totalDepth']
            print(' ')
            sys.stdout.flush()
            self.services.update_state()

            # B) RUN FTRIDYN

            print('call F-TRIDYN:')
            print(' ')
            for i in range(len(self.plasmaSpecies)): #self.plasmaSpecies.iteritems():
                prj=self.plasmaSpecies[i]
                maxDepth=[]
                if (self.fluxFraction[i] > 0.0) and not (all(k==0 for k in self.weightAngle[i])):
                    print(('\t running F-Tridyn for {0} with flux fraction = {1}'.format(prj, self.fluxFraction[i])))
                    print(('\t and not all angle weights are zero; max angleWeight is {}\n'.format(max(self.weightAngle[i]))))
                    sys.stdout.flush()
                    
                    #component/method calls now include arguments (variables)
                    self.services.call(ftridyn, 'init', timeStamp, dTime=time, fPrj=prj, fTargetList=targetList, ftParameters=self.ftridyn , fEnergyIn=self.energyIn[i], fAngleIn=self.angleIn[i], fWeightAngle=self.weightAngle[i], ft_folder=self.FT_OUTPUT_FOLDER, input_file=self.ft_input_file[i], otherInFiles=[self.FT_SURFACE_FILE,self.ftx_lay_file[i]], energy_file_name=self.FT_energy_file_name[i], orig_energy_files_path=self.eadist_output_path[i], orig_energy_files_pattern=self.eadist_output_file[i], output_file=outFile, print_test=self.print_test)
                    sys.stdout.flush()
                    
                    self.services.call(ftridyn, 'step', timeStamp, ftParameters=self.ftridyn, fEnergyIn=self.energyIn[i], fAngleIn=self.angleIn[i], fWeightAngle=self.weightAngle[i], fPrj=prj, output_file=outFile, print_test=self.print_test)
                    sys.stdout.flush()
                    self.services.stage_state()


                    # C) POSTPROCESSING OF prj -> W

                    # 1) access ft folder (read file containing path to FT output directory, as outputPath=...):
                    ft_dic=param_handler.read(self.FT_OUTPUT_PWD_PATH)
                    for key,value in ft_dic.items():
                        self.ftridyn[key] = value

                    print(' ')
                    print('F-TRIDYN run COMPLETED for ', prj)
                    print('\t from {0}:'.format(self.FT_OUTPUT_PWD_PATH))
                    print('\t \t the path to output of FTRIDYN is')
                    print('\t \t {0}'.format(self.ftridyn['outputPath']))
                    print(' ')

                    print('analyze and format output for ', prj, ':\n')
                    
                    #2) #get maximum projectile range to ensure bins are added correctly in 'translate_ftridyn_to_xolotl'

                    ft_output_prj_file=self.ft_output_prj_file[i]
                    ft_output_file=self.ft_output_file[i]
                    angleFolder=self.ftridyn['outputPath']+'/'+self.FT_OUTPUT_FOLDER+'/ANGLE'

                    for j in range(len(self.angleIn[i])):
                        if (self.weightAngle[i][j] > 0.0):
                            filePrj=angleFolder+str(self.angleIn[i][j])+'/'+ft_output_prj_file
                            num_lines_prj = sum(1 for line in open(filePrj))
                            if num_lines_prj == 0:
                                print('\t WARNING: no ions were implanted at this angle (prj file empty)')
                            elif num_lines_prj == 1:
                                depth, bla=numpy.loadtxt(filePrj, usecols = (2,3) , unpack=True)
                                maxDepth.append(depth)
                            elif num_lines_prj > 1:
                                depth, bla=numpy.loadtxt(filePrj, usecols = (2,3) , unpack=True)
                                maxDepth.append(max(depth))
                    
                    if len(maxDepth)>0:
                        print("something was implanted:")
                        maxRange=max(maxDepth)
                        self.maxRangeXolotl[i]=maxRange/10.0 #range in nm for Xolotl 
                        print(('\t maximum projectile range for {} is {} [A]'.format(prj, maxRange)))
                        print(' ')
                        #ft_output_file=self.ft_output_file[i]
                        #get implantation profile
                        #pass values as dictionary
                        ft_implProfiles_dictionary={}
                        ft_implProfiles_dictionary['ftridynOnePrjOutput']=ft_output_prj_file
                        ft_implProfiles_dictionary['ftridynFolder']=angleFolder
                        ft_implProfiles_dictionary['angle']=self.angleIn[i]
                        ft_implProfiles_dictionary['weightAngle']=self.weightAngle[i]
                        ft_implProfiles_dictionary['prjRange']=maxRange
                        ft_implProfiles_dictionary['print_test']=self.print_test

                        #different grid keywords depending on Xolotl version:
                        if (self.driver['xolotl_v']==1): #'grid' in self.xp.parameters:
                            ft_implProfiles_dictionary['nBins']=self.xp.parameters['grid'][0]
                        elif (self.driver['xolotl_v']==2): #'gridParam' in  self.xp.parameters:
                            ft_implProfiles_dictionary['nBins']=self.xp.parameters['gridParam'][0]
                        else:
                            ft_implProfiles_dictionary['nBins']=200
                            print('\t WARNING: xolotl_v not recognized; assume nBins=200')
                                
                        ft_implProfiles_dictionary['logFile']=outFile


                        pkl_impl_file='ft_implProfiles.pkl' #cwd+'/ft_implProfiles.pkl' ; rm cwd from paths
                        pickle.dump(ft_implProfiles_dictionary, open(pkl_impl_file, "wb" ) )

                        print('\t translate_ft2xol:')
                        sys.stdout.flush()
                        try:
                            self.TRANSLATE_FT2XOL
                            ft_implProfile_script=self.TRANSLATE_FT2XOL
                            print('\t Launch user-defined python script :  ')
                            print('\t', ft_implProfile_script)
                            
                        except: #DEFAULT PATH: 
                            ft_implProfile_script = 'translate_ftridyn_to_xolotl.py'
                            print('\t Launch default python script: ')
                            print('\t ',ft_implProfile_script)

                        sys.stdout.flush()

                        task_id_impl = self.services.launch_task(1,self.services.get_working_dir(),
                                                                 ft_implProfile_script, logfile='tridynPlotting.log')
                        ret_val_impl = self.services.wait_task(task_id_impl)

                    else: #if len(maxDepth)==0
                        print("\t nothing was implanted. maxRange and profile = 0 ")
                        maxRange=0.0
                        self.maxRangeXolotl[i]=0.0
                        outputFTFile=open(self.FT_OUTPUT_PROFILE_TEMP, "w")
                        outputFTFile.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ")
                        outputFTFile.close()
                        sys.stdout.flush()
                        #print "END OF THIS SIMULATION" & return

                    print(' ')
                    sys.stdout.flush()
                    
                    #3) get the sputtering yield (or use fixed value)                   
                    print('get sputtering and reflection yields:')
                    #pass values as dictionary
                    ft_getYields_dictionary={}
                    ft_getYields_dictionary['ftridynOneOutOutput']=ft_output_file
                    ft_getYields_dictionary['ftridynFolder']=angleFolder
                    ft_getYields_dictionary['angle']=self.angleIn[i]
                    ft_getYields_dictionary['weightAngle']=self.weightAngle[i]
                    ft_getYields_dictionary['fNImpacts']=self.ftridyn['nImpacts']
                    ft_getYields_dictionary['logFile']=outFile
                    ft_getYields_dictionary['print_test']=self.print_test
                    
                    pkl_gy_file='ft_getYields.pkl' # cwd+'/ft_getYields.pkl' ; rm cwd from paths
                    #we need to close this file later (to open it again and read yields), so use alternative to
                    #pickle.dump(ft_getYields_dictionary, open(pkl_gy_file, "wb" ) )
                    with open(pkl_gy_file, "wb") as pf:
                        pickle.dump(ft_getYields_dictionary, pf)
                    pf.close()
                    sys.stdout.flush()
                    
                    try:
                        self.GET_YIELDS
                        ft_getYields_script=self.GET_YIELDS
                        print('\t Launch user-defined python script : ')
                        print('\t ', ft_getYields_script)
                        
                    except: #DEFAULT PATH:
                        ft_getYields_script = 'get_yields.py'
                        print('\t Launch default python script ')
                        print('\t ', ft_getYields_script)

                    sys.stdout.flush()

                    task_id_gy = self.services.launch_task(1,self.services.get_working_dir(),
                                                           ft_getYields_script, logfile='get_yields.log')
                    ret_val_gy = self.services.wait_task(task_id_gy)

                    if os.path.exists(pkl_gy_file):
                        with open(pkl_gy_file, "rb") as pf:
                            getYields_dic = pickle.load(pf)
                            if self.print_test:
                                print('\t \t TEST: get_yields function returned dictionary:')
                                print('\t \t', getYields_dic)
                                sys.stdout.flush()
                            if 'yields' in getYields_dic.keys():
                                yields=getYields_dic['yields']
                            else:
                                yields=[0.0, 0.0]
                        pf.close()
                        print('\t reading the pickle file, get_yields returned [total SpY, total RY] = ', yields)
                    else:
                        print('\t WARNING! could not read yield from get_yields pickle file. Set to zero')
                        yields=[0.0, 0.0]
                    
                    #overwrite spY value if mode is 'calculate'
                    if self.spYieldMode[i]=='calculate':
                        self.spYield[i]=float(yields[0])
                    if self.rYieldMode[i]=='calculate':
                        self.rYield[i]=float(yields[1])

                    #4) save tridyn.dat
                    #append output to allTridyn.dat for each species
                    ft_output_profile_final=self.FT_OUTPUT_PROFILE_FINAL+'_'+prj
                    tempfile = open(self.FT_OUTPUT_PROFILE_TEMP,"r")
                    f = open(ft_output_profile_final, "a")
                    f.write('%s%s \n' %(tempfile.read().rstrip('\n'),self.maxRangeXolotl[i]))                    
                    f.close()
                    tempfile.close()

                    ## keep copies of tridyn.dat in timeFolder
                    ft_output_profile_temp_prj=self.FT_OUTPUT_PROFILE_TEMP+'_'+prj #for each species
                    shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, timeFolder+'/'+ft_output_profile_temp_prj)
                    #shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP,ft_output_profile_temp_prj)
                
                    #5) MOVE FOLDERS TO DIRECTORY WITH TIME-STAMP & RENAME FOR (Tg,Prj) SPECIES  
                
                    shutil.move(self.ftridyn['outputPath']+'/'+self.FT_OUTPUT_FOLDER,timeFolder+'/'+self.FT_OUTPUT_FOLDER+'_'+prj+'W')                
                    self.services.update_state()
                    print(' ')
                    print('... done with F-TRIDYN for {}'.format(prj))
                    print(' ')
                    sys.stdout.flush()

                #if flux fraction == 0 or all weight angles == 0.0:
                else:
                    print(('Skip running FTridyn for {0}'.format(prj))) #, as fraction in plasma is {1}\n'.format(prj, self.gitr['fluxFraction'][i]))
                    if self.fluxFraction[i]==0:
                        print('\t as fraction in plasma is {}\n'.format(self.fluxFraction[i]))
                    if all(k==0 for k in self.weightAngle[i]):
                        print('\t as all angle weights are zero\n')
                    self.spYield[i]=0.0
                    self.rYield[i]=1.0
                    self.maxRangeXolotl[i]=0.0
                    outputFTFile=open(self.FT_OUTPUT_PROFILE_TEMP, "w")
                    outputFTFile.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ")
                    outputFTFile.close()

                    #keep copies of tridyn.dat in timeFolder    
                    ft_output_profile_temp_prj=self.FT_OUTPUT_PROFILE_TEMP+'_'+prj #for each species
                    shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, timeFolder+'/'+ft_output_profile_temp_prj)
                    #shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP,ft_output_profile_temp_prj)                    
                    sys.stdout.flush()
                    
            #end of for loop:
            sys.stdout.flush()

            print('... done with F-TRIDYN for all species.')
            print(' ')
            print('species-independent analysis and consolidate output :')
            print(' ')
            
            ######species independent ############

            #6) write sputtering yields to file so they can be used by Xolotl

            yieldString=str(time)
            print('Sputtering and Reflection Yields due to:')
            
            for i in range(len(self.plasmaSpecies)): #self.plasmaSpecies.iteritems():  
                prj=self.plasmaSpecies[i]
                print(('\t{0} :  spY = {1} and rY = {2} '.format(prj,self.spYield[i], self.rYield[i])))
                prjYieldsString=prj+' ' +str(self.spYield[i])+' '+str(self.rYield[i])+' '
                yieldString+=prjYieldsString
                
            #write sp Yields to file (temp) and append to spYield output (final)
            spTempfile = open(self.FTX_SPUT_YIELDS_FILE_TEMP,"w+")
            spTempfile.write(yieldString)
            spFile = open(self.FTX_SPUT_YIELDS_FILE_FINAL, "a+")
            spFile.write(yieldString)
            spFile.write('\n')
            spFile.close()
            spTempfile.close()
            
            shutil.copyfile(self.FTX_SPUT_YIELDS_FILE_TEMP,timeFolder+'/'+self.FTX_SPUT_YIELDS_FILE_TEMP)
            shutil.copyfile(self.FTX_SPUT_YIELDS_FILE_FINAL,timeFolder+'/'+self.FTX_SPUT_YIELDS_FILE_FINAL) #perhaps unnecessary
            
            #if len(maxDepth)==0 and max(self.spYield)==0:

            print(' ')
            print('\t self.maxRangeXolotl = ', self.maxRangeXolotl ,' and max of self.maxRangeXolotl =', max(self.maxRangeXolotl))
            print('\t spYield = ', self.spYield, ' and max(self.spYield) = ', max(self.spYield))
            
            if max(self.maxRangeXolotl)==0 and max(self.spYield)==0:
                print("nothing was implanted or sputtered")
                print("likely all weights are zero")
                print("END OF THIS SIMULATION")
                sys.stdout.flush()
                return
            sys.stdout.flush()
            
            #7) write format tridyn.dat to include W redep in Xolotl
            
            print(' ')
            if os.path.exists(self.FT_OUTPUT_PROFILE_TEMP):
                os.remove(self.FT_OUTPUT_PROFILE_TEMP)

            # allow user to specify the tridynDat_model used to write the tridyn.dat file independently of the xolotl_v option
            # if the option is not specified, fall back to the xolotl_v option
            if 'tridynDat_model' in self.driver.keys():
                tridynDat_model = self.driver['tridynDat_model']
                if not (tridynDat_model == 1 or tridynDat_model == 2):
                    print('\t WARNING: tridynDat_model version not recognized: ')
                    print(f'\t \t tridynDat_model = {tridynDat_model}')
                    tridynDat_model = self.driver['xolotl_v']
                    print(f'\t Trying to continue with tridynDat_model = {tridynDat_model} (same as xolotl_v)')
            else:
                tridynDat_model = self.driver['xolotl_v']
                
            write_tridynDat.write_tridynDat(outFile=self.FT_OUTPUT_PROFILE_TEMP, tridynDat_model=tridynDat_model,
                                            plasmaSpecies=self.plasmaSpecies, timeFolder=timeFolder, maxRangeXolotl=self.maxRangeXolotl,
                                            fluxFraction=self.fluxFraction, rYield=self.rYield, xp_parameters=self.xp.parameters, print_test=self.print_test)
            sys.stdout.flush()

            if os.path.exists(self.FT_OUTPUT_PROFILE_TEMP):
                print('succesfully created ', self.FT_OUTPUT_PROFILE_TEMP)
            
            #save tridyn.dat (combined for all species) with time-stamp:
            print('\t to keep ', self.FT_OUTPUT_PROFILE_TEMP, ' for every loop, save as ', self.FT_OUTPUT_PROFILE_TEMP+'_t'+str(time))
            shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, self.FT_OUTPUT_PROFILE_TEMP+'_t'+str(time))
            print(' ')
            
            #compress output
            if self.driver['ZIP_FTRIDYN_OUTPUT']=='True':
                print(('zip F-TRIDYNs output: {} \n'.format(timeFolder)))
                zippedTimeFolder=timeFolder+'.zip'
                zip_output='zipOutputTimeFolder.txt'
                zipString='zip -r %s %s >> %s ' %(zippedTimeFolder, timeFolder, zip_output)
                subprocess.call([zipString], shell=True)
                shutil.rmtree(timeFolder)

            else:
                print(('leaving {} uncompressed \n'.format(timeFolder)))

            print(' ')
            sys.stdout.flush()
            self.services.update_state()

            
            #make sure at least one of the species was implanted:
            if all(k==0 for k in self.maxRangeXolotl):
                print("WARNING: none of the species was implanted")
                print("EXIT SIMULATION")
                sys.stdout.flush()
                return

            #if something was implanted:
            #Xolotl parameter modifications that need to be done at every loop:
            
            ######################################
            ############## RUN XOLOTL ############
            ######################################
            print('-----------')
            print('  XOLOTL:  ')
            print('-----------')
            print(' ')
            sys.stdout.flush()

            #list modules right before running Xolotl:
            #if self.print_test:
            #    modList=os.system('module list')
            #    print('\t TEST: we are running Xolotl with the following modules loaded: ', modList)
            #    print('\t load cray-hdf5-parallel')
            #os.system('module load cray-hdf5-parallel')
            #if self.print_test:
            #    modList=os.system('module list')
            #    print('\t TEST: we are running Xolotl with the following modules loaded: ', modList)
            #    print(' ')
            #sys.stdout.flush()
            
            #calculate effective sputtering yield; i.e., weighted by relative flux of W-to-He
            totalSpYield=0
            for i in range(len(self.plasmaSpecies)): #self.plasmaSpecies.iteritems():
                prj=self.plasmaSpecies[i]
                totalSpYield+=(float(self.fluxFraction[i])*float(self.spYield[i])) #self.fluxFraction[i]
            print(('total weighted sputtering yield = {} (passed to Xolotl)'.format(totalSpYield)))
            
            for i in range(len(self.plasmaSpecies)):
                prj=self.plasmaSpecies[i]
                print(('\t contribution of {0} to total sputtering yield = {1} '.format( prj, float(self.fluxFraction[i])*float(self.spYield[i]))))
            self.xp.parameters['sputtering'] = totalSpYield            
            print(' ')
            sys.stdout.flush()
            
            #time and time-step related parameters - moved to beginning of while loop
            #final_time = time+self.time['LOOP_TIME_STEP']
            #rounded_final_time = round(final_time, self.time_decimal)
            ##if (self.print_uq):
            #print(f"Rounding final time from {final_time} to {rounded_final_time} before running Xolotl")
            #self.xp.parameters['petscArgs']['-ts_final_time']=rounded_final_time

            print('\n')
            print('run XOLOTL: ')
            print(('\t from t = {}'.format(time)))
            print(('\t to t = {}'.format(self.xp.parameters['petscArgs']['-ts_final_time'])))
            print(('\t and time-step = {} '.format( self.time['LOOP_TIME_STEP'])))

            #note on network file:
            # XOLOTL_NETWORK_FILE is the file ; xp.parameters['networkFile'] is the dictionary value (for params.txt)
            # better not mix them, because one can exist without the other (e.g., if networkFile line not needed in params.txt) 
            
            if self.driverMode == 'INIT':
                if 'networkFile' in self.xp.parameters:
                    print('\t init mode: networkFile defined in Xolotl input parameters of config file')
                    if os.path.exists(self.INPUT_DIR+'/'+self.NETWORK_FILE):
                        if os.path.getsize(self.INPUT_DIR+'/'+self.NETWORK_FILE) != 0:
                            print('\t \t network file exists and is not empty')
                            print('\t \t load networkFile file \n')                
                        else:
                            print('\t \t network file exists but its empty')
                            print('\t \t skip loading networkFile file')
                            print('\t \t remove networkFile from Xolotl parameters')
                            del self.xp.parameters['networkFile']
                    else:
                        print('\t \t but network file does not exist')
                        print('\t \t skip loading networkFile file')
                        print('\t \t remove networkFile from Xolotl parameters')
                        del self.xp.parameters['networkFile']
                else:
                    #no network file in input params -- do not try to load
                    print('\t init mode: no network defined in Xolotl input parameters of config file')
                    print('\t \t WARNING: will create (not load) the network \n')

            elif self.driverMode == 'RESTART' or self.driverMode == 'NEUTRAL':
                #add (or replace) networkFile line to parameter file
                print('\t restart (or neutral) mode: modify xolotl parameters that might change at every loop, including adding the networkFile \n')
                self.xp.parameters['networkFile'] = self.XOLOTL_NETWORK_FILE

                ## check if we need to keep netParam in parameter file of restart
                if self.driver['xolotl_v'] == 2:
                    print('\t xolotl v',self.driver['xolotl_v'] , ' requires netParam and grouping even in restart')
                    print('\t will not delete them') # (regardless of value of netParam_restart and grouping_restart)')
                elif self.driver['xolotl_v'] == 1:
                    print('\t xolotl v',self.driver['xolotl_v'], ' does not require netParam / grouping in restart (info contained in netFile)')
                    if 'netParam' in self.xp.parameters:
                        print('\t \t remove netParam') # (regardless of value of netParam_restart)')
                        del self.xp.parameters['netParam'] 
                    else:
                        print('\t \t netParam does not exist in the xolotl dictionary. No need to delete it')
                    if 'grouping' in self.xp.parameters:
                        print('\t \t remove grouping') # (regardless of value of grouping_restart)')
                        del self.xp.parameters['grouping']
                    else:
                        print('\t \t grouping does not exist in the xolotl dictionary. No need to delete it')
                print(' ')
            sys.stdout.flush()
                    
            #determine if he_conc true/false ; if true, add '-he_conc' to petsc arguments 
            if self.driver['XOLOTL_HE_CONC']=='Last':
                if time+1.5*self.time['LOOP_TIME_STEP']>end_time:  #*1.5, to give marging of error
                    self.petsc_heConc=True
                    print('printing He concentrations in the last loop')
                elif time<(end_time-self.time['LOOP_TIME_STEP']):
                    self.petsc_heConc=False
            elif self.driver['XOLOTL_HE_CONC']=='True':
                print('\t he_conc printed in this (and every) loop')
                self.petsc_heConc=True
            elif self.driver['XOLOTL_HE_CONC']=='False':
                print('\t he_conc not printed in this (or any other) loop')
                self.petsc_heConc=False

            if self.petsc_heConc:
                self.xp.parameters['petscArgs']['-helium_conc'] = ''

            sys.stdout.flush()

            #-check_collapse option in petcs args:
            #exit status is printed to solverStatus.txt  'good' (successful run); 'collapsed' (ts below threshold) ; 'diverged' otherwise        
            #Xolotl is launched again (same paramter and network files) until run successfully, or up to maxCollapseLoop tries,

            # flag to keep track of the success of keepLastTS (because sometimes the network file is corrupted)
            keep_last_ts_success = False
            keep_last_ts_loop_number = 0

            # make sure we save a copy of the network file in case we need to restart because keepLastTS failed
            temp_network_file = os.path.join(os.path.dirname(self.XOLOTL_NETWORK_FILE), "tempNetworkFile.h5")
            shutil.copyfile(self.XOLOTL_NETWORK_FILE, temp_network_file)
            print(f"\u2B95 Successfully copied network file {self.XOLOTL_NETWORK_FILE} to temporary network file {temp_network_file}")

            network_size_file=cwd+'/worker_netFile_size.txt' #this one needs the full path
            
            # loop until keepLastTS is successful
            while not keep_last_ts_success:
                keep_last_ts_loop_number += 1

                if keep_last_ts_loop_number > 1:
                    print(f"\u2B95 keep_last_ts_success is False, starting new loop (loop number is {keep_last_ts_loop_number})")

                # start from a fresh network file
                shutil.copyfile(temp_network_file, self.XOLOTL_NETWORK_FILE)
                print(f"\u2B95 Successfully copied network file {temp_network_file} to temporary network file {self.XOLOTL_NETWORK_FILE}")

                n_overgrid_loops=0
                
                while self.xolotlExitStatus=='collapsed' or self.xolotlExitStatus=='overgrid':

                    if self.xolotlExitStatus=='collapsed':
                        self.collapsedLoops+=1
                        print(' ')

                    #set a maximum number of tries
                    if self.collapsedLoops<=int(self.driver['MAX_COLLAPSE_LOOPS']):
                        #this first call can probably be removed, and just keep while over netFile size
                        self.services.call(xolotl, 'init', timeStamp, dTime=time, time_decimal=self.time_decimal, xParameters=self.xp.parameters, network_size_file=network_size_file, output_file=outFile, print_test=self.print_test)
                        if os.path.exists(network_size_file):
                            print('found network_size_file')
                            worker_network_sizeFile=open(network_size_file, 'r')
                            worker_netFile_size=int(worker_network_sizeFile.readline())
                            worker_network_sizeFile.close
                            print('\t workers network file size is ', worker_netFile_size)
                            driver_netFile_size=int(os.path.getsize(self.XOLOTL_NETWORK_FILE))
                            print('\t driver network file size is ', driver_netFile_size)
                            temp_netFile_size=int(os.path.getsize(temp_network_file))
                            print('\t temp network file size is ', temp_netFile_size)
                            inputDir_netFile_size=int(os.path.getsize(self.INPUT_DIR+'/'+self.XOLOTL_NETWORK_FILE))
                            print('\t input dir network file size is ', inputDir_netFile_size)
                            sys.stdout.flush()
                            #while worker_netFile_size != driver_netFile_size:
                            while worker_netFile_size != inputDir_netFile_size:
                                print ('input dir and worker network sizes differ')
                                print(' copy input dir network file here, update state and try running worker init again')
                                shutil.copyfile(self.INPUT_DIR+'/'+self.XOLOTL_NETWORK_FILE, self.XOLOTL_NETWORK_FILE)
                                self.services.update_state()
                                self.services.call(xolotl, 'init', timeStamp, dTime=time, time_decimal=self.time_decimal, xParameters=self.xp.parameters, network_size_file=network_size_file, output_file=outFile, print_test=self.print_test)
                                if os.path.exists(network_size_file):
                                    worker_network_sizeFile=open(network_size_file, 'r')
                                    worker_netFile_size=int(worker_network_sizeFile.readline())
                                    worker_network_sizeFile.close
                                    print('\t updated workers network file size is ', worker_netFile_size)
                                    driver_netFile_size=int(os.path.getsize(self.XOLOTL_NETWORK_FILE))
                                    print('\t updated driver network file size is ', driver_netFile_size)
                                    temp_netFile_size=int(os.path.getsize(temp_network_file))
                                    print('\t updated temp network file size is ', temp_netFile_size)
                                    inputDir_netFile_size=int(os.path.getsize(self.INPUT_DIR+'/'+self.XOLOTL_NETWORK_FILE))
                                    print('\t updated input dir network file size is ', inputDir_netFile_size)
                                    sys.stdout.flush()
                                else:
                                    print('WARNING: could not find network size file')
                                    print('\t continue at your own risk ')
                            print('network file sizes are the same. continue with worker:step')
                        else:
                            print('WARNING: could not find network size file')
                            print('\t continue at your own risk ')
                        sys.stdout.flush()
                            
                        self.services.call(xolotl, 'step', timeStamp, dTime=time, time_decimal=self.time_decimal, xHe_conc=self.petsc_heConc, xParameters=self.xp.parameters, output_file=outFile, dZipOutput=self.driver['ZIP_XOLOTL_OUTPUT'], n_overgrid_loops=n_overgrid_loops, xolotl_num_tries=self.driver['XOLOTL_NUM_TRIES'], print_test=self.print_test)

                        sys.stdout.flush()
                        self.services.stage_state()

                        statusFile=open(self.XOLOTL_EXIT_STATUS, "r")
                        self.xolotlExitStatus=statusFile.read().rstrip('\n')

                        print(('\t Xolotl ended simulation with status : {}\n'.format(self.xolotlExitStatus)))

                        if self.xolotlExitStatus=='good':
                            print(('Xolotl successfully executed after {} tries'.format(self.collapsedLoops)))
                            print('\t continue with IPS simulation \n')

                        elif self.xolotlExitStatus=='diverged':
                            print('ERROR: XOLOTL SOLVER DIVERGED ')
                            print('\t END IPS SIMULATION \n')
                            quit()
                        elif self.xolotlExitStatus=='overgrid':
                            n_overgrid_loops+=1
                            print('WARNING: XOLOTL OVERGRID ')
                            print('\t RUNNING transgerGrid...')
                            #xolotlStop already copied as _overgrid within xolotl_comp; no need to copy here again
                            sys.stdout.flush()
                            
                            #also need to update the values of grid and voidPortion to match what's in the network file:
                            #Most likely this will only work if using 'grid', not 'gridParam' in network file!!
                            shutil.copyfile(self.XOLOTL_NETWORK_FILE,self.XOLOTL_NETWORK_FILE+'_overgrid_'+str(n_overgrid_loops))
                            shutil.move('xolotlStop.h5', self.XOLOTL_NETWORK_FILE)
                            [newGridSize, newVoidPortion] = transferGrid.transferGrid(self.XOLOTL_NETWORK_FILE,print_test=self.print_test)
                            sys.stdout.flush()
                            
                            if os.path.exists(self.XOLOTL_NETWORK_FILE):
                                #make sure xolotlStop exists
                                shutil.copyfile(self.XOLOTL_NETWORK_FILE,'xolotlStop.h5')
                            else:
                                print('\t WARNING: networkFile does not exists after transferGrid')
                                print('\t keep old one')
                                shutil.copyfile(self.XOLOTL_NETWORK_FILE+'_overgrid_'+str(n_overgrid_loops),self.XOLOTL_NETWORK_FILE)
                            print('\t in file ', self.XOLOTL_NETWORK_FILE)    
                            print('\t \t updated the values of grid to ', newGridSize)
                            print('\t \t updated the values of voidPortion to ', newVoidPortion)
                            print('\t DONE RUNNING transferGrid!')
                            sys.stdout.flush()
                            if self.driver['xolotl_v']==1:
                                self.xp.parameters['grid'][0] = newGridSize
                            elif self.driver['xolotl_v']==2:
                                self.xp.parameters['gridParam'][0] = newGridSize
                            self.xp.parameters['voidPortion'] = newVoidPortion
                            self.services.update_state()
                            
                        elif self.xolotlExitStatus=='collapsed':
                            print('\t WARNING: simulation exited loop with status collapse')
                            print(('\t try number {0} out of {1}\n'.format(self.collapsedLoops,self.maxCollapseLoops)))                    

                        else:
                            print('\t WARNING: Xolotl exit status UNKOWN -- IPS simulation continues \n')
                        sys.stdout.flush()
                        
                    else: #reached maximum number of tries for collapsing time steps
                        print('\t ERROR: reached maximum number of tries for collapsing time steps without a successful run')
                        print('END IPS SIMULATION \n')
                        quit()
                sys.stdout.flush()
                
                #UPDATE: Xolotl generates HDF5 file of TRIDYN, copied as 'last_TRIDYN_toBin.h5'
                #thus no need to copy it; and binTRIDYN will transform it to text file, 'last_TRIDYN.dat' 
                print('bin Xolotls output:')
                sys.stdout.flush()
                if (self.driver['xolotl_v']==1):
                    binTRIDYN.v1(inFile='last_TRIDYN_toBin.h5', outFile='last_TRIDYN.dat', print_test=self.print_test) #instead of binTRIDYN.binTridyn()                 
                    print('...succesfully ran binTRIDYN for xolotl v1')
                elif(self.driver['xolotl_v']==2):
                    binTRIDYN.v2(inFile='last_TRIDYN_toBin.h5', outFile='last_TRIDYN.dat', print_test=self.print_test) #formerly binTRIDYN_tempGrid
                    print('...succesfully ran binTRIDYN for xolotl v2')
                print(' ')
                
                #store xolotls profile output for each loop (not plasma state)          
                #UPDATE: 'toBin' is .h5 format instead of .dat
                currentXolotlOutputFileToBin='last_TRIDYN_toBin_%f.h5' %time
                shutil.copyfile('last_TRIDYN_toBin.h5', currentXolotlOutputFileToBin)
                currentXolotlOutputFile='last_TRIDYN_%f.dat' %time
                shutil.copyfile('last_TRIDYN.dat', currentXolotlOutputFile)
                
                #append output:
                #retention
                if (self.print_uq):
                    print(f">>> PR: the current directory is {os.getcwd()}")
                exists_str = "exists" if os.path.isfile(self.XOLOTL_RETENTION_TEMP) else "does not exist !!!"
                if (self.print_uq):
                    print(f">>> PR: XOLOTL_RETENTION_TEMP = '{self.XOLOTL_RETENTION_TEMP}', this file {exists_str}")
                tempfileRet = open(self.XOLOTL_RETENTION_TEMP,"r")
                fRet = open(self.XOLOTL_RETENTION_FINAL, "a")
                fRet.write(tempfileRet.read())
                fRet.close()
                tempfileRet.close()
                exists_str = "exists" if os.path.isfile(self.XOLOTL_RETENTION_FINAL) else "does not exist !!!"
                if (self.print_uq):
                    print(f">>> PR: XOLOTL_RETENTION_FINAL = '{self.XOLOTL_RETENTION_FINAL}', this file {exists_str}")
                
                #surface
                exists_str = "exists" if os.path.isfile(self.XOLOTL_SURFACE_TEMP) else "does not exist !!!"
                if (self.print_uq):
                    print(f">>> PR: XOLOTL_SURFACE_TEMP = '{self.XOLOTL_SURFACE_TEMP}', this file {exists_str}")
                tempfileSurf = open(self.XOLOTL_SURFACE_TEMP,"r")
                fSurf = open(self.XOLOTL_SURFACE_FINAL, "a")
                fSurf.write(tempfileSurf.read())
                fSurf.close()
                tempfileSurf.close()
                exists_str = "exists" if os.path.isfile(self.XOLOTL_SURFACE_FINAL) else "does not exist !!!"
                if (self.print_uq):
                    print(f">>> PR: XOLOTL_SURFACE_FINAL = '{self.XOLOTL_SURFACE_FINAL}', this file {exists_str}")

                    print(f">>> PR: applying fix here:")
                for (tmp_file, final_file) in [(self.XOLOTL_RETENTION_TEMP, self.XOLOTL_RETENTION_FINAL), (self.XOLOTL_SURFACE_TEMP, self.XOLOTL_SURFACE_FINAL)]:
                    if (self.print_uq):
                        print(f">>> PR: got tmp_file '{tmp_file}' and final_file '{final_file}' ")
                    with open(final_file, "a") as out_file:
                        if (self.print_uq):
                            print(f">>> PR: successfully opened final_file '{final_file}'")
                        with open(tmp_file, "r") as in_file:
                            if (self.print_uq):
                                print(f">>> PR: successfully opened tmp_file '{tmp_file}'")
                            out_file.write("\n")
                            out_file.write(in_file.read())
                            if (self.print_uq):
                                print(f">>> PR: successfully written final_file '{final_file}'")
                if (self.print_uq):
                    if os.path.isfile(self.XOLOTL_RETENTION_FINAL) and os.path.isfile(self.XOLOTL_SURFACE_FINAL):
                        print(f">>> PR: both XOLOTL_RETENTION_FINAL '{self.XOLOTL_RETENTION_FINAL}' and XOLOTL_SURFACE_FINAL '{self.XOLOTL_SURFACE_FINAL}' exist")
                    else:
                        print(f">>> PR: either XOLOTL_RETENTION_FINAL '{self.XOLOTL_RETENTION_FINAL}' or XOLOTL_SURFACE_FINAL '{self.XOLOTL_SURFACE_FINAL}' does not exist !!!")
                    print(f">>> PR: done applying fix")

                    print(f">>> PR: cwd is '{cwd}'")

                #save network file with a different name to use in the next time step
                currentXolotlNetworkFile='xolotlStop_%f.h5' %time
                shutil.copyfile('xolotlStop.h5',currentXolotlNetworkFile)
                ## try using keepLastTS to produce netfile with only info from the last TS
                print('produce new network file using xolotlStop:')
                try:
                    iF= 'xolotlStop.h5' #cwd+'/xolotlStop.h5' ; rm cwd from paths
                    oF= self.XOLOTL_NETWORK_FILE #cwd+'/'+self.XOLOTL_NETWORK_FILE ; rm cwd from paths
                    os.remove(oF) #can not exist & it's copied as w/ time-stamp above
                    if self.print_test:
                        print('\t run keepLastTS with: ')
                        print('\t \t inFile = ', iF)
                        print( '\t \t outFile = ', oF)
                    sys.stdout.flush()
                    keepLastTS.keepLastTS(inFile=iF, outFile=oF, print_test=self.print_test)
                    if self.print_test:
                        print('\t ... keepLastTS returned succesfully')
                    print('done writing a new network file')
                    print(' ')
                    sys.stdout.flush()
                    keep_last_ts_success = True
                #if fails, use old method of copying entire xolotlStop as networkFile
                except Exception as e:       
                    print(e)
                    print('\t running keepLastTS failed')
                    # print('\t just copy xolotlStop as networkFile')
                    # shutil.copyfile('xolotlStop.h5',self.XOLOTL_NETWORK_FILE)
                    print("\u2B95 keepLastTS failed, trying this loop again")
                # print('done writing a new network file')
                # print(' ')
                # sys.stdout.flush()
            
            # copy last_TRIDYN.dat and netowrk file to input dir as well:
            print('copy last_TRIDYN.dat and ', self.XOLOTL_NETWORK_FILE, 'to: ')
            print('\t',  self.INPUT_DIR )
            shutil.copyfile('last_TRIDYN.dat', self.INPUT_DIR+'/last_TRIDYN.dat')
            source_ = self.XOLOTL_NETWORK_FILE
            dest_ = self.INPUT_DIR+'/'+self.XOLOTL_NETWORK_FILE
            with open(source_, "r") as f1:
                with open(dest_, "w") as f2:
                    shutil.copyfile(source_, dest_)
            print(' ')


            #done running this loop and saving output
            if (self.print_test):
                print('done running this loop and saving output \n')
                
            ## update driver mode after the 1st loop, from INIT to RESTART
            if self.driverMode != 'NEUTRAL':
                self.driverMode = 'NEUTRAL'
                print(('switched driverMode to {} \n'.format(self.driverMode)))

            sys.stdout.flush()
            #using while instead of loop to accomodate variable drive time step 
            #--> update time explicitely before (possibly) increasing time step
            time+=self.time['LOOP_TIME_STEP']
            print(' ')
            print(('after loop {}, check for updates in time steps '.format(self.time['LOOP_N'])))

            #update Xolotl and driver time steps if needed
            if self.time['LOOP_TS_FACTOR'] != 1:
                if (self.time['LOOP_N']%self.time['LOOP_TS_NLOOPS']==0):
                    print('\t update driver time step and start_stop ')
                    self.time['LOOP_TIME_STEP']*=self.time['LOOP_TS_FACTOR']
                    self.xp.parameters['petscArgs']['-start_stop']*=self.time['LOOP_TS_FACTOR']
                    print(('\t multiplied time step and start_stop by {} '.format(self.time['LOOP_TS_FACTOR']))) 
                    print(('\t for a new time step = {0} and start_stop ={1} \n'.format(self.time['LOOP_TIME_STEP'], self.xp.parameters['petscArgs']['-start_stop'])))

                else:
                    print(('\t no update to driver time step ({0}) or start_stop ({1}) \n'.format(self.time['LOOP_TIME_STEP'] , self.xp.parameters['petscArgs']['-start_stop'])))
            else:
                print(('\t driver time step (({0}) and start_stop ({1}) unchanged (factor=1) \n'.format( self.time['LOOP_TIME_STEP'], self.xp.parameters['petscArgs']['-start_stop'])))

            #no need for this here anymore ; check done in the very beginning of the while loop
            #if time+self.time['LOOP_TIME_STEP']>end_time: 
            #    self.time['LOOP_TIME_STEP']=end_time-time
            #    self.xp.parameters['petscArgs']['-start_stop']=end_time/10.0
            #    if (self.print_uq):
            #        print(f"PR: >>> Updated -start_stop to {self.xp.parameters['petscArgs']['-start_stop']:.4f}")
            #    print(' ')
            #    print('\t time step longer than needed for last loop ')
            #    print(('\t adapting driver time step to {} to reach exactly endTime '.format(self.time['LOOP_TIME_STEP'])))
            #    self.xp.parameters['petscArgs']['-start_stop']=self.time['LOOP_TIME_STEP']/10.0
            #    print(('\t and Xolotls data is saved every (start_stop) = {} \n'.format( self.xp.parameters['petscArgs']['-start_stop'])))
            #else:
            #    self.xp.parameters['petscArgs']['-start_stop']=(time+self.time['LOOP_TIME_STEP'])/10.0
            #    if (self.print_uq):
            #        print(f"PR: >>> Updated -start_stop to {self.xp.parameters['petscArgs']['-start_stop']:.4f}")

            if self.driver['XOLOTL_MAXTS_FACTOR'] != 1:
                if (self.time['LOOP_N']%self.driver['XOLOTL_MAXTS_NLOOPS']==0):
                    print(('\t change in Xolotls maximum time step after loop {} '.format( self.time['LOOP_N'])))
                    self.xp.parameters['petscArgs']['-ts_adapt_dt_max']*=self.driver['XOLOTL_MAXTS_FACTOR']
                    print(('\t multiply time step by {0}, for a new time step = {1}'.format(self.driver['XOLOTL_MAXTS_FACTOR'] , self.xp.parameters['petscArgs']['-ts_adapt_dt_max'])))
                    if self.xp.parameters['petscArgs']['-ts_adapt_dt_max'] > self.driver['XOLOTL_MAX_TS']:
                        self.xp.parameters['petscArgs']['-ts_adapt_dt_max']=self.driver['XOLOTL_MAX_TS']
                        print(('\t Xolotls time-step reached the maximum allowed; set to limit value, {}'.format(self.driver['XOLOTL_MAX_TS'])))
                else:
                    print(( '\t continue with xolotl dt max {}'.format(self.xp.parameters['petscArgs']['-ts_adapt_dt_max'])))
                print(' ')
            else:
                print(('\t Xolotls max time step ({}) unchanged (factor=1)'.format(self.xp.parameters['petscArgs']['-ts_adapt_dt_max'])))

            print(' ')
            
            sys.stdout.flush()
            self.services.update_state()
            
        if (self.print_test):
            print('\n')
            print('Done with looping over time')
            print("Look driver's output files, as defined in config file:")
            try:
                self.OUTPUT_FILES
                print('\t Output files are defined as:')
                print('\t', self.OUTPUT_FILES)
            except:
                print('\t No output files defined in config file')
            print(' ')
            sys.stdout.flush()
            
    def finalize(self, timeStamp=0.0,**keywords):

        print('  FT-X driver:finalize called')
        print('\t with keywords: ',keywords)
        
        print('\t output file of the FT-X workflow:')
        if 'LOG_FILE' in keywords:
            logFile=keywords['LOG_FILE']
            outFile= cwd+'/'+logFile
            print('\t \t log file defined in keywords: ')
            print('\t \t ', outFile)
            outF=open(outFile , 'a')
            sys.stdout = outF
        else:
            try:
                self.LOG_FILE
                logFile = self.LOG_FILE
                outFile=logFile
                print('\t \t log file defined in config file', outFile)
                outF=open(outFile , 'a')
                sys.stdout = outF
            except:
                print('\t \t No log file defined; using default sys.stdout')
                outFile=None
        sys.stdout.flush()
        
        #can we add compressing output here? e.g., last_TRIDYN, xolotlStop...
        #and remove large output files? e.g., FTRIDYN.zip

        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')
        self.services.call(ftridyn, 'finalize', timeStamp)
        self.services.call(xolotl, 'finalize', timeStamp)
