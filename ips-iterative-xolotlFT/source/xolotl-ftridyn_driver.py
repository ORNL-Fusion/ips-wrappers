#! /usr/bin/env python

from  component import Component
import sys
import os
import os.path
import subprocess
import numpy
import shutil
import translate_xolotl_to_ftridyn
#import translate_ftridyn_to_xolotl_launch
#import get_yields_launch
import binTRIDYN
import param_handler
import traceback
import transferGrid
import pickle

class xolotlFtridynDriver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeStamp=0.0, **keywords):

        cwd = self.services.get_working_dir()

        print(' ')
        print(('FT-X driver:init called with: ',keywords))
        print('\t output file of the FT-X workflow:')
        if 'LOG_FILE' in keywords:
            logFile=keywords['LOG_FILE']
            outFile=cwd+'/'+logFile
            print('\t \t log file defined in keywords: ', outFile)
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
                outfile=None

        print('xolotl-ftridyn_driver: init \n')

        #stage input files
        print(('staging input files {} \n'.format(self.INPUT_FILES)))
        self.services.stage_input_files(self.INPUT_FILES)
        print('\t ...input files staged succesfully')
        print(('input directory for this simulation is {} \n'.format( self.INPUT_DIR)))

        plasma_state_file = self.services.get_config_param('PLASMA_STATE_FILES')
        plasma_state_list = plasma_state_file.split()
        for index in range(len(plasma_state_list)):
            open(plasma_state_list[index], 'a').close()                
        #A MORE ELEGANT WAY --  FOR THE FUTURE
            #for file in plasma_state_list:
            #    open(file, 'a').close()
            
        self.services.update_plasma_state()
        self.services.stage_plasma_state()

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

        except Exception as e:
            print(e)
            print('no [TIME_PARAMETERS] defined in config file')
            sys.stdout.flush()
            
            print('INIT, STEP, END TIME and other paramters (LPPS_TS, START_MODE) are defined in driver section of config file')
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
                
        #Check if/what's given in input file
        #overwrite config file values
        try:
            self.TIME_FILE
            print('\t use input from time file defined in config file (overwirte config file values)')
            print(self.INPUT_DIR+'/'+self.TIME_FILE)
            self.time=param_handler.read(self.INPUT_DIR+'/'+self.TIME_FILE)


            ## for testing purposes, a more thorough print line:
            ## use these three afterwards
            #for k,v, in self.time.items():
            #    print(('\t {0} : {1}'.format(k, v)))
            #print(' ')

            for k,v in self.time.items():
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
            
        except Exception as e:
            print(e)
            print('no TIME_FILE defined. Use time-parameters defined in config file')
            sys.stdout.flush()

        
        print('time characteristics are')
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

        
        #### XOLOTL PARAMETERS ##### 

        print('\n')
        print('XOLOTL paramters: \n')

        #get dimension to read the correct input parameter template
        dim=int(self.XOLOTL_INPUT_PARAMETERS['dimensions'])
        #xolotl_param_template=self.XOLOTL_PARAM_TEMPLATE_PATH+'/paramXolotl_'+str(dim)+'D.txt'
        xolotl_param_template=self.INPUT_DIR+'/paramXolotl_'+str(dim)+'D.txt' 
        print(('\t reading Xolotl default parameters from {} \n'.format(xolotl_param_template)))

        self.xp = param_handler.xolotl_params()
        self.xp.read(xolotl_param_template)
        print(('\t running Xolotl in {} D \n'.format(dim)))

        #overwrite default Xolotl parameters that are specified in ips.config
        print('modify XOLOTL paramters with parameters read from the simulations config file: \n')

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

        #if not coupling, delete -tridyn from petsc arguments to not print TRIDYN_*.dat files
        #if self.driver['FTX_COUPLING']=='False':
        #    del self.xp.parameters['petscArgs']['-tridyn']

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

        print(('processes included in Xolotl are: {}\n'.format(processString.strip())))
        self.xp.parameters['process']=processString.strip()

        sys.stdout.flush()
        
        ### PLASMA RELATED PARAMETERS ###
        ## code agnostic (gitr, solps... all handled together)
        
        #keep track if plasma input exists (file or config file) with string: F = false, T = true
        plasmaExists='F'
        self.plasma = {}
        try:
            self.PLASMA_INPUT_PARAMETERS
            plasmaExists='T'
            print('Read PLASMA inputs defined in config file.')

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

        except Exception as e:
            print(e)
            print('no PLASMA_INPUT_PARAMETERS defined in config file. ')
            sys.stdout.flush()
                    
        try:
            self.PLASMA_OUTPUT_FILE
            plasmaExists='T'
            print('use input from PLASMA, as its output file is defined in config file')
            print(('read PLASMA parameters (from file) : {}\n'.format(self.INPUT_DIR+'/'+self.PLASMA_OUTPUT_FILE)))
            print('these inputs will overwrite values given in config file ')
            self.plasma=param_handler.read(self.INPUT_DIR+'/'+self.PLASMA_OUTPUT_FILE)

            for k,v, in self.plasma.items():
                print(('{0} : {1}'.format(k, v)))
            print(' ')
            
        except Exception as e:
            print(e)
            print('no PLASMA_OUTPUT_FILE defined in config file. ')
            sys.stdout.flush()
            
        if (plasmaExists == 'T'):
            if 'plasmaOutputDir' in self.plasma:
                print('\t PLASMA output will be read from {} \n'.format(self.plasma['plasmaOutputDir']+'_'+prj))

            ## OVERWRITE FTX PARAMETERS WITH INPUTS FROM PLASMA 
            
            if 'flux' in self.plasma:
                self.xp.parameters['flux']=self.plasma['flux']*1.0e-18
                print(('replaced flux in Xolotl (in ion/nm2s = 1e-18 ion/m2s) by value given by PLASMA {} (in ion/m2s) \n'.format(self.plasma['flux'])))
            else:
                print(('no flux specified in PLASMA, so using values in Xolotl {} \n'.format(self.xp.parameters['flux'])))

            #xp.parameters['tempModel']=='twoLine':
            if 'tempHandler' in self.xp.parameters or 'tempHandler' in self.plasma: 
                if 'heat' in self.plasma:
                    print('TEST: use heat defined by PLASMA')
                    self.xp.parameters['tempHandler']='heat'
                    if len(self.plasma['heat'])>1:
                        self.xp.parameters['tempParam']=[self.plasma['heat'][0]*1.0e-18, self.plasma['heat'][1]] #m2 -> /nm2
                    else:
                        self.xp.parameters['tempParam']=[self.plasma['heat']*1.0e-18, 300.0]
                    print('use heat given by PLASMA ', self.xp.parameters['tempParam'], '(here used in W/nm2s = 1e-18 W/m2s)')
                elif 'heat' in self.xp.parameters:
                    print('TEST: use heat defined by Xolotl')
                    self.xp.parameters['tempHandler']='heat'
                    if len(self.xp.parameters['heat'])>1:
                        self.xp.parameters['tempParam']=[self.xp.parameters['heat'][0], self.xp.parameters['heat'][1]]
                    else:
                        self.xp.parameters['tempParam']=[self.xp.parameters['heat'], 300.0]
                    print('use heat given by Xolotl ', self.xp.parameters['tempParam']) #if given by Xolotl, assume it's in Xolotl's units                    
                elif 'tempHandler' in self.plasma:
                    print('no heat flux provided. Use temperature model defined in PLASMAs input (assume tempParam is given too):')
                    print('\t tempHandler=', self.plasma['tempHandler'], ' and tempParam =', self.plasma['tempParam'])
                elif 'tempHandler' in self.xp.parameters:
                    print('no heat flux provided. Use temperature model defined in FTX config file (assume tempParam is given too):')
                    print('\t tempHandler=', self.xp.parameters['tempHandler'], ' and tempParam =', self.xp.parameters['tempParam'])                    
                else:
                    print('\t WARNING: no heat or other temperature model given in PLASMA or Xolotl. Use constant temp at 300K')
                    self.xp.parameters['tempHandler']='constant'
                    self.xp.parameters['tempParam']=300.0
                if 'startTemp' in self.xp.parameters:
                    del self.xp.parameters['startTemp']
                    print('\t and removed startTemp from xolotls parameters')


            #xp.parameters['tempModel']=='oneLine':
            #might need to check for more keywords if there're other one-line temp models 
            elif ('heat' in self.plasma) or ('heat' in self.xp.parameter) or ('startTemp' in self.plasma) or ('startTemp' in self.xp.parameters):
                if 'heat' in self.plasma:
                    if len(self.plasma['heat'])>1:
                        self.xp.parameters['heat']=[self.plasma['heat'][0]*1.0e-18,self.plasma['heat'][1]]
                    else: #no bulkT given in xolotl
                        self.xp.parameters['heat']=[self.plasma['heat'][0]*1.0e-18,300]
                    print('use heat flux given by PLASMA ', self.xp.parameters['heat'] ,' (here in W/nm2s = 1e-18 W/m2s)')
                    if 'startTemp' in self.xp.parameters:
                        del self.xp.parameters['startTemp']
                        print('\t and removed startTemp from xolotl parameters')
                elif 'heat' in self.xp.parameters:
                    print('no heat flux specified in PLASMA, so using values in Xolotl ', self.xp.parameters['heat'])
                    if 'startTemp' in self.xp.parameters:
                        del self.xp.parameters['startTemp']
                        print('\t and removed startTemp from xolotls parameters')
                elif 'startTemp' in self.plasma:
                    self.xp.parameters['startTemp']=self.plasma['startTemp']
                    print('no heat flux provided by PLASMA or Xolotl; use fixed temperature specified by PLASMA ', self.plasma['startTemp'])
                elif 'startTemp' in self.xp.parameters:
                    print('no heat flux provided by PLASMA or Xolotl; use fixed temperature specified by Xolotl (in config or default): ',self.xp.parameters['startTemp'])
                else:
                    print('no heat or temperature defined in PLASMAs output or Xolotls input. use default value T = 300K')
                    self.xp.parameters['startTemp'] = 300.0
            else:
                print('\t WARNING: no surface temperature / heat flux model defined. will assume constant surface temeprature at 300K')
                self.xp.parameters['startTemp']=300
                if 'heat' in self.xp.parameters:
                    del self.xp.parameters['heat']
                
            sys.stdout.flush()
            
        else:
            print('no PLASMA input file or paramters defined in config file')
            sys.stdout.flush()
            
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
                    print(('\t PLASMA output of angular distributions will be read from {} \n'.format(self.plasma[plasma_output_dir])))
                    self.angleFile.append(self.plasma[plasma_output_dir].strip()+'/'+self.PLASMA_ANGLE_FILE.strip()) #self.gitr['gitrOutputDir'].strip()
                    self.aWeightFile.append(self.plasma[plasma_output_dir].strip()+'/'+self.PLASMA_AWEIGHT_FILE.strip()) #self.gitr['gitrOutputDir'].strip()

                except Exception as e2:
                    print(e2)
                    print('no input from PLASMA')
                
                print(('\t reading angles and weights for {0} from {1} {2};\n'.format(prj, self.angleFile[i], self.aWeightFile[i])))
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
                    print(('\t PLASMA output of energy distributions will be read from {} \n'.format(self.plasma[plasma_output_dir])))
                    self.FT_energy_file_name.append(self.ft_energy_input_file[i]) #"He_W0001.ED1"
                    #where all the energy distribution files are located
                    self.eadist_output_path.append(self.plasma[plasma_output_dir].strip())
                    try:
                        self.PLASMA_EADIST_FILE 
                        self.eadist_output_file.append(self.PLASMA_EADIST_FILE)
                        print('\t using energy distribution file format given in PLASMA eadist file', self.PLASMA_EADIST_FILE)
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
                    
            #initialize maxRangeXolotl list
            self.maxRangeXolotl.append(0.0)


        #stage initial network File (INIT mode) OR restart files (RESTART mode)
        if (self.driverMode=='INIT'):
            if os.path.exists(self.INPUT_DIR+'/'+self.NETWORK_FILE):
                print(('\t INIT mode: stage initial network file {}\n'.format(self.NETWORK_FILE)))
                self.services.stage_input_files(self.NETWORK_FILE)
                print('\t \t ...initial network file staged succesfully {}\n')
            else:
                print(('\t WARNING: INIT mode: could not find initial network file {}\n'.format(self.NETWORK_FILE)))
                print(('\t WARNING: INIT mode: in input file directory {}\n'.format(self.INPUT_DIR)))
                print(('\t WARNING: INIT mode: SKIP staging  {}\n'.format(self.NETWORK_FILE)))
                

        elif (self.driverMode=='RESTART'):
            restart_files = self.services.get_config_param('RESTART_FILES') 
            print(('\t RESTART mode: stage restart files {} \n'.format(restart_files)))
            restart_list = restart_files.split()
            for index in range(len(restart_list)): 
                print('\t staging input file ', restart_list[index]) 
                self.services.stage_input_files(restart_list[index])
            print('\t \t ...restart files staged succesfully')

        sys.stdout.flush()
        self.services.update_plasma_state()


    def step(self, timeStamp=0.0,**keywords):

        cwd = self.services.get_working_dir()

        print('\n')
        print(('FT-X driver:step called with: ', keywords))
        print('\t output file of the FT-X workflow:')
        if 'LOG_FILE' in keywords:
            logFile=keywords['LOG_FILE']
            outFile=cwd+'/'+logFile
            print('\t \t log file defined in keywords: ', outFile)
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

        print(' ')
        print('xolotl-ftridyn_driver: step \n')


        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')
        
        self.services.stage_plasma_state() 

        #check that loop doesnt go over the end time
        time=self.time['INIT_TIME']
        end_time=self.time['END_TIME']

        if time >= end_time:
            print('init time ', time , ' >= end time', end_time)
            print("EXIT SIMULATION")
            sys.stdout.flush()
            return


        if time+self.time['LOOP_TIME_STEP']>end_time:
            self.time['LOOP_TIME_STEP']=end_time-time
            print(' ')
            print('\t WARNING: time step given in config file longer than needed for last loop ')
            print(('\t before starting time-loop, adapt driver time step to {} to reach exactly endTime '.format( self.time['LOOP_TIME_STEP'])))
            self.xp.parameters['petscArgs']['-start_stop']=self.time['LOOP_TIME_STEP']/10.0
            print(('\t accordingly, Xolotls data is saved every (start_stop) = {} '.format( self.xp.parameters['petscArgs']['-start_stop'])))
        else:
            print('\t before starting time-loop, checked that time step given in config file is not longer than needed to reach the end of the simulation')

        print('\n')
        sys.stdout.flush()

        

        #for time in numpy.arange(self.initTime,self.endTime,self.timeStep):
        while time<end_time: #self.driver['END_TIME']:

            self.services.stage_plasma_state()
            print(('driver time (in loop) {} \n'.format(time)))
            self.services.update_plasma_state()

            #keep all files to be saved (not plasma state) in folder with time stamp
            timeFolder='t'+str(time)
            if not os.path.exists(timeFolder):
                os.makedirs(timeFolder)
            print(('output of this time-loop will be saved in {} \n'.format(timeFolder)))

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
            print('F-TRIDYN:\n')
            # A) GET INPUT THAT MIGHT CHANGE EVERT LOOP READY

            #determine parameters related to init/restart
            iW=self.plasmaSpecies.index('W')
            if (self.driverMode == 'INIT'):
                print('\t init mode yes\n')
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
                print('\t init mode no \n') 
                self.ftridyn['iQ0']=-1
                self.ftridyn['nDataPts'] = translate_xolotl_to_ftridyn.xolotlToLay(totalDepth=self.ftridyn['totalDepth'],logFile=outFile)
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
                        
                print(('\t passing to F-Tridyn the list of targets t{} \n'.format(targetList)))

                #Xolotl only outputs He_W0001.LAY; but it's same substrate composition for running all projectiles
                iHe=self.plasmaSpecies.index('He')
                for i in range(len(self.plasmaSpecies)): #self.plasmaSpecies.iteritems():
                    prj=self.plasmaSpecies[i]
                    if prj!='He': #do not self-copy
                        print(('\t copy {0} as {1} '.format( self.ftx_lay_file[iHe], self.ftx_lay_file[i]))) 
                        shutil.copyfile(self.ftx_lay_file[iHe],self.ftx_lay_file[i])
                
                if (self.ftridyn['totalDepth']==0.0):
                    print('\t Totaldepth from last_TRIDYN.dat \n') 
                    self.ftridyn['nTT']=10*numpy.max(numpy.loadtxt('last_TRIDYN.dat')[:,0])
                else:
                    print(('\t totalDepth fixed to {} \n'.format(self.ftridyn['totalDepth'])))
                    self.ftridyn['nTT']=self.ftridyn['totalDepth']

            sys.stdout.flush()
            self.services.update_plasma_state()

            # B) RUN FTRIDYN

            print(' ')
            for i in range(len(self.plasmaSpecies)): #self.plasmaSpecies.iteritems():
                prj=self.plasmaSpecies[i]
                maxDepth=[]
                if (self.fluxFraction[i] > 0.0) and not (all(k==0 for k in self.weightAngle[i])):
                    print(('running F-Tridyn for {0} with flux fraction = {1}\n'.format(prj, self.fluxFraction[i])))
                    print(('\t and not all angle weights are zero; max angleWeight is {}\n'.format(max(self.weightAngle[i]))))
                    sys.stdout.flush()
                    
                    #component/method calls now include arguments (variables)
                    self.services.call(ftridyn, 'init', timeStamp, dTime=time, fPrj=prj, fTargetList=targetList, ftParameters=self.ftridyn , fEnergyIn=self.energyIn[i], fAngleIn=self.angleIn[i], fWeightAngle=self.weightAngle[i], ft_folder=self.FT_OUTPUT_FOLDER, input_file=self.ft_input_file[i], otherInFiles=[self.FT_SURFACE_FILE,self.ftx_lay_file[i]], energy_file_name=self.FT_energy_file_name[i], orig_energy_files_path=self.eadist_output_path[i], orig_energy_files_pattern=self.eadist_output_file[i], output_file=outFile)
                    sys.stdout.flush()
                    
                    self.services.call(ftridyn, 'step', timeStamp, ftParameters=self.ftridyn, fEnergyIn=self.energyIn[i], fAngleIn=self.angleIn[i], fWeightAngle=self.weightAngle[i], output_file=outFile)
                    sys.stdout.flush()
                    self.services.stage_plasma_state()


                    # C) POSTPROCESSING OF prj -> W

                    # 1) access ft folder (read file containing path to FT output directory, as outputPath=...):
                    t=param_handler.read(self.FT_OUTPUT_PWD_PATH)
                    for key,value in t.items():
                        self.ftridyn[key] = value

                    print(('from {0} , \t the path to output of FTRIDYN is {1} \n'.format(self.FT_OUTPUT_PWD_PATH, self.ftridyn['outputPath'])))

                    #2) #get maximum projectile range to ensure bins are added correctly in 'translate_ftridyn_to_xolotl'

                    ft_output_prj_file=self.ft_output_prj_file[i]
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
                        ft_output_file=self.ft_output_file[i]
                        #get implantation profile
                        #pass values as dictionary
                        ft_implProfiles_dictionary={}
                        ft_implProfiles_dictionary['ftridynOnePrjOutput']=ft_output_prj_file
                        ft_implProfiles_dictionary['ftridynFolder']=angleFolder
                        ft_implProfiles_dictionary['angle']=self.angleIn[i]
                        ft_implProfiles_dictionary['weightAngle']=self.weightAngle[i]
                        ft_implProfiles_dictionary['prjRange']=maxRange
                        #different grid keywords depending on Xolotl version:
                        if 'grid' in self.xp.parameters:
                            ft_implProfiles_dictionary['nBins']=self.xp.parameters['grid'][0]
                        elif 'gridParam' in  self.xp.parameters:
                            ft_implProfiles_dictionary['nBins']=self.xp.parameters['gridParam'][0]
                        else:
                            ft_implProfiles_dictionary['nBins']=200
                            print('\t WARNING: no grid or gridParam provided; assume nBins=200')
                                
                        ft_implProfiles_dictionary['logFile']=outFile


                        pkl_impl_file=cwd+'/ft_implProfiles.pkl'  ## define name in config file, here give abs path
                        pickle.dump(ft_implProfiles_dictionary, open(pkl_impl_file, "wb" ) )

                        try:
                            self.TRANSLATE_FT2XOL
                            ft_implProfile_script=self.TRANSLATE_FT2XOL
                            print('Launch task of running python ', ft_implProfile_script)
                            
                        except: #DEFAULT PATH: 
                            ft_implProfile_script = self.BIN_PATH+'/ips-iterative-xolotlFT/python_scripts_for_coupling/translate_ftridyn_to_xolotl.py'
                            print('Launch task of running default python script ', ft_implProfile_script)

                        sys.stdout.flush()

                        task_id_impl = self.services.launch_task(1,self.services.get_working_dir(),
                                                            'python', ft_implProfile_script, logfile='tridynPlotting.log')
                        ret_val_impl = self.services.wait_task(task_id_impl)
                        sys.stdout.flush()

                    else: #if len(maxDepth)==0
                        print("nothing was implanted. maxRange and profile = 0 ")
                        maxRange=0.0
                        self.maxRangeXolotl[i]=0.0
                        #tridyn.dat = zeros
                        outputFTFile=open(self.FT_OUTPUT_PROFILE_TEMP, "w")
                        outputFTFile.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ")
                        outputFTFile.close()
                        sys.stdout.flush()
                        #print "END OF THIS SIMULATION" & return


                    #3) get the sputtering yield (or use fixed value)                   
                    print(' ')
                    #pass values as dictionary
                    ft_getYields_dictionary={}
                    ft_getYields_dictionary['ftridynOneOutOutput']=ft_output_file
                    ft_getYields_dictionary['ftridynFolder']=angleFolder
                    ft_getYields_dictionary['angle']=self.angleIn[i]
                    ft_getYields_dictionary['weightAngle']=self.weightAngle[i]
                    ft_getYields_dictionary['fNImpacts']=self.ftridyn['nImpacts']
                    ft_getYields_dictionary['logFile']=outFile
                    
                    pkl_gy_file=cwd+'/ft_getYields.pkl'  ## define name in config file, here give abs path
                    #we need to close this file later (to open it again and read yields), so use alternative to
                    #pickle.dump(ft_getYields_dictionary, open(pkl_gy_file, "wb" ) )
                    with open(pkl_gy_file, "wb") as pf:
                        pickle.dump(ft_getYields_dictionary, pf)
                    pf.close()
                    sys.stdout.flush()
                    
                    try:
                        self.GET_YIELDS
                        ft_getYields_script=self.GET_YIELDS
                        print('Launch task of running python ', ft_getYields_script)
                        
                    except: #DEFAULT PATH:
                        ft_getYields_script = self.BIN_PATH+'/ips-iterative-xolotlFT/python_scripts_for_coupling/get_yields.py'
                        print('Launch task of running default python script ', ft_getYields_script)

                    sys.stdout.flush()

                    task_id_gy = self.services.launch_task(1,self.services.get_working_dir(),
                                                        'python', ft_getYields_script, logfile='get_yields.log')
                    ret_val_gy = self.services.wait_task(task_id_gy)

                    if os.path.exists(pkl_gy_file):
                        with open(pkl_gy_file, "rb") as pf:
                            getYields_dic = pickle.load(pf)
                            yields=getYields_dic['yields']
                        pf.close()
                        print('reading the pickle file, get_yields returned [total SpY, total RY] = ', yields)
                    else:
                        print('WARNING! could not read yield from get_yields pickle file. Set to zero')
                        yields=[0.0, 0.0]
                    
                    #overwrite spY value if mode is 'calculate'
                    if self.spYieldMode[i]=='calculate':
                        self.spYield[i]=float(yields[0])
                    if self.rYieldMode[i]=='calculate':
                        self.rYield[i]=float(yields[1])

                    #4) save tridyn.dat
                    #append output to allTridynNN.dat for each species (and save to what folder?)        
                    ft_output_profile_final=self.FT_OUTPUT_PROFILE_FINAL+'_'+prj
                    tempfile = open(self.FT_OUTPUT_PROFILE_TEMP,"r")
                    f = open(ft_output_profile_final, "a")                    
                    f.write('%s%s \n' %(tempfile.read().rstrip('\n'),self.maxRangeXolotl[i]))                    
                    f.close()
                    tempfile.close()
                    
                    #keep copies of tridyn.dat
                    ft_output_profile_temp_prj=self.FT_OUTPUT_PROFILE_TEMP+'_'+prj #for each species
                    shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP,ft_output_profile_temp_prj)
                    shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, timeFolder+'/'+ft_output_profile_temp_prj)
            
                
                    #5) MOVE FOLDERS TO DIRECTORY WITH TIME-STAMP & RENAME FOR (Tg,Prj) SPECIES  
                
                    shutil.move(self.ftridyn['outputPath']+'/'+self.FT_OUTPUT_FOLDER,timeFolder+'/'+self.FT_OUTPUT_FOLDER+'_'+prj+'W')                
                    self.services.update_plasma_state()
                    print(('\t done with F-TRIDYN for {} '.format(prj)))
                    print('\n')
                    sys.stdout.flush()

                #if flux fraction == 0 or all weight angles == 0.0:
                #if (self.gitr['fluxFraction'][i] > 0.0) and not (all(k==0 for k in self.weightAngle[i])):
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

                    #keep copies of tridyn.dat    
                    ft_output_profile_temp_prj=self.FT_OUTPUT_PROFILE_TEMP+'_'+prj #for each species
                    shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP,ft_output_profile_temp_prj)
                    shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, timeFolder+'/'+ft_output_profile_temp_prj)

                sys.stdout.flush()
            #end of for loop:

            
            ######species independent ############

            #6) write sputtering yields to file so they can be used by Xolotl

            yieldString=str(time)
            print('Sputtering and Reflection Yields due to:')
            
            for i in range(len(self.plasmaSpecies)): #self.plasmaSpecies.iteritems():  
                prj=self.plasmaSpecies[i]
                print(('\t{0} :  spY = {1} and rY = {2} '.format(prj,self.spYield[i], self.rYield[i])))
                prjYieldsString=prj+' ' +str(self.spYield[i])+' '+str(self.rYield[i])
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
            
            print('self.maxRangeXolotl = ', self.maxRangeXolotl ,' and max of self.maxRangeXolotl =', max(self.maxRangeXolotl))
            print('spYield = ', self.spYield, ' and max(self.spYield) = ', max(self.spYield))
            
            if max(self.maxRangeXolotl)==0 and max(self.spYield)==0:
                print("nothing was implanted or sputtered")
                print("likely all weights are zero")
                print("END OF THIS SIMULATION")
                sys.stdout.flush()
                return

            #7) write format tridyn.dat to include W redep in Xolotl
            
            print(' ')
            if os.path.exists(self.FT_OUTPUT_PROFILE_TEMP):
                os.remove(self.FT_OUTPUT_PROFILE_TEMP)
            combinedFile = open(self.FT_OUTPUT_PROFILE_TEMP, 'a')

            ##New attempt at checking which species should be written into tridyn.dat:
            ##    - first check for flux fraction ;
            ##             if 0, skip writing into tridyn.dat because there's no information on impact energy and angle (no FT)
            ##                   regardless of network: lines missing for species that exist in the network does NOT cause any issues (just nothing implanted)
            ##             if >0, check network --> 
            ##    - use the network param values when available,
            ##    - assume the standard (species included if fluxFraction>0, not included if fluxFraction<0) if the network isn't given explicitely
            ##      I.e., now we try to merge two checks: it works in cases where more species passed by GITR, but not handled by Xolotl; or if fluxFraction=0.
            
            for i in range(len(self.plasmaSpecies)):
                prj=self.plasmaSpecies[i]                
                ft_output_profile_temp_prj=timeFolder+'/'+self.FT_OUTPUT_PROFILE_TEMP+'_'+prj
                profile=open(ft_output_profile_temp_prj, "r")
                tridynString=profile.read().rstrip('\n')
                combinedTridynString=str(tridynString)+str(self.maxRangeXolotl[i])
                print(('for {0}, fraction in plasma = {1} , and reflection = {2} '.format(prj,self.fluxFraction[i], self.rYield[i])))
                print(('\t effective fraction (in plasma * (1-reflection)) = {} '.format(self.fluxFraction[i]*(1-self.rYield[i]))))
                sys.stdout.flush()

                if ( (self.fluxFraction[i] > 0)):
                    #the format of tridyn.dat is different for the pulsed & UQ executables of Xolotl:
                    if prj!='W': #then the name in the tridyn.dat line is the same as prj
                        ##if He,  check He's position in netParam, i.e., index i=0  
                        if prj=='He':
                            if ('netParam' in self.xp.parameters):
                                if (self.xp.parameters['netParam'][i]==0):
                                    print('\t Xolotl netowrk exists for ' , prj, 'given in plasmaSpecies, but entry in netParam is zero ; will skip in tridyn.dat')
                                else:
                                    print('\t For ' , prj , 'netparam = ' ,self.xp.parameters['netParam'][i] , ' is used in Xolotl ; write line for ', prj , ' in tridyn.dat')
                                    combinedFile.write("%s %s %s\n" %(prj,str(1),str(self.fluxFraction[i]*(1-self.rYield[i])))) 
                                    combinedFile.write("%s\n" %(combinedTridynString))
                            else:
                                print('\t WARNING: netparam not given in Xolotl ; write line for ', prj , ' in tridyn.dat')
                                print('\t \t this might give an ERROR if species isnt part of Xolotls network')
                                combinedFile.write("%s %s %s\n" %(prj,str(1),str(self.fluxFraction[i]*(1-self.rYield[i]))))
                                combinedFile.write("%s\n" %(combinedTridynString))
                        ##if D or T,  check position in netParam, i.e., index i=2,3 -> i-1 = 1 or 2 (no W in netParam)
                        elif prj=='D' or prj=='T':
                            if ('netParam' in self.xp.parameters):
                                if (self.xp.parameters['netParam'][i-1]==0):
                                    print('\t Xolotl network exists for ' , prj, 'given in plasmaSpecies, but entry in netParam is zero ; will skip in tridyn.dat')
                                else:                                    
                                    print('\t For ' , prj , 'netparam = ' ,self.xp.parameters['netParam'][i-1] , ' is used in Xolotl ; write line for ', prj , ' in tridyn.dat')
                                    combinedFile.write("%s %s %s\n" %(prj,str(1),str(self.fluxFraction[i]*(1-self.rYield[i])))) 
                                    combinedFile.write("%s\n" %(combinedTridynString))
                            else:
                                print('\t WARNING: netparam not given in Xolotl ; write line for ', prj , ' in tridyn.dat')
                                print('\t \t this might give an ERROR if species isnt part of Xolotls network')
                                combinedFile.write("%s %s %s\n" %(prj,str(1),str(self.fluxFraction[i]*(1-self.rYield[i]))))
                                combinedFile.write("%s\n" %(combinedTridynString))
                        ##W is called interstitial in tridyn.dat, i.e., "I" & always exists in the network (no need to check) 
                    elif prj=='W':
                        print('\t WARNING: the wrapper assumes that interstitials always exists in Xolotls network')
                        print('\t \t it will write line for ', prj , ' in tridyn.dat')
                        print('\t \t this might give an ERROR if interstitials arent part of Xolotls network')
                        combinedFile.write("%s %s %s\n" %('I',str(1),str(self.fluxFraction[i]*(1-self.rYield[i]))))
                        combinedFile.write("%s\n" %(combinedTridynString))
                    else:
                        print('\t WARNING: species ', prj, 'cannot be handled by Xolotl yet.')
                        print('\t \t it has been used so far (for spY, etc), but will skip writing into tridyn.dat')

                else:
                    print('\t Flux fraction for ', prj, 'is zero')
                    print('\t \t for now, skip writing anything, even is prj exists in network (no checks in place)')
                    print('\t \t Xolotl will run, with no ', prj, ' implanted')
                    
            profile.close()
            sys.stdout.flush()
            combinedFile.close()

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

            sys.stdout.flush()
            self.services.update_plasma_state()

            ######################################
            ############## RUN XOLOTL ############ 
            ###################################### 

            #make sure at least one of the species was implanted:
            if all(k==0 for k in self.maxRangeXolotl):
                print("WARNING: none of the species was implanted")
                print("EXIT SIMULATION")
                sys.stdout.flush()
                return

            #if something was implanted:
            #Xolotl parameter modifications that need to be done at every loop:

            #list modules right before running Xolotl:
            modList=os.system('module list')
            print('TEST: we are running Xolotl with the following modules loaded: ', modList)

            print('load cray-hdf5-parallel')
            os.system('module load cray-hdf5-parallel')

            modList=os.system('module list')
            print('TEST: we are running Xolotl with the following modules loaded: ', modList)
            sys.stdout.flush()
            
            #calculate effective sputtering yield; i.e., weighted by relative flux of W-to-He
            totalSpYield=0
            for i in range(len(self.plasmaSpecies)): #self.plasmaSpecies.iteritems():
                prj=self.plasmaSpecies[i]
                print(('contribution of {0} to total sputtering yield = {1} '.format( prj, float(self.fluxFraction[i])*float(self.spYield[i]))))
                totalSpYield+=(float(self.fluxFraction[i])*float(self.spYield[i])) #self.fluxFraction[i]
            print(('total weighted sputtering yield = {} (passed to Xolotl)\n'.format(totalSpYield)))
            self.xp.parameters['sputtering'] = totalSpYield            
            
            #time and time-step related parameters
            self.xp.parameters['petscArgs']['-ts_final_time']=time+self.time['LOOP_TIME_STEP']
            print('\n')
            print('XOLOTL: ')
            print(('\t Run from t = {}'.format(time)))
            print(('\t to t = {}'.format(self.xp.parameters['petscArgs']['-ts_final_time'])))
            print(('\t and time-step = {} '.format( self.time['LOOP_TIME_STEP'])))

            if self.driverMode == 'INIT':
                if os.path.exists(self.INPUT_DIR+'/'+self.NETWORK_FILE):
                    print('\t init mode: modify xolotl parameters that might change at every loop, and load networkFile file \n')                
                    self.xp.parameters['networkFile'] = self.XOLOTL_NETWORK_FILE                    
                else:
                    #no network file in input dir -- do not add to param dictionary, so it's not in param file and doesn't try to load
                    print('\t init mode: modify xolotl parameters that might change at every loop \n')
                    print('\t \t WARNING: no network file in input dir; will create (not load) the network \n')
            elif self.driverMode == 'RESTART':
                #add (or replace) networkFile line to parameter file
                print('\t restart mode: modify xolotl parameters that might change at every loop, including adding the networkFile \n')
                self.xp.parameters['networkFile'] = self.XOLOTL_NETWORK_FILE
                if 'netParam' in self.xp.parameters:
                    print('\t \t and delete netParam from xolotl parameters (i.e., from the parameter file)')
                    del self.xp.parameters['netParam']
                else:
                    print('\t \t netParam does not exist in the xolotl parameters. No need to delete it')
                if 'grouping' in self.xp.parameters:
                    print('\t \t and delete grouping from xolotl parameters (i.e., from the parameter file)')
                    del self.xp.parameters['grouping']
                else:
                    print('\t \t grouping does not exist in the xolotl parameters. No need to delete it')

                    
            #determine if he_conc true/false ; if true, add '-he_conc' to petsc arguments 
            if self.driver['XOLOTL_HE_CONC']=='Last':
                if time+1.5*self.time['LOOP_TIME_STEP']>end_time: #self.driver['END_TIME']:  #*1.5, to give marging of error
                    self.petsc_heConc=True
                    print('printing He concentrations in the last loop')
                elif time<(end_time-self.time['LOOP_TIME_STEP']): #self.driver['END_TIME']-self.driver['LOOP_TIME_STEP']):
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

            n_overgrid_loops=0
            
            while self.xolotlExitStatus=='collapsed' or self.xolotlExitStatus=='overgrid':

                if self.xolotlExitStatus=='collapsed':
                    self.collapsedLoops+=1
                    print(' ')

                #set a maximum number of tries
                if self.collapsedLoops<=int(self.driver['MAX_COLLAPSE_LOOPS']):

                    self.services.call(xolotl, 'init', timeStamp, dTime=time, xParameters=self.xp.parameters, output_file=outFile) #, xFtCoupling=self.driver['FTX_COUPLING'])
                    self.services.call(xolotl, 'step', timeStamp, dTime=time, xHe_conc=self.petsc_heConc, xParameters=self.xp.parameters, output_file=outFile, dZipOutput=self.driver['ZIP_XOLOTL_OUTPUT'], n_overgrid_loops=n_overgrid_loops)

                    sys.stdout.flush()
                    self.services.stage_plasma_state()

                    statusFile=open(self.XOLOTL_EXIT_STATUS, "r")
                    self.xolotlExitStatus=statusFile.read().rstrip('\n')

                    print(('\t Xolotl ended simulation with status {} \n'.format(self.xolotlExitStatus)))

                    if self.xolotlExitStatus=='good':
                        print(('\t Xolotl successfully executed after {} tries'.format(self.collapsedLoops)))
                        print('continue with IPS simulation \n')

                    elif self.xolotlExitStatus=='diverged':
                        print('\t ERROR: XOLOTL SOLVER DIVERGED ')
                        print('END IPS SIMULATION \n')
                        quit()
                    elif self.xolotlExitStatus=='overgrid':
                        n_overgrid_loops+=1
                        print('\t WARNING: XOLOTL OVERGRID ')
                        print('\t RUNNING transgerGrid...')
                        #xolotlStop already copied as _overgrid within xolotl_comp; no need to copy here again

                        #also need to update the values of grid and voidPortion to match what's in the network file:
                        #Most likely this will only work if using 'grid', not 'gridParam' in network file!!
                        shutil.copyfile(self.xp.parameters['networkFile'],self.xp.parameters['networkFile']+'_overgrid_'+str(n_overgrid_loops))
                        shutil.move('xolotlStop.h5', self.xp.parameters['networkFile'])
                        [newGridSize, newVoidPortion] = transferGrid.transferGrid(self.xp.parameters['networkFile'])
                        if os.path.exists(self.xp.parameters['networkFile']):
                            #make sure xolotlStop exists
                            shutil.copyfile(self.xp.parameters['networkFile'],'xolotlStop.h5')
                        else:
                            print('\t WARNING: networkFile does not exists after transferGrid')
                            print('\t keep old one')
                            shutil.copyfile(self.xp.parameters['networkFile']+'_overgrid_'+str(n_overgrid_loops),self.xp.parameters['networkFile'])
                        print('\t in file ', self.xp.parameters['networkFile'])    
                        print('\t \t updated the values of grid to ', newGridSize)
                        print('\t \t updated the values of voidPortion to ', newVoidPortion)
                        print('\t DONE RUNNING transferGrid!')
                        sys.stdout.flush()
                        self.xp.parameters['grid'][0] = newGridSize
                        self.xp.parameters['voidPortion'] = newVoidPortion
                        self.services.update_plasma_state()
                        
                    elif self.xolotlExitStatus=='collapsed':
                        print('\t WARNING: simulation exited loop with status collapse')
                        print(('\t try number {0} out of {1}\n'.format(self.collapsedLoops,self.maxCollapseLoops)))                    

                    else:
                        print('\t WARNING: Xolotl exit status UNKOWN -- IPS simulation continues \n')

                else: #reached maximum number of tries for collapsing time steps
                    print('\t ERROR: reached maximum number of tries for collapsing time steps without a successful run')
                    print('END IPS SIMULATION \n')
                    quit()

            #UPDATE: Xolotl generates HDF5 file of TRIDYN, copied as 'last_TRIDYN_toBin.h5'
            #thus no need to copy it; and binTRIDYN will transform it to text file, 'last_TRIDYN.dat' 
            print('bin Xolotls output')
            binTRIDYN.binTridyn()
            print('...succesfull')
            
            #store xolotls profile output for each loop (not plasma state)          
            #UPDATE: 'toBin' is .h5 format instead of .dat
            currentXolotlOutputFileToBin='last_TRIDYN_toBin_%f.h5' %time
            shutil.copyfile('last_TRIDYN_toBin.h5', currentXolotlOutputFileToBin)
            currentXolotlOutputFile='last_TRIDYN_%f.dat' %time
            shutil.copyfile('last_TRIDYN.dat', currentXolotlOutputFile)

            #copy last_TRIDYN.dat and netowrk file to input dir as well:
            print('coy last_TRIDYN.dat and ', self.xp.parameters['networkFile'], 'to ', self.INPUT_DIR )
            shutil.copyfile('last_TRIDYN.dat', self.INPUT_DIR+'/last_TRIDYN.dat')
            shutil.copyfile(self.xp.parameters['networkFile'], self.INPUT_DIR+'/'+self.xp.parameters['networkFile'])
            
            #append output:
            #retention
            tempfileRet = open(self.XOLOTL_RETENTION_TEMP,"r")
            fRet = open(self.XOLOTL_RETENTION_FINAL, "a")
            fRet.write(tempfileRet.read())
            fRet.close()
            tempfileRet.close()
            
            #surface
            tempfileSurf = open(self.XOLOTL_SURFACE_TEMP,"r")
            fSurf = open(self.XOLOTL_SURFACE_FINAL, "a")
            fSurf.write(tempfileSurf.read())
            fSurf.close()
            tempfileSurf.close()

            #save network file with a different name to use in the next time step
            currentXolotlNetworkFile='xolotlStop_%f.h5' %time
            shutil.copyfile('xolotlStop.h5',currentXolotlNetworkFile)
            shutil.copyfile('xolotlStop.h5',self.XOLOTL_NETWORK_FILE)

            #update driver mode after the 1st loop, from INIT to RESTART
            if self.driverMode == 'INIT':
                self.driverMode = 'RESTART'
                print(('switched driverMode to {} \n'.format(self.driverMode)))

            if self.driverMode != 'NEUTRAL':
                self.driverMode = 'NEUTRAL'
                print(('switched startMode to {} \n'.format(self.driverMode)))

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
                print(('\t driver time step (({0}) ane start_stop ({1}) unchanged (factor=1) \n'.format( self.time['LOOP_TIME_STEP'], self.xp.parameters['petscArgs']['-start_stop'])))

            if time+self.time['LOOP_TIME_STEP']>end_time: 
                self.time['LOOP_TIME_STEP']=end_time-time 
                print(' ')
                print('\t time step longer than needed for last loop ')
                print(('\t adapting driver time step to {} to reach exactly endTime '.format(self.time['LOOP_TIME_STEP'])))
                self.xp.parameters['petscArgs']['-start_stop']=self.time['LOOP_TIME_STEP']/10.0
                print(('\t and Xolotls data is saved every (start_stop) = {} '.format( self.xp.parameters['petscArgs']['-start_stop'])))

            print('\n')

            if self.time['XOLOTL_MAXTS_FACTOR'] != 1:
                if (self.time['LOOP_N']%self.driver['XOLOTL_MAXTS_NLOOPS']==0):
                    print(('\t change in Xolotls maximum time step after loop {} '.format( self.time['LOOP_N'])))
                    self.xp.parameters['petscArgs']['-ts_adapt_dt_max']*=self.time['XOLOTL_MAXTS_FACTOR']
                    print(('\t multiply time step by {0}, for a new time step = {1} \n'.format(self.driver['XOLOTL_MAXTS_FACTOR'] , self.xp.parameters['petscArgs']['-ts_adapt_dt_max'])))
                    if self.xp.parameters['petscArgs']['-ts_adapt_dt_max'] > self.time['XOLOTL_MAX_TS']:
                        self.xp.parameters['petscArgs']['-ts_adapt_dt_max']=self.time['XOLOTL_MAX_TS']
                        print(('\t Xolotls time-step reached the maximum allowed; set to limit value, {} \n'.format(self.driver['XOLOTL_MAX_TS'])))
                else:
                    print(( '\t continue with xolotl dt max {} \n'.format(self.xp.parameters['petscArgs']['-ts_adapt_dt_max'])))
            else:
                print(('\t Xolotls max time step ({}) unchanged (factor=1) \n'.format(self.xp.parameters['petscArgs']['-ts_adapt_dt_max'])))

            print(' ')
            sys.stdout.flush()
            self.services.update_plasma_state()

    def finalize(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: finalize')
        #can we add compressing output here? e.g., last_TRIDYN, xolotlStop...
        #and remove large output files? e.g., FTRIDYN.zip

        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')
        self.services.call(ftridyn, 'finalize', timeStamp)
        self.services.call(xolotl, 'finalize', timeStamp)
