#! /usr/bin/env python

from  component import Component
import sys
import os
import os.path
import subprocess
import numpy
import shutil
import translate_xolotl_to_ftridyn
import translate_ftridyn_to_xolotl
import get_yields
import binTRIDYN
import param_handler
import traceback

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
                #logFile = 'log.stdOut'
                #self.outFile=cwd+'/'+logFile
                

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
        
        print('reading DRIVER parameters from ips config file: \n')
        
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

        self.driverMode=self.driver['START_MODE']
        print('\n')
        print(('running IPS from t = {0} , to t = {1} , in steps of dt = {2}'.format(self.driver['INIT_TIME'], self.driver['END_TIME'], self.\
driver['LOOP_TIME_STEP'])))


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
            else:                
                print(('\t reading string input parameter {0} = {1}'.format(k, v )))
                self.xp.parameters[k]=v
                
        print('\n')
        print('replacing PETSC arguments from ips config file:\n')
        for k,v in self.XOLOTL_INPUT_PETSCARGS.items():
            if param_handler.is_int(v):
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

        #print "THIS IS A TEST: final Xolotl parameters: "
        #print self.xp.parameters


        #INCLUDE IF NEEDED (already debugged)
        #explicitly delete 'grouping' from xolotl params if grouping is OFF (groupHeV > MaxVSize)

        #if 'grouping' in self.xp.parameters:
        #    print 'this is a test: check if we should keep grouping'

        #    netParamList=self.xp.parameters['netParam'].split()      #maxVSize=netParamList[3]
        #    print 'this is test: netParamList is ', netParamList        
        #    groupingList=self.xp.parameters['grouping'].split()      #groupHeV=groupingList[0]
        #    print 'this is test: groupingList is ', groupingList

        #    if float(groupingList[0]) > float(netParamList[3]):
        #        del self.xp.parameters['grouping']
        #        print 'this is a test: grouping deleted as groupHeV=', groupingList[0], ' and maxVSize=', netParamList[3]
        #else:
        #    print 'this is a test: no grouping exists in xp.parameters dictionary'



        ### GITR RELATED PARAMETERS ###

        print(('read GITR parameters (from file) : {}\n'.format(self.INPUT_DIR+'/'+self.GITR_OUTPUT_FILE)))
        self.gitr = {} #xolotl_param_handler.xolotl_params()
        self.gitr=param_handler.read(self.INPUT_DIR+'/'+self.GITR_OUTPUT_FILE)
        #print('read GITR parameters (from file) : {}\n'.format(self.INPUT_DIR+'/'+self.GITR_OUTPUT_FILE))

        for k,v, in self.gitr.items():
            print(('{0} : {1}'.format(k, v)))
        print(' ')

        for k,v in self.GITR_INPUT_PARAMETERS.items(): #self.gitr.parameters.iteritems():
            if param_handler.is_int(v):
                print(( '\t reading integer GITR input parameter {0} = {1}'.format( k, v)))
                self.gitr[k]=int(v)
            elif param_handler.is_float(v):
                print(( '\t reading float GITR input parameter {0} = {1}'.format( k, v)))
                self.gitr[k]=float(v)
            elif param_handler.is_list(v):
                values = []
                for val in v.split(' '):
                    if param_handler.is_int(val):
                        values.append(int(val))
                    elif param_handler.is_float(val):
                        values.append(float(val))
                    else:
                        values.append(val)
                print(('\t reading list GITR input parameter {0} = {1}'.format( k, values)))
                self.gitr[k]=values                
            else:
                print(('\t reading string GITR input parameter {0} = {1} '.format( k, v))) 
                self.gitr[k]=v
                
        #if 'gitrOutputDir' in self.gitr:
        #    print('\t GITR output will be read from {} \n'.format(self.gitr['gitrOutputDir']+'_'+prj))

        ##NOT SURE IF THERE ARE MORE XOLOTL OR FT PARAMETERS THAT NEED TO OVERWRITTEN
        if 'flux' in self.gitr:
            self.xp.parameters['flux']=self.gitr['flux']
            print(('replaced flux in Xolotl by value given by GITR {}\n'.format(self.gitr['flux'])))  
        else:
            print(('no flux specified in GITR, so using values in Xolotl {} \n'.format(self.xp.parameters['flux'])))

        if 'heat' in self.gitr:
            self.xp.parameters['heat']=self.gitr['heat']
            print(('replaced heat flux in Xolotl by value given by GITR {}\n'.format(self.gitr['heat'])))
        elif 'heat' in self.xp.parameters:
            print(('no heat flux specified in GITR, so using values in Xolotl {} \n'.format(self.xp.parameters['heat'])))
        elif 'startTemp' in self.gitr:
            self.xp.parameters['startTemp']=self.gitr['startTemp']
            print(('no heat flux provided by GITR or Xolotl; use fixed temperature specified by GITR {}\n'.format(self.gitr['startTemp'])))
        else:
            print(('no heat flux provided by GITR or Xolotl; use fixed temperature specified by Xolotl {}\n'.format(self.xp.parameters['startTemp'])))

        #### FTRIDYN PARAMETERS ##### 
        ##LOOP OVER LIST OF PLASMA SPECIES SPECIES #########
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
        self.GITR_eadist_output_path=[]
        self.GITR_eadist_output_file=[]
        

        inputEnergy=self.gitr['inputEnergy']#.split(' ')
        inputAngle=self.gitr['inputAngle']#.split(' ')
        inputSpYield=self.ftridyn['inputSpYield'].split(' ')
        inputRYield=self.ftridyn['inputRYield'].split(' ')
        #inputFluxFraction=self.GITR_INPUT_PARAMETERS['inputFluxFraction'].split(' ')        

        for i in range(len(self.gitr['plasmaSpecies'])): #self.plasmaSpecies.iteritems():
            prj=self.gitr['plasmaSpecies'][i]

            self.energyIn.append(float(inputEnergy[i]))
            self.inAngle.append(float(inputAngle[i]))
            self.spYield.append(float(inputSpYield[i]))
            self.rYield.append(float(inputRYield[i]))
            #self.fluxFraction.append(float(inputFluxFraction[i]))
            print(('\t index {0}, species {1} '.format( i, self.gitr['plasmaSpecies'][i])))
            print(('\t energy {0}, angle {1} '.format( self.energyIn[i] ,self.inAngle[i])))
            print(('\t spYield {0} , rYield {1} , fluxFraction {2} \n'.format(self.spYield[i], self.rYield[i], self.gitr['fluxFraction'][i] )))

            if self.inAngle[i] < 0 :
                ##ADD +'_'+prj in the middle if the full path name for ITER cases, as multiple plasma species will follow distributions
                gitr_output_dir='gitrOutputDir'+'_'+prj
                print(('\t GITR output of angular distributions will be read from {} \n'.format(self.gitr[gitr_output_dir])))
                self.angleFile.append(self.gitr[gitr_output_dir].strip()+'/'+self.GITR_ANGLE_FILE.strip()) #self.gitr['gitrOutputDir'].strip()
                self.aWeightFile.append(self.gitr[gitr_output_dir].strip()+'/'+self.GITR_AWEIGHT_FILE.strip()) #self.gitr['gitrOutputDir'].strip()
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

            if self.energyIn[i] < 0:
                gitr_output_dir='gitrOutputDir'+'_'+prj
                print(('\t GITR output of energy distributions will be read from {} \n'.format(self.gitr[gitr_output_dir])))
                ##ADD +'_'+prj in the middle if the full path name for ITER cases, as multiple plasma species will follow distributions
                self.FT_energy_file_name.append(self.ft_energy_input_file[i]) #"He_W0001.ED1"
                #where all the energy distribution files are located
                #self.GITR_eadist_output_path.append(self.gitr['gitrOutputDir'].strip())#+'_'+prj)
                self.GITR_eadist_output_path.append(self.gitr[gitr_output_dir].strip())
                self.GITR_eadist_output_file.append(['dist','.dat'])#(self.GITR_EADIST_FILE)
            else:
                self.FT_energy_file_name.append('')
                self.GITR_eadist_output_path.append('')
                self.GITR_eadist_output_file.append([' ',' '])#('')

            #initialize maxRangeXolotl list
            self.maxRangeXolotl.append(0.0)


        #stage initial network File (INIT mode) OR restart files (RESTART mode)
        if (self.driver['START_MODE']=='INIT'):
            #print 'check if theres a network file in the input directory'
            #print self.INPUT_DIR+'/'+self.NETWORK_FILE
            #print os.path.exists(self.INPUT_DIR+'/'+self.NETWORK_FILE)
            if os.path.exists(self.INPUT_DIR+'/'+self.NETWORK_FILE):
                print(('\t INIT mode: stage initial network file {}\n'.format(self.NETWORK_FILE)))
                self.services.stage_input_files(self.NETWORK_FILE)
                print('\t \t ...initial network file staged succesfully {}\n')
            else:
                print(('\t WARNING: INIT mode: could not find initial network file {}\n'.format(self.NETWORK_FILE)))
                print(('\t WARNING: INIT mode: in input file directory {}\n'.format(self.INPUT_DIR)))
                print(('\t WARNING: INIT mode: SKIP staging  {}\n'.format(self.NETWORK_FILE)))
                

        elif (self.driver['START_MODE']=='RESTART'):
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
                #logFile = 'log.stdOut'
                #self.outFile=cwd+'/'+logFile
                #outF=open(self.outFile , 'a')

        print(' ')
        print('xolotl-ftridyn_driver: step \n')


        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')
        
        self.services.stage_plasma_state() 

        time=self.driver['INIT_TIME']
        sys.stdout.flush()

        #for time in numpy.arange(self.initTime,self.endTime,self.timeStep):
        while time<self.driver['END_TIME']:

            self.services.stage_plasma_state()
            print(('driver time (in loop) {} \n'.format(time)))
            self.services.update_plasma_state()

            #keep all files to be saved (not plasma state) in folder with time stamp
            timeFolder='t'+str(time)
            if not os.path.exists(timeFolder):
                os.makedirs(timeFolder)
            print(('output of this time-loop will be saved in {} \n'.format(timeFolder)))

            self.driver['LOOP_N']+=1
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
            iW=self.gitr['plasmaSpecies'].index('W')
            if (self.driverMode == 'INIT'):
                print('\t init mode yes\n')
                self.ftridyn['iQ0']=0

                targetList=[]
                targetList.append(self.gitr['plasmaSpecies'][iW]) #only W in the first loop
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
                targetList.append(self.gitr['plasmaSpecies'][iW])
                for prj in ['He', 'D', 'T']: #generate_ftridyn_input expects max 4 target species; declare all, even if empty
                    if prj in self.gitr['plasmaSpecies']:
                        i=self.gitr['plasmaSpecies'].index(prj)
                        if self.gitr['fluxFraction'][i]>0.0: #species exists and fraction > 0
                            targetList.append(prj)
                        else:                        
                            targetList.append('') #leave empty
                print(('\t passing to F-Tridyn the list of targets t{} \n'.format(targetList)))

                #Xolotl only outputs He_W0001.LAY; but it's same substrate composition for running all projectiles
                iHe=self.gitr['plasmaSpecies'].index('He')
                for i in range(len(self.gitr['plasmaSpecies'])): #self.plasmaSpecies.iteritems():
                    prj=self.gitr['plasmaSpecies'][i]
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
            for i in range(len(self.gitr['plasmaSpecies'])): #self.plasmaSpecies.iteritems():
                prj=self.gitr['plasmaSpecies'][i]
                maxDepth=[]
                if (self.gitr['fluxFraction'][i] > 0.0) and not (all(k==0 for k in self.weightAngle[i])):
                    print(('running F-Tridyn for {0} with flux fraction = {1}\n'.format(prj, self.gitr['fluxFraction'][i])))
                    print(('\t and not all angle weights are zero; max angleWeight is {}\n'.format(max(self.weightAngle[i]))))

                    #component/method calls now include arguments (variables)

                    self.services.call(ftridyn, 'init', timeStamp, dTime=time, fPrj=prj, fTargetList=targetList, ftParameters=self.ftridyn , fEnergyIn=self.energyIn[i], fAngleIn=self.angleIn[i], fWeightAngle=self.weightAngle[i], ft_folder=self.FT_OUTPUT_FOLDER, input_file=self.ft_input_file[i], otherInFiles=[self.FT_SURFACE_FILE,self.ftx_lay_file[i]], energy_file_name=self.FT_energy_file_name[i], orig_energy_files_path=self.GITR_eadist_output_path[i], orig_energy_files_pattern=self.GITR_eadist_output_file[i], output_file=outFile)

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
                        print("something was implanted.")
                        maxRange=max(maxDepth)
                        self.maxRangeXolotl[i]=maxRange/10.0 #range in nm for Xolotl 
                        print(' ')
                        print(('\t maximum projectile range for {} is {} [A]'.format(prj, maxRange)))
                        ft_output_file=self.ft_output_file[i]
                        #get implantation profile
                        translate_ftridyn_to_xolotl.ftridyn_to_xolotl(ftridynOnePrjOutput=ft_output_prj_file, ftridynFolder=angleFolder, angle=self.angleIn[i], weightAngle=self.weightAngle[i], prjRange=maxRange, nBins=self.xp.parameters['grid'][0], logFile=outFile) #gAngleDistrib=self.angleDistrFile[i]
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
                    yields=get_yields.sputtering_and_reflection(ftridynOneOutOutput=ft_output_file, ftridynFolder=angleFolder, fNImpacts=self.ftridyn['nImpacts'], angle=self.angleIn[i], weightAngle=self.weightAngle[i], logFile=outFile)
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

                #if flux fraction == 0 or all weight angles == 0.0:
                #if (self.gitr['fluxFraction'][i] > 0.0) and not (all(k==0 for k in self.weightAngle[i])):
                else:
                    print(('Skip running FTridyn for {0}'.format(prj))) #, as fraction in plasma is {1}\n'.format(prj, self.gitr['fluxFraction'][i]))
                    if self.gitr['fluxFraction'][i]==0:
                        print('\t as fraction in plasma is {}\n'.format(self.gitr['fluxFraction'][i]))
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

            for i in range(len(self.gitr['plasmaSpecies'])): #self.plasmaSpecies.iteritems():  
                prj=self.gitr['plasmaSpecies'][i]
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
            
            #tridyn.dat always in order: He, W, D, T --> maintain that here, regardless of order in plasmaSpecies
            #for prj in ['He', 'W', 'D', 'T']:
                #i=self.gitr['plasmaSpecies'].index(prj)
            for i in range(len(self.gitr['plasmaSpecies'])):
                prj=self.gitr['plasmaSpecies'][i]                
                ft_output_profile_temp_prj=timeFolder+'/'+self.FT_OUTPUT_PROFILE_TEMP+'_'+prj
                profile=open(ft_output_profile_temp_prj, "r")
                tridynString=profile.read().rstrip('\n')
                combinedTridynString=str(tridynString)+str(self.maxRangeXolotl[i])
                print(('for {0}, fraction in plasma = {1} , and reflection = {2} '.format(prj,self.gitr['fluxFraction'][i], self.rYield[i])))
                print(('\t effective fraction (in plasma * (1-reflection)) = {} '.format(self.gitr['fluxFraction'][i]*(1-self.rYield[i]))))
                combinedFile.write("%s\n" %(str(self.gitr['fluxFraction'][i]*(1-self.rYield[i])))) #self.fluxFraction[i]
                combinedFile.write("%s\n" %(combinedTridynString))
                profile.close()
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

            #calculate effective sputtering yield; i.e., weighted by relative flux of W-to-He
            totalSpYield=0
            for i in range(len(self.gitr['plasmaSpecies'])): #self.plasmaSpecies.iteritems():
                prj=self.gitr['plasmaSpecies'][i]
                print(('contribution of {0} to total sputtering yield = {1} '.format( prj, float(self.gitr['fluxFraction'][i])*float(self.spYield[i]))))
                totalSpYield+=(float(self.gitr['fluxFraction'][i])*float(self.spYield[i])) #self.fluxFraction[i]
            print(('total weighted sputtering yield = {} (passed to Xolotl)\n'.format(totalSpYield)))
            self.xp.parameters['sputtering'] = totalSpYield            
            
            #time and time-step related parameters
            self.xp.parameters['petscArgs']['-ts_final_time']=time+self.driver['LOOP_TIME_STEP']
            print('\n')
            print('XOLOTL: ')
            print(('\t Run from t = {}'.format(time)))
            print(('\t to t = {}'.format(self.xp.parameters['petscArgs']['-ts_final_time'])))
            print(('\t and time-step = {} '.format( self.driver['LOOP_TIME_STEP'])))

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
                print('\t restart mode: modify xolotl parameters that might change at every loop, including networkFile \n')
                self.xp.parameters['networkFile'] = self.XOLOTL_NETWORK_FILE

            #determine if he_conc true/false ; if true, add '-he_conc' to petsc arguments 
            if self.driver['XOLOTL_HE_CONC']=='Last':
                if time+1.5*self.driver['LOOP_TIME_STEP']>self.driver['END_TIME']:  #*1.5, to give marging of error
                    self.petsc_heConc=True
                    print('printing He concentrations in the last loop')
                elif time<(self.driver['END_TIME']-self.driver['LOOP_TIME_STEP']):
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

            while self.xolotlExitStatus=='collapsed':
                
                self.collapsedLoops+=1
                print(' ')

                #set a maximum number of tries
                if self.collapsedLoops<=int(self.driver['MAX_COLLAPSE_LOOPS']):

                    self.services.call(xolotl, 'init', timeStamp, dTime=time, xParameters=self.xp.parameters, output_file=outFile) #, xFtCoupling=self.driver['FTX_COUPLING'])
                    self.services.call(xolotl, 'step', timeStamp, dTime=time, xHe_conc=self.petsc_heConc, xParameters=self.xp.parameters, output_file=outFile, dZipOutput=self.driver['ZIP_XOLOTL_OUTPUT'])

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
                        print('\t ERROR: XOLOTL OVERGRID ')
                        print('END IPS SIMULATION \n')
                        quit()
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
            #shutil.copyfile('last_TRIDYN.dat', 'last_TRIDYN_toBin.dat')
            print('bin Xolotls output')
            binTRIDYN.binTridyn()
            print('...succesfull')
            
            #store xolotls profile output for each loop (not plasma state)          
            #UPDATE: 'toBin' is .h5 format instead of .dat
            #currentXolotlOutputFileToBin='last_TRIDYN_toBin_%f.dat' %time
            #shutil.copyfile('last_TRIDYN_toBin.dat', currentXolotlOutputFileToBin)
            currentXolotlOutputFileToBin='last_TRIDYN_toBin_%f.h5' %time
            shutil.copyfile('last_TRIDYN_toBin.h5', currentXolotlOutputFileToBin)
            currentXolotlOutputFile='last_TRIDYN_%f.dat' %time
            shutil.copyfile('last_TRIDYN.dat', currentXolotlOutputFile)


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

            if self.driver['START_MODE'] != 'NEUTRAL':
                self.driver['START_MODE'] = 'NEUTRAL'
                print(('switched startMode to {} \n'.format(self.driver['START_MODE'])))

            #using while instead of loop to accomodate variable drive time step 
            #--> update time explicitely before (possibly) increasing time step
            time+=self.driver['LOOP_TIME_STEP']
            print(' ')
            print(('after loop {}, check for updates in time steps '.format(self.driver['LOOP_N'])))

            #update Xolotl and driver time steps if needed
            if self.driver['LOOP_TS_FACTOR'] != 1:
                if (self.driver['LOOP_N']%self.driver['LOOP_TS_NLOOPS']==0):
                    print('\t update driver time step and start_stop ')
                    self.driver['LOOP_TIME_STEP']*=self.driver['LOOP_TS_FACTOR']
                    self.xp.parameters['petscArgs']['-start_stop']*=self.driver['LOOP_TS_FACTOR']
                    print(('\t multiplied time step and start_stop by {} '.format(self.driver['LOOP_TS_FACTOR']))) 
                    print(('\t for a new time step = {0} and start_stop ={1} \n'.format(self.driver['LOOP_TIME_STEP'], self.xp.parameters['petscArgs']['-start_stop'])))

                else:
                    print(('\t no update to driver time step ({0}) or start_stop ({1}) \n'.format(self.driver['LOOP_TIME_STEP'] , self.xp.parameters['petscArgs']['-start_stop'])))
            else:
                print(('\t driver time step (({0}) ane start_stop ({1}) unchanged (factor=1) \n'.format( self.driver['LOOP_TIME_STEP'], self.xp.parameters['petscArgs']['-start_stop'])))

            if time+self.driver['LOOP_TIME_STEP']>self.driver['END_TIME']:
                self.driver['LOOP_TIME_STEP']=self.driver['END_TIME']-time
                print(' ')
                print('\t time step longer than needed for last loop ')
                print(('\t adapting driver time step to {} to reach exactly endTime '.format( self.driver['LOOP_TIME_STEP'])))
                self.xp.parameters['petscArgs']['-start_stop']=self.driver['LOOP_TIME_STEP']/10.0
                print(('\t and Xolotls data is saved every (start_stop) = {} '.format( self.xp.parameters['petscArgs']['-start_stop'])))

            print('\n')

            if self.driver['XOLOTL_MAXTS_FACTOR'] != 1:
                if (self.driver['LOOP_N']%self.driver['XOLOTL_MAXTS_NLOOPS']==0):
                    print(('\t change in Xolotls maximum time step after loop {} '.format( self.driver['LOOP_N'])))
                    self.xp.parameters['petscArgs']['-ts_adapt_dt_max']*=self.driver['XOLOTL_MAXTS_FACTOR']
                    print(('\t multiply time step by {0}, for a new time step = {1} \n'.format(self.driver['XOLOTL_MAXTS_FACTOR'] , self.xp.parameters['petscArgs']['-ts_adapt_dt_max'])))
                    if self.xp.parameters['petscArgs']['-ts_adapt_dt_max'] > self.driver['XOLOTL_MAX_TS']:
                        self.xp.parameters['petscArgs']['-ts_adapt_dt_max']=self.driver['XOLOTL_MAX_TS']
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
