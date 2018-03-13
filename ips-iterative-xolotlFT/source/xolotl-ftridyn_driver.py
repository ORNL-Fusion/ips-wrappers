#! /usr/bin/env python

from  component import Component
import sys
import os
import subprocess
import numpy
import shutil
import translate_xolotl_to_ftridyn
import translate_ftridyn_to_xolotl
import binTRIDYN
import param_handler

class xolotlFtridynDriver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: init')

        #stage input files
        print 'staging input files', self.INPUT_FILES
        self.services.stage_input_files(self.INPUT_FILES)
        print '\t ...input files staged succesfully'    

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
        
        print '\n'
        print 'reading DRIVER parameters from ips config file: ' #this is a test:

        self.driver={}
        for k,v in self.DRIVER_INPUT_PARAMETERS.iteritems():
            if param_handler.is_int(v):
                print '\t integer input parameter ', k, ' = ', v
                self.driver[k]=int(v)
            elif param_handler.is_float(v):
                print '\t float input parameter ', k, ' = ' , v
                self.driver[k]=float(v)
            else:
                print '\t other ', type(v), ' input parameter ', k, ' = ' , v
                self.driver[k]=v

        print '\n'
        self.driverMode=self.driver['START_MODE']
        print 'running IPS from t = ', self.driver['INIT_TIME'] , ' to t = ', self.driver['END_TIME'], ' in steps of dt = ', self.driver['LOOP_TIME_STEP']

        #### XOLOTL PARAMETERS ##### 

        print '\n'
        print 'XOLOTL paramters: '

        #get dimension to read the correct input parameter template
        dim=int(self.XOLOTL_INPUT_PARAMETERS['dimensions'])
        #xolotl_param_template=self.XOLOTL_PARAM_TEMPLATE_PATH+'/paramXolotl_'+str(dim)+'D.txt'
        xolotl_param_template=self.INPUT_DIR+'/paramXolotl_'+str(dim)+'D.txt' 
        print '\t reading Xolotl default parameters from', xolotl_param_template
        self.xp = param_handler.xolotl_params()
        self.xp.read(xolotl_param_template)
        print '\t running Xolotl in ' , dim, 'D'
        print ' '

        #overwrite default Xolotl parameters that are specified in ips.config
        print 'modify XOLOTL paramters '
        print 'reading XOLOTL parameters from ips config file:' #this is a test:

        for k,v in self.XOLOTL_INPUT_PARAMETERS.iteritems():
            if param_handler.is_int(v):
                print '\t integer input parameter ', k, ' = ' , self.xp.parameters[k] , ' with ' , v
                self.xp.parameters[k]=int(v)
            elif param_handler.is_float(v):                
                print '\t float input parameter ', k, ' = ' , self.xp.parameters[k] , ' with ' , v
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
                print '\t reading list input parameter ', k, ' = ' , values  #this is a test
                self.xp.parameters[k]=values
            else:
                print '\t reading string input parameter ', k, ' = ' , v  #this is a test
                self.xp.parameters[k]=v

        print ' '
        print 'replacing PETSC arguments from ips config file:' #this is a test:
        for k,v in self.XOLOTL_INPUT_PETSCARGS.iteritems():
            if param_handler.is_int(v):
                print '\t integer ', k, ' = ' , self.xp.parameters['petscArgs'][k] , ' with ' , v
                self.xp.parameters['petscArgs'][k]=int(v)
            elif param_handler.is_float(v):
                print '\t float ', k, ' = ' , self.xp.parameters['petscArgs'][k] , ' with ' , v
                self.xp.parameters['petscArgs'][k]=float(v)
            else:
                print '\t other ', type(v), ' argument ', k, ' = ' , self.xp.parameters['petscArgs'][k] , ' with ' , v
                self.xp.parameters['petscArgs'][k]=v

        #if not coupling, delete -tridyn from petsc arguments to not print TRIDYN_*.dat files
        if self.driver['FTX_COUPLING']=='False':
            del self.xp.parameters['petscArgs']['-tridyn']

        #CONTROL WHICH PROCESSES ARE ON:
        #delete Xolotl processes that are specified as false in ips.config
        print ' '
        xp_processes={}
        for k,v in self.XOLOTL_INPUT_PROCESSES.iteritems(): #specified in ips.config
            xp_processes[k]=v

        processList=self.xp.parameters['process'] #all processes, specified in param template
        processString=''

        for p in range(len(processList)):
            key=processList[p]
            if (key in xp_processes) and (xp_processes[key]=='false'):
                print 'delete ', key, ' from Xolotl params, as it was set to false in ips.config' #this is a test
            else:
                processString+=key+' '
        print 'processes included in Xolotl are: ', processString.strip() #this is a test
        self.xp.parameters['process']=processString.strip()


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
        
        self.gitr = {} #xolotl_param_handler.xolotl_params()
        #self.gitr=param_handler.read(self.INPUT_DIR+'/'+self.GITR_OUTPUT_FILE)
        self.gitr=param_handler.read(self.SUBMIT_DIR+'/'+self.GITR_OUTPUT_FILE)

        print 'read (from file) : ', self.SUBMIT_DIR+'/'+self.GITR_OUTPUT_FILE
        for k,v, in self.gitr.iteritems():
            print k, ' : ' , v
        print ' '

        for k,v in self.GITR_INPUT_PARAMETERS.iteritems(): #self.gitr.parameters.iteritems():
            if param_handler.is_int(v):
                print '\t reading integer GITR input parameter ', k, ' = ' , v #this is a test
                self.gitr[k]=int(v)
            elif param_handler.is_float(v):
                print '\t reading float GITR input parameter ', k, ' = '  , v #this is a test
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
                print '\t reading list GITR input parameter ', k, ' = ' , values  #this is a test
                self.gitr[k]=values                
            else:
                print '\t reading string GITR input parameter ', k, ' = ' , v  #this is a test
                self.gitr[k]=v
        print ' '

        print '\t GITR output will be read from ' , self.gitr['gitrOutputDir']  #this is a test

        ##NOT SURE IF THERE ARE MORE XOLOTL OR FT PARAMETERS THAT NEED TO OVERWRITTEN
        if 'flux' in self.gitr:
            self.xp.parameters['flux']=self.gitr['flux']
            print 'replaced flux in Xolotl by value given by GITR, ', self.gitr['flux']  #this is a test
        else:
            print 'no flux specified in GITR, so using values in Xolotl, ', self.xp.parameters['flux'] #this is a test


        #### FTRIDYN PARAMETERS ##### 
        ##LOOP OVER LIST OF PLASMA SPECIES SPECIES #########
        print '\n'
        print 'reading FTRIDYN parameters from ips config file: ' #this is a test:

        self.ftridyn={}
        for k,v in self.FTRIDYN_INPUT_PARAMETERS.iteritems():
            if param_handler.is_int(v):
                print '\t integer input parameter ', k, ' = ', v
                self.ftridyn[k]=int(v)
            elif param_handler.is_float(v):
                print '\t float input parameter ', k, ' = ' , v
                self.ftridyn[k]=float(v)
            else:
                print '\t other ', type(v), ' input parameter ', k, ' = ' , v
                self.ftridyn[k]=v
        print '\n'

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
        self.yieldMode=[]        
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
            print '\t index ',i, 'species ',  self.gitr['plasmaSpecies'][i] #self.plasmaSpecies[i]
            print '\t energy ' , self.energyIn[i] , ' angle ' , self.inAngle[i] 
            print '\t spYield ' , self.spYield[i] , ' rYield ' , self.rYield[i], ' fluxFraction ', self.gitr['fluxFraction'][i] #self.fluxFraction[i]  
            

            if self.inAngle[i] < 0 :
                ##ADD +'_'+prj in the middle if the full path name for ITER cases, as multiple plasma species will follow distributions
                self.angleFile.append(self.gitr['gitrOutputDir'].strip()+'/'+self.GITR_ANGLE_FILE.strip())
                self.aWeightFile.append(self.gitr['gitrOutputDir'].strip()+'/'+self.GITR_AWEIGHT_FILE.strip())
                print '\t reading angles and weights for ', prj , ' from; ', self.angleFile[i], self.aWeightFile[i]
                a = numpy.loadtxt(self.angleFile[i], usecols = (0) , unpack=True)
                w = numpy.loadtxt(self.aWeightFile[i], usecols = (0) , unpack=True)
                self.angleIn.append(a)
                self.weightAngle.append(w)
            else:
                self.angleIn.append([self.inAngle[i]])
                self.weightAngle.append([1.0])
                self.angleFile.append('')
                self.aWeightFile.append('')
                print '\t ',prj, ' angle value as defined by user' #test angles are assigned correctly

            print '\n'

            #AND MAYBE SOMETHING SIMILAR WITH ENERGIES?

            if self.spYield[i]<0 or self.rYield[i]<0:
                self.yieldMode.append('calculate')                
            else:
                self.yieldMode.append('fixed')

            #FTRIDYN FILES
            #prepare input files; i.e., those transferred from FT init (generateInput) to FT step (run code)
            #leave 'others' empty for a pure FT run

            if self.energyIn[i] < 0:
                ##ADD +'_'+prj in the middle if the full path name for ITER cases, as multiple plasma species will follow distributions
                #print 'get information about energy and angle from: ',self.gitr['gitrOutputDir']#.strip()+'_'+prj
                self.FT_energy_file_name.append(self.ft_energy_input_file[i]) #"He_W0001.ED1"
                self.GITR_eadist_output_path.append(self.gitr['gitrOutputDir'].strip())#+'_'+prj) #where all the energy distribution files are located
                self.GITR_eadist_output_file.append(['dist','.dat'])#(self.GITR_EADIST_FILE)
            else:
                self.FT_energy_file_name.append('')
                self.GITR_eadist_output_path.append('')
                self.GITR_eadist_output_file.append([' ',' '])#('')
            
            #initialize maxRangeXolotl list
            self.maxRangeXolotl.append(0.0)


        #MAYBE THIS CAN ALSO BE WRITTEN MORE ELEGANTLY
        if (self.driver['START_MODE']=='RESTART'):
            restart_files = self.services.get_config_param('RESTART_FILES')
            restart_list = restart_files.split()
            for index in range(len(restart_list)):
                filepath='../../restart_files/'+restart_list[index]
                shutil.copyfile(filepath,restart_list[index])

        self.services.update_plasma_state()

    def step(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: step')

        cwd = self.services.get_working_dir()

        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')
        
        self.services.stage_plasma_state() 

        time=self.driver['INIT_TIME']

        #for time in numpy.arange(self.initTime,self.endTime,self.timeStep):
        while time<self.driver['END_TIME']:

            self.services.stage_plasma_state()
            print 'driver time (in loop)  %f' %(time)
            self.services.update_plasma_state()

            #keep all files to be saved (not plasma state) in folder with time stamp
            timeFolder='t'+str(time)
            if not os.path.exists(timeFolder):
                os.makedirs(timeFolder)
            print '\n output of this time-loop will be saved in ', timeFolder

            self.driver['LOOP_N']+=1
            self.collapsedLoops=0 #reset 
            self.xolotlExitStatus='collapsed'

            print '\n set xolotl exit status back to collapsed'

            ###################################### 
            ############## run FTridyn ############
            #### for each (Tg,Prj) combination ####
            ###################################### 

            print '\n'
            print 'F-TRIDYN:'

            # A) GET INPUT THAT MIGHT CHANGE EVERT LOOP READY

            #determine parameters related to init/restart
            iW=self.gitr['plasmaSpecies'].index('W')
            if (self.driverMode == 'INIT'):
                print('\t init mode yes')
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
                print('\t init mode no')
                self.ftridyn['iQ0']=-1
                self.ftridyn['nDataPts'] = translate_xolotl_to_ftridyn.xolotlToLay(totalDepth=self.ftridyn['totalDepth'])
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
                print '\t passing to F-Tridyn the list of targets ' , targetList


                #Xolotl only outputs He_W0001.LAY; but it's same substrate composition for running all projectiles
                iHe=self.gitr['plasmaSpecies'].index('He')
                for i in range(len(self.gitr['plasmaSpecies'])): #self.plasmaSpecies.iteritems():
                    prj=self.gitr['plasmaSpecies'][i]
                    if prj!='He': #do not self-copy
                        print '\t copy ', self.ftx_lay_file[iHe], ' as ', self.ftx_lay_file[i] #this is a test:
                        shutil.copyfile(self.ftx_lay_file[iHe],self.ftx_lay_file[i])
                
                if (self.ftridyn['totalDepth']==0.0):
                    print '\t Totaldepth from last_TRIDYN.dat'
                    self.ftridyn['nTT']=10*numpy.max(numpy.loadtxt('last_TRIDYN.dat')[:,0])
                else:
                    print '\t totalDepth fixed to ', self.ftridyn['totalDepth']
                    self.ftridyn['nTT']=self.ftridyn['totalDepth']


            self.services.update_plasma_state()

            # B) RUN FTRIDYN

            for i in range(len(self.gitr['plasmaSpecies'])): #self.plasmaSpecies.iteritems():
                prj=self.gitr['plasmaSpecies'][i]

                if self.gitr['fluxFraction'][i] > 0.0: #self.fluxFraction[i]>0.0:                     
                    print 'running F-Tridyn for ', prj  , ' with flux fraction = ', self.gitr['fluxFraction'][i], ' > 0.0' #this is a test

                    #component/method calls now include arguments (variables)
                    self.services.call(ftridyn, 'init', timeStamp, dTime=time, fPrj=prj, fTargetList=targetList, ftParameters=self.ftridyn , fEnergyIn=self.energyIn[i], fAngleIn=self.angleIn[i], fWeightAngle=self.weightAngle[i], ft_folder=self.FT_OUTPUT_FOLDER, input_file=self.ft_input_file[i], otherInFiles=[self.FT_SURFACE_FILE,self.ftx_lay_file[i]], energy_file_name=self.FT_energy_file_name[i], orig_energy_files_path=self.GITR_eadist_output_path[i], orig_energy_files_pattern=self.GITR_eadist_output_file[i])

                    self.services.call(ftridyn, 'step', timeStamp, fEnergyIn=self.energyIn[i], fAngleIn=self.angleIn[i], fWeightAngle=self.weightAngle[i])

                    self.services.stage_plasma_state()


                    # C) POSTPROCESSING OF prj -> W

                    # 1) uncompress zipped file (do not remane yet)
                    #FTRIDYN folder is compressed to be added to plasma state; regardless of 'zipOutput' parameter's value

                    zippedFile=self.FT_OUTPUT_FOLDER+'.zip'
                    unzip_output='unzipOutput'+prj+'.txt'
                    unzipString='unzip %s -d %s >> %s' %(zippedFile,cwd,unzip_output)
                    subprocess.call([unzipString], shell=True)


                    #2) #get maximum projectile range to ensure bins are added correctly in 'translate_ftridyn_to_xolotl'

                    ft_output_prj_file=self.ft_output_prj_file[i]
                    angleFolder=self.FT_OUTPUT_FOLDER+'/ANGLE'

                    maxDepth=[]
                    for j in range(len(self.angleIn[i])):
                        if (self.weightAngle[i][j] > 0.0):
                            filePrj=angleFolder+str(self.angleIn[i][j])+'/'+ft_output_prj_file
                            depth, bla=numpy.loadtxt(filePrj, usecols = (2,3) , unpack=True)
                            maxDepth.append(max(depth))
                    
                    maxRange=max(maxDepth)
                    self.maxRangeXolotl[i]=maxRange/10.0 #range in nm for Xolotl 
                    print ' '
                    print '\t maximum projectile range for', prj , ' is ', maxRange, ' [A]'


                    #3) get the sputtering yield (or use fixed value)                   
                
                    print '\n'
                    #script always needed to reformat output for xolotl
                    #outputs sputtering and reflection yields
                    ft_output_file=self.ft_output_file[i]
                    yields=translate_ftridyn_to_xolotl.ftridyn_to_xolotl(ftridynOnePrjOutput=ft_output_prj_file, ftridynOneOutOutput=ft_output_file, ftridynFolder=angleFolder, fNImpacts=self.ftridyn['nImpacts'], angle=self.angleIn[i], weightAngle=self.weightAngle[i], prjRange=maxRange, nBins=self.xp.parameters['grid'][0]) #gAngleDistrib=self.angleDistrFile[i]
                    
                    #overwrite spY value if mode is 'calculate'
                    if self.yieldMode[i]=='calculate':
                        self.spYield[i]=float(yields[0])
                        self.rYield[i]=float(yields[1])

                    #4) save tridyn.dat

                    #append output to allTridynNN.dat for each species (and save to what folder?)        
                    ft_output_profile_final=self.FT_OUTPUT_PROFILE_FINAL+'_'+prj
                    tempfile = open(self.FT_OUTPUT_PROFILE_TEMP,"r")
                    f = open(ft_output_profile_final, "a")                    
                    f.write('%s %s \n' %(tempfile.read().rstrip('\n'),self.maxRangeXolotl[i]))                    
                    f.close()
                    tempfile.close()
            
                    #keep copies of tridyn.dat
                    ft_output_profile_temp_prj=self.FT_OUTPUT_PROFILE_TEMP+'_'+prj #for each species
                    shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP,ft_output_profile_temp_prj)
                    shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, timeFolder+'/'+ft_output_profile_temp_prj)
            
                
                    #5) COPY FOLDERS WITH TIME-STAMP & RENAME FOR (Tg,Prj) SPECIES  
                
                    shutil.move(self.FT_OUTPUT_FOLDER,timeFolder+'/'+self.FT_OUTPUT_FOLDER+'_'+prj+'W')                
                    self.services.update_plasma_state()

                    print '\n'

                #if flux fraction == 0:
                else:

                    print 'Skip running FTridyn for ' , prj , ' as fraction in plasma is', self.gitr['fluxFraction'][i] , '\n' #self.fluxFraction[i] , '\n'
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

                    
            #end of for loop:
            ######species independent ############

            #6) write sputtering yields to file so they can be used by Xolotl

            yieldString=str(time)
            print 'Sputtering and Reflection Yields due to: '
            for i in range(len(self.gitr['plasmaSpecies'])): #self.plasmaSpecies.iteritems():  
                prj=self.gitr['plasmaSpecies'][i]
                print '\t %s :  spY = %s and rY = %s ' %(prj,str(self.spYield[i]), str(self.rYield[i]))
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
            
            #7) write format tridyn.dat to include W redep in Xolotl
            
            if os.path.exists(self.FT_OUTPUT_PROFILE_TEMP):
                os.remove(self.FT_OUTPUT_PROFILE_TEMP)
            combinedFile = open(self.FT_OUTPUT_PROFILE_TEMP, 'a')
            
            #tridyn.dat always in order: He, W, D, T --> maintain that here, regardless of order in plasmaSpecies
            for prj in ['He', 'W', 'D', 'T']:
                i=self.gitr['plasmaSpecies'].index(prj)
                ft_output_profile_temp_prj=timeFolder+'/'+self.FT_OUTPUT_PROFILE_TEMP+'_'+prj
                profile=open(ft_output_profile_temp_prj, "r")
                tridynString=profile.read().rstrip('\n')
                combinedTridynString=str(tridynString)+str(self.maxRangeXolotl[i])
                print 'for ', prj, ' fraction in plasma = ', self.gitr['fluxFraction'][i] , ' and reflection = ', self.rYield[i]
                print '\t effective fraction (in plasma * (1-reflection)) = ', self.gitr['fluxFraction'][i]*(1-self.rYield[i])
                combinedFile.write("%s\n" %(str(self.gitr['fluxFraction'][i]*(1-self.rYield[i])))) #self.fluxFraction[i]
                combinedFile.write("%s\n" %(combinedTridynString))
                profile.close()
            combinedFile.close()

            #compress output
            if self.driver['ZIP_OUTPUT']=='True':
                print 'zip output: ', timeFolder
                zippedTimeFolder=timeFolder+'.zip'
                zip_output='zipOutputTimeFolder.txt'
                zipString='zip -r %s %s >> %s ' %(zippedTimeFolder, timeFolder, zip_output)
                subprocess.call([zipString], shell=True)

                shutil.rmtree(timeFolder)

            else:
                print 'leaving ', timeFolder , 'uncompressed'

            self.services.update_plasma_state()

            ######################################
            ############## run Xolotl ############ 
            ###################################### 

            #Xolotl paramter modifications that need to be done at every loop

            #calculate effective sputtering yield; i.e., weighted by relative flux of W-to-He
            totalSpYield=0

            for i in range(len(self.gitr['plasmaSpecies'])): #self.plasmaSpecies.iteritems():
                prj=self.gitr['plasmaSpecies'][i]
                print 'contribution of ', prj , ' to total sputtering yield = ', float(self.gitr['fluxFraction'][i])*float(self.spYield[i])
                totalSpYield+=(float(self.gitr['fluxFraction'][i])*float(self.spYield[i])) #self.fluxFraction[i]
            print 'total weighted sputtering yield = ', totalSpYield , ' (passed to Xolotl)'

            self.xp.parameters['sputtering'] = totalSpYield            
            self.xp.parameters['petscArgs']['-ts_final_time']=time+self.driver['LOOP_TIME_STEP']
            
            print '\n'
            print 'XOLOTL:'
            print '\t Run from t = ', time
            print '\t to t = ', self.xp.parameters['petscArgs']['-ts_final_time']
            print '\t and time-step = ', self.driver['LOOP_TIME_STEP']
            print '\n'

            if self.driverMode == 'INIT':
                print '\t init mode: modify xolotl parameters that might change at every loop'
            elif self.driverMode == 'RESTART':
                #add (or replace) networkFile line to parameter file
                print '\t restart mode: modify xolotl parameters that might change at every loop, including networkFile line '
                self.xp.parameters['networkFile'] = self.XOLOTL_NETWORK_FILE

            #determine if he_conc true/false ; if true, add '-he_conc' to petsc arguments 
            if self.driver['XOLOTL_HE_CONC']=='Last':
                if time+1.5*self.driver['LOOP_TIME_STEP']>self.driver['END_TIME']:  #*1.5, to give marging of error                                                    
                    self.petsc_heConc=True
                    print 'printing He concentrations in the last loop'
                elif time<(self.driver['END_TIME']-self.driver['LOOP_TIME_STEP']):
                    self.petsc_heConc=False
            elif self.driver['XOLOTL_HE_CONC']=='True':
                print '\t he_conc printed in this (and every) loop' #this is a test:
                self.petsc_heConc=True
            elif self.driver['XOLOTL_HE_CONC']=='False':
                print '\t he_conc not printed in this (or any other) loop' #this is a test:
                self.petsc_heConc=False

            if self.petsc_heConc:
                self.xp.parameters['petscArgs']['-helium_conc'] = ''

            #-check_collapse option in petcs args:
            #exit status is printed to solverStatus.txt  'good' (successful run); 'collapsed' (ts below threshold) ; 'diverged' otherwise        
            #Xolotl is launched again (same paramter and network files) until run successfully, or up to maxCollapseLoop tries,

            while self.xolotlExitStatus=='collapsed':
                
                self.collapsedLoops+=1

                #set a maximum number of tries
                if self.collapsedLoops<=int(self.driver['MAX_COLLAPSE_LOOPS']):

                    self.services.call(xolotl, 'init', timeStamp, dTime=time, xFtCoupling=self.driver['FTX_COUPLING'], xParameters=self.xp.parameters)
                    self.services.call(xolotl, 'step', timeStamp, dTime=time, dZipOutput=self.driver['ZIP_OUTPUT'], xHe_conc=self.petsc_heConc, xParameters=self.xp.parameters)

                    self.services.stage_plasma_state()

                    #print 'reading status from ', str(self.EXIT_STATUS_XOLOTL)
                    statusFile=open(self.XOLOTL_EXIT_STATUS, "r")
                    self.xolotlExitStatus=statusFile.read().rstrip('\n')
                
                    print '\t Xolotl ended simulation with status', self.xolotlExitStatus

                    if self.xolotlExitStatus=='good':
                        print '\t Xolotl successfully executed after ', self.collapsedLoops, ' tries \n'
                        print 'continue with IPS simulation \n'

                    elif self.xolotlExitStatus=='diverged':
                        print '\t ERROR: XOLOTL SOLVER DIVERGED \n'
                        print 'END IPS SIMULATION \n'
                        quit()

                    elif self.xolotlExitStatus=='collapsed':
                        print '\t WARNING: simulation exited loop with status collapse'
                        print '\t try number ', self.collapsedLoops, ' out of ', self.maxCollapseLoops, '\n'

                    else:
                        print '\t WARNING: Xolotl exit status UNKOWN -- IPS simulation continues \n'
                
                else: #reached maximum number of tries for collapsing time steps
                    print '\t ERROR: reached maximum number of tries for collapsing time steps without a successful run \n'
                    print 'END IPS SIMULATION \n'
                    quit()


            shutil.copyfile('last_TRIDYN.dat', 'last_TRIDYN_toBin.dat')
            binTRIDYN.binTridyn()
            
            #store xolotls profile output for each loop (not plasma state)          
            currentXolotlOutputFileToBin='last_TRIDYN_toBin_%f.dat' %time
            shutil.copyfile('last_TRIDYN_toBin.dat', currentXolotlOutputFileToBin)
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
                print 'switched driverMode to ', self.driverMode

            if self.driver['START_MODE'] != 'NEUTRAL':
                self.driver['START_MODE'] = 'NEUTRAL'
                print 'switched startMode to ',  self.driver['START_MODE']

            #using while instead of loop to accomodate variable drive time step 
            #--> update time explicitely before (possibly) increasing time step
            time+=self.driver['LOOP_TIME_STEP']

            #update Xolotl and driver time steps if needed
            if self.driver['LOOP_TS_FACTOR'] != 1:
                if (self.driver['LOOP_N']%self.driver['LOOP_TS_NLOOPS']==0):
                    print '\n',"change in driver's time step (and start_stop) after loop ", self.driver['LOOP_N']
                    self.driver['LOOP_TIME_STEP']*=self.driver['LOOP_TS_FACTOR']
                    self.xp.parameters['petscArgs']['-start_stop']*=self.driver['LOOP_TS_FACTOR']
                    print 'multiplied time step and start_stop by ', self.driver['LOOP_TS_FACTOR']
                    print '\t for a new time step = ', self.driver['LOOP_TIME_STEP'], ' and start_stop = ', self.xp.parameters['petscArgs']['-start_stop']
                else:
                    print '\n', 'in loop ', self.driver['LOOP_N'] ,' no update to '
                    print '\t driver time step (', self.driver['LOOP_TIME_STEP'] , ') or start_stop (' , self.xp.parameters['petscArgs']['-start_stop'], ')'
            else:
                print '\n', 'in loop ', self.driver['LOOP_N'],' no change (factor=1) to ' 
                print '\t driver time step ', self.driver['LOOP_TIME_STEP'], ' and start_stop ' , self.xp.parameters['petscArgs']['-start_stop'] , ')'


            if time+self.driver['LOOP_TIME_STEP']>self.driver['END_TIME']:
                self.driver['LOOP_TIME_STEP']=self.driver['END_TIME']-time
                print 'time step longer than needed for last loop '
                print 'adapting driver time step to ', self.driver['LOOP_TIME_STEP'] ,' to reach exactly endTime'
                self.xp.parameters['petscArgs']['-start_stop']=self.driver['LOOP_TIME_STEP']/10.0
                print 'and Xolotls data is saved every (start_stop) ', self.xp.parameters['petscArgs']['-start_stop']


            if self.driver['XOLOTL_MAXTS_FACTOR'] != 1:
                if (self.driver['LOOP_N']%self.driver['XOLOTL_MAXTS_NLOOPS']==0):
                    print '\n',"change in Xolotl's maximum time step after loop ", self.driver['LOOP_N']
                    self.xp.parameters['petscArgs']['-ts_adapt_dt_max']*=self.driver['XOLOTL_MAXTS_FACTOR']
                    print 'multiply time step by ', self.driver['XOLOTL_MAXTS_FACTOR'] , ' for a new time step = ', self.xp.parameters['petscArgs']['-ts_adapt_dt_max']
                else:
                    print '\n', 'in loop ', self.driver['LOOP_N'] , ' continue with xolotl dt max ', self.xp.parameters['petscArgs']['-ts_adapt_dt_max']

            else:
                print 'Xolotls max time step unchanged (factor=1)'

            self.services.update_plasma_state()

    def finalize(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: finalize')

        #can we add compressing output here? e.g., last_TRIDYN, xolotlStop...
        #and remove large output files? e.g., FTRIDYN.zip

        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')
        self.services.call(ftridyn, 'finalize', timeStamp)
        self.services.call(xolotl, 'finalize', timeStamp)
