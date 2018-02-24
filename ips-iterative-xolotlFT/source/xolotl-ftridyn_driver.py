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
import xolotl_param_handler

class xolotlFtridynDriver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: init')

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
            if xolotl_param_handler.is_int(v):
                print '\t integer input parameter ', k, ' = ', v
                self.driver[k]=int(v)
            elif xolotl_param_handler.is_float(v):
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
        xolotl_param_template=self.XOLOTL_PARAM_TEMPLATE_PATH+'/paramXolotl_'+str(dim)+'D.txt'
        print '\t reading Xolotl default parameters from', xolotl_param_template
        self.xp = xolotl_param_handler.xolotl_params()
        self.xp.read(xolotl_param_template)
        print '\t running Xolotl in ' , dim, 'D'
        print ' '

        #overwrite default Xolotl parameters that are specified in ips.config
        print 'modify XOLOTL paramters '
        print 'reading XOLOTL parameters from ips config file:' #this is a test:

        for k,v in self.XOLOTL_INPUT_PARAMETERS.iteritems():
            if xolotl_param_handler.is_int(v):
                print '\t integer input parameter ', k, ' = ' , self.xp.parameters[k] , ' with ' , v
                self.xp.parameters[k]=int(v)
            elif xolotl_param_handler.is_float(v):                
                print '\t float input parameter ', k, ' = ' , self.xp.parameters[k] , ' with ' , v
                self.xp.parameters[k]=float(v)
            else:
                print '\t other ', type(v), ' input parameter ', k, ' = ' , self.xp.parameters[k] , ' with ' , v
                self.xp.parameters[k]=v

        print ' '
        print 'replacing PETSC arguments from ips config file:' #this is a test:
        for k,v in self.XOLOTL_INPUT_PETSCARGS.iteritems():
            if xolotl_param_handler.is_int(v):
                print '\t integer ', k, ' = ' , self.xp.parameters['petscArgs'][k] , ' with ' , v
                self.xp.parameters['petscArgs'][k]=int(v)
            elif xolotl_param_handler.is_float(v):
                print '\t float ', k, ' = ' , self.xp.parameters['petscArgs'][k] , ' with ' , v
                self.xp.parameters['petscArgs'][k]=float(v)
            else:
                print '\t other ', type(v), ' argument ', k, ' = ' , self.xp.parameters['petscArgs'][k] , ' with ' , v
                self.xp.parameters['petscArgs'][k]=v

        #if not coupling, delete -tridyn from petsc arguments to not print TRIDYN_*.dat files
        if self.driver['FTX_COUPLING']=='False':
            del self.xp.parameters['petscArgs']['-tridyn']

        ### YET TO FIGURE OUT ###
        #delete Xolotl 'process'es that are specified as False in ips.config
        #xolotl_processes=xp.parameters['process']
        #for i in range(len(xolotl_processes)):            
        #    if not self.xolotl_processes[i]:
        #        xp.parameters['process'].replace(xProcess,'')        
        #        del xp.parameters['process'].


        ### GITR RELATED PARAMETERS ###
        #turn this into for
        #for k,v in self.GITR_INPUT_PARAMETERS.iteritems():
            #if xolotl_param_handler.is_int(v):
            #    print '\t this is a test: replacing integer GITR input parameter ', k, ' = ' , self.xp.parameters[k] , ' with ' , v
            #    self.gitr.parameters[k]=int(v)
            #elif xolotl_param_handler.is_float(v):
            #    print '\t this is a test: replacing float GITR input parameter ', k, ' = ' , self.xp.parameters[k] , ' with ' , v
            #    self.gitr.parameters[k]=float(v)
            #else:
            #    print '\t this is a test: replacing other ', type(v), ' GITR input parameter ', k, ' = ' , self.xp.parameters[k] , ' with ' , v
            #    self.gitr.parameters[k]=v                


        print '\n'
        print 'read from ips.config GITR_INPUT_PARAMETERS' #this is a test:
        print '\t plasmaSpecies =', self.GITR_INPUT_PARAMETERS['inputFluxFraction']
        print '\t plasmaSpecies =', self.GITR_INPUT_PARAMETERS['plasmaSpecies']
        
        self.plasmaSpecies=self.GITR_INPUT_PARAMETERS['plasmaSpecies']

        #### FTRIDYN PARAMETERS ##### 
        ##LOOP OVER LIST OF PLASMA SPECIES SPECIES #########
        print '\n'
        print 'reading FTRIDYN parameters from ips config file: ' #this is a test:

        self.ftridyn={}
        for k,v in self.FTRIDYN_INPUT_PARAMETERS.iteritems():
            if xolotl_param_handler.is_int(v):
                print '\t integer input parameter ', k, ' = ', v
                self.ftridyn[k]=int(v)
            elif xolotl_param_handler.is_float(v):
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
        self.fluxFraction=[]
        #other parameters
        self.spYield=[]
        self.spYieldMode=[]        
        self.maxRangeXolotl=[]
        #files
        self.angleDistrFile=[]
        self.FT_energy_file_name=[]
        self.GITR_energy_output_path=[]
        self.GITR_energy_output_file=[]

        inputEnergy=self.ftridyn['inputEnergy'].split(' ')
        inputAngle=self.ftridyn['inputAngle'].split(' ')
        inputSpYield=self.ftridyn['inputSpYield'].split(' ')
        inputFluxFraction=self.GITR_INPUT_PARAMETERS['inputFluxFraction'].split(' ')        

        i=0
        for prj in self.plasmaSpecies: 

            self.energyIn.append(float(inputEnergy[i]))
            self.inAngle.append(float(inputAngle[i]))
            self.spYield.append(float(inputSpYield[i]))
            self.fluxFraction.append(float(inputFluxFraction[i]))
            print '\t index ',i, 'species ', self.plasmaSpecies[i]
            print '\t energy ' , self.energyIn[i] , ' angle ' , self.inAngle[i] 
            print '\t spYield ' , self.spYield[i] , ' fluxFraction ', self.fluxFraction[i]  

            if self.inAngle[i] < 0 :
                self.angleDistrFile.append(self.GITR_OUTPUT_DIR+'_'+prj +'/'+self.GITR_ANGLE_DISTRIB_FILE)
                print '\t angle distribution file for ', prj , ' found; ', self.angleDistrFile[i] #test angles are assigned correctly
                a , w = numpy.loadtxt(self.angleDistrFile[i], usecols = (0,1) , unpack=True)
                self.angleIn.append(a)
                self.weightAngle.append(w)
            else:
                self.angleIn.append([self.inAngle[i]])
                self.weightAngle.append([1.0])
                self.angleDistrFile.append(' ')
                print '\t ',prj, ' angle value as defined by user' #test angles are assigned correctly

            print '\n'

            #AND MAYBE SOMETHING SIMILAR WITH ENERGIES?

            if self.spYield[i]<0:
                self.spYieldMode.append('calculate')
            else:
                self.spYieldMode.append('fixed')

            #FTRIDYN FILES
            #prepare input files; i.e., those transferred from FT init (generateInput) to FT step (run code)
            #leave 'others' empty for a pure FT run

            if self.energyIn[i] < 0:
                self.FT_energy_file_name.append(self.FT_ENERGY_INPUT_FILE[i]) #"He_W0001.ED1"                              
                self.GITR_energy_output_path.append(self.GITR_OUTPUT_DIR+'_'+prj) #where all the energy distribution files are located
                self.GITR_energy_output_file.append(['dist','.dat'])
            else:
                self.FT_energy_file_name.append(' ')
                self.GITR_energy_output_path.append(' ')
                self.GITR_energy_output_file.append([' ',' '])
            
            #initialize maxRangeXolotl list
            self.maxRangeXolotl.append(0.0)

            i+=1


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
            if (self.driverMode == 'INIT'):
                print('\t init mode yes')
                self.ftridyn['iQ0']=0
                self.ftridyn['nDataPts'] = 100 #same as default value in generate_ftridyn_input
                if (self.ftridyn['totalDepth']==0.0):
                    self.ftridyn['nTT']=self.ftridyn['initialTotalDepth']
                else:
                    self.ftridyn['nTT']=self.ftridyn['totalDepth']

            else:
                print('\t init mode no')
                self.ftridyn['iQ0']=-1
                self.ftridyn['nDataPts'] = translate_xolotl_to_ftridyn.xolotlToLay(totalDepth=self.ftridyn['totalDepth'])

                #Xolotl only outputs He_W0001.LAY; but it's same substrate composition for running all projectiles
                i=0
                for prj in self.plasmaSpecies:
                    if prj!='He': #do not self-copy
                        print '\t copy ', self.FTX_LAY_FILE[0], ' as ', self.FTX_LAY_FILE[i] #this is a test:
                        shutil.copyfile(self.FTX_LAY_FILE[0],self.FTX_LAY_FILE[i])
                    i+=1
                
                if (self.ftridyn['totalDepth']==0.0):
                    print '\t Totaldepth from last_TRIDYN.dat'
                    self.ftridyn['nTT']=10*numpy.max(numpy.loadtxt('last_TRIDYN.dat')[:,0])
                else:
                    print '\t totalDepth fixed to ', self.ftridyn['totalDepth']
                    self.ftridyn['nTT']=self.ftridyn['totalDepth']


            self.services.update_plasma_state()

            # B) RUN FTRIDYN FOR He-> W

            i=0 ; #keep track of which projectile, as energy, angle, etc. as given as list, e.g.: [E(He), E(W), E(D), E(T)]
            for prj in self.plasmaSpecies:
                
                if self.fluxFraction[i]>0.0: 
                    
                    #component/method calls now include arguments (variables)
                    self.services.call(ftridyn, 'init', timeStamp, dTime=time, fPrj=prj, fTg='W', ftParameters=self.ftridyn , fEnergyIn=self.energyIn[i], fAngleIn=self.angleIn[i], fWeightAngle=self.weightAngle[i], ft_folder=self.FT_OUTPUT_FOLDER, input_file=self.FT_INPUT_FILE[i], otherInFiles=[self.FT_SURFACE_FILE,self.FTX_LAY_FILE[i]], energy_file_name=self.FT_energy_file_name[i], orig_energy_files_path=self.GITR_energy_output_path[i], orig_energy_files_pattern=self.GITR_energy_output_file[i])

                    self.services.call(ftridyn, 'step', timeStamp, fPrj=prj, fTg='W', fEnergyIn=self.energyIn[i], fAngleIn=self.angleIn[i], fWeightAngle=self.weightAngle[i])

                    self.services.stage_plasma_state()

                    # C) POSTPROCESSING OF prj -> W

                    # 1) uncompress zipped file (do not remane yet)
                    #FTRIDYN folder is compressed to be added to plasma state; regardless of 'zipOutput' parameter's value

                    zippedFile=self.FT_OUTPUT_FOLDER+'.zip'
                    unzip_output='unzipOutput'+prj+'.txt'
                    unzipString='unzip %s -d %s >> %s' %(zippedFile,cwd,unzip_output)
                    subprocess.call([unzipString], shell=True)


                    #2) #get maximum projectile range to ensure bins are added correctly in 'translate_ftridyn_to_xolotl'

                    ft_output_prj_file=self.FT_OUTPUT_PRJ_FILE[i]
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
                    ft_output_file=self.FT_OUTPUT_FILE[i]
                    spYieldCalc=translate_ftridyn_to_xolotl.ftridyn_to_xolotl(ftridynOnePrjOutput=ft_output_prj_file, ftridynOneOutOutput=ft_output_file, ftridynFolder=angleFolder, fNImpacts=self.ftridyn['nImpacts'], gAngleDistrib=self.angleDistrFile[i], angle=self.angleIn[i], prjRange=maxRange, nBins=self.xp.parameters['grid'][0])
                
                    #overwrite spY value if mode is 'calculate'
                    if self.spYieldMode[i]=='calculate':
                        self.spYield[i]=spYieldCalc
                        

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

                    print 'Skip running FTridyn for ' , prj , ' as fraction in plasma is', self.fluxFraction[i] , '\n'
                    self.spYield[i]=0.0
                    self.maxRangeXolotl[i]=0.0
                    outputFTFile=open(self.FT_OUTPUT_PROFILE_TEMP, "w")
                    outputFTFile.write("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0")
                    outputFTFile.close()

                    #keep copies of tridyn.dat    
                    ft_output_profile_temp_prj=self.FT_OUTPUT_PROFILE_TEMP+'_'+prj #for each species
                    shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP,ft_output_profile_temp_prj)
                    shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, timeFolder+'/'+ft_output_profile_temp_prj)
                    
                i+=1
                    
            #end of for loop:
            ######species independent ############

            #6) write sputtering yields to file so they can be used by Xolotl
            i=0
            spYString=str(time)
            print 'Sputtering Yield due to '
            for prj in self.plasmaSpecies:                
                print '\t ', prj , ' = ', self.spYield[i]
                spYString+=' ' +str(self.spYield[i])
                i+=1
            
            #write sp Yields to file (temp) and append to spYield output (final)
            spTempfile = open(self.FTX_SPUT_YIELDS_FILE_TEMP,"w+")
            spTempfile.write(spYString)
            spFile = open(self.FTX_SPUT_YIELDS_FILE_FINAL, "a+")
            spFile.write(spYString)
            spFile.write('\n')
            spFile.close()
            spTempfile.close()
            
            shutil.copyfile(self.FTX_SPUT_YIELDS_FILE_TEMP,timeFolder+'/'+self.FTX_SPUT_YIELDS_FILE_TEMP)
            shutil.copyfile(self.FTX_SPUT_YIELDS_FILE_FINAL,timeFolder+'/'+self.FTX_SPUT_YIELDS_FILE_FINAL) #perhaps unnecessary
            
            #7) write format tridyn.dat to include W redep in Xolotl
            
            if os.path.exists(self.FT_OUTPUT_PROFILE_TEMP):
                os.remove(self.FT_OUTPUT_PROFILE_TEMP)
            combinedFile = open(self.FT_OUTPUT_PROFILE_TEMP, 'a')
            
            i=0
            for prj in self.plasmaSpecies:
                ft_output_profile_temp_prj=timeFolder+'/'+self.FT_OUTPUT_PROFILE_TEMP+'_'+prj
                profile=open(ft_output_profile_temp_prj, "r")
                tridynString=profile.read().rstrip('\n')
                combinedTridynString=str(tridynString)+str(self.maxRangeXolotl[i])
                if prj=='He': #no fraction needed for He 
                    combinedFile.write("%s" %(combinedTridynString))
                else: #prj!='He' #include plasma fraction
                    combinedFile.write("\n%s" %(str(self.fluxFraction[i])))
                    combinedFile.write("\n%s" %(combinedTridynString))
                profile.close()
                i+=1            
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
            i=0
            totalSpYield=0
            for prj in self.plasmaSpecies:
                totalSpYield+=(float(self.fluxFraction[i])*self.spYield[i])
                i+=1

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

                    #print 'this is a test: passing xolotl dictionary \n', self.xp.parameters

                    self.services.call(xolotl, 'init', timeStamp, dTime=time, xFtCoupling=self.driver['FTX_COUPLING'], xParameters=self.xp.parameters)
                    self.services.call(xolotl, 'step', timeStamp, dTime=time, dZipOutput=self.driver['ZIP_OUTPUT'], xHe_conc=self.petsc_heConc)

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
                    print 'multiplied time step by ', self.driver['LOOP_TS_FACTOR'], ' for a new time step = ', self.driver['LOOP_TIME_STEP']
                else:
                    print '\n', 'in loop ', self.driver['LOOP_N'] ,' no update to the driver time step ', self.driver['LOOP_TIME_STEP']
            else:
                print 'loop time step unchanged (factor=1)', self.driver['LOOP_TIME_STEP']


            if time+self.driver['LOOP_TIME_STEP']>self.driver['END_TIME']:
                self.driver['LOOP_TIME_STEP']=self.driver['END_TIME']-time
                print 'time step longer than needed for last loop '
                print 'adapting driver time step to ', self.driver['LOOP_TIME_STEP'] ,' to reach exactly endTime'

            self.xp.parameters['petscArgs']['-start_stop']=self.driver['LOOP_TIME_STEP']/10.0
            print 'and Xolotls data is saved every (start_stop) ', self.xp.parameters['petscArgs']['-start_stop']


            if self.driver['XOLOTL_MAXTS_FACTOR'] != 1:
                if (self.driver['LOOP_N']%self.driver['XOLOTL_MAXTS_NLOOPS']==0):
                    print '\n',"change in Xolotl's maximum time step after loop ", self.driver['LOOP_N']
                    print 'type of petscArg -ts_adapt_dt_max is ', type(self.xp.parameters['petscArgs']['-ts_adapt_dt_max'])
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
