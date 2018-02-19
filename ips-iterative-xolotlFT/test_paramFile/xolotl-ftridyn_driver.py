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
        
        print 'using the parameter template file', self.XOLOTL_PARAM_TEMPLATE

        self.services.update_plasma_state()
        self.services.stage_plasma_state()


        #### DRIVER PARAMETERS #####
        
        print '\n'
        print 'this is a test: reading DRIVER parameters from ips config file: '

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
        print 'reading Xolotl default parameters '
        self.xp = xolotl_param_handler.xolotl_params()
        self.xp.read(self.XOLOTL_PARAM_TEMPLATE)

        print 'running Xolotl in ' , self.xp.parameters['dimensions'], 'D'

        #overwrite default Xolotl parameters that are specified in ips.config

        print 'modify Xolotl paramters '
        print 'this is a test: reading XOLOTL parameters from ips config file:'

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

        print 'this is a test: replacing PETSC arguments from ips config file:'
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
        self.wRedepFluxFraction=float(self.GITR_INPUT_PARAMETERS['wRedepFluxFraction'])
        print '\n'
        print 'this is a test: read GITR parameter wRedepFluxFraction, with value ' ,self.wRedepFluxFraction

        #### FTRIDYN PARAMETERS ##### 

        print '\n'
        print 'this is a test: reading FTRIDYN parameters from ips config file: '

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

        if self.ftridyn['inAngleHe'] < 0 :
            self.angleDistrFileHe = self.GITR_OUTPUT_DIR_He +'/'+self.GITR_ANGLE_DISTRIB_FILE
            print '\t angle distribution file for He found; ', self.angleDistrFileHe #test angles are assigned correctly
            self.angleInHe, self.weightAngleHe = numpy.loadtxt(self.angleDistrFileHe, usecols = (0,1) , unpack=True)
        else:
            self.angleInHe=[self.ftridyn['inAngleHe']]
            self.weightAngleHe = [1.0]
            self.angleDistrFileHe=' '
            print '\t He angle value as defined by user' #test angles are assigned correctly
            
        if self.ftridyn['inAngleW'] < 0 :
            self.angleDistrFileW = self.GITR_OUTPUT_DIR_W +'/'+self.GITR_ANGLE_DISTRIB_FILE
            print '\t angle distribution file for W found; ', self.angleDistrFileW #test angles are assigned correctly
            self.angleInW, self.weightAngleW = numpy.loadtxt(self.angleDistrFileW, usecols = (0,1) , unpack=True)
        else:
            self.angleInW=[self.ftridyn['inAngleW']]
            self.weightAngleW = [1.0]
            self.angleDistrFileW=' '
            print '\t W angle value as defined by user' #test angles are assigned correctly

        print '\n'

        #AND MAYBE SOMETHING SIMILAR WITH ENERGIES?

        if self.ftridyn['spYieldHe']<0:
            self.ftridyn['spYieldModeHe']='calculate'
        else:
            self.ftridyn['spYieldModeHe']='fixed'

        if self.ftridyn['spYieldW']<0:
            self.ftridyn['spYieldModeW']='calculate'
        else:
            self.ftridyn['spYieldModeW']='fixed'

        #FTRIDYN FILES
        #prepare input files; i.e., those transferred from FT init (generateInput) to FT step (run code)
        #leave 'others' empty for a pure FT run
        self.other_ft_input_files_He=[self.FT_SURFACE_FILE, self.FTX_LAY_FILE_He]#self.FT_OTHER_INPUT_FILES_He
        self.other_ft_input_files_W=[self.FT_SURFACE_FILE, self.FTX_LAY_FILE_W]#self.FT_OTHER_INPUT_FILES_W

        if self.ftridyn['inEnergyHe'] < 0:
            self.FT_energy_file_name_He = self.FT_ENERGY_INPUT_FILE_He #"He_W0001.ED1"                              
            self.GITR_energy_output_path_He=self.GITR_OUTPUT_DIR_He #where all the energy distribution files are located
            self.GITR_energy_output_file_He=['dist','.dat']
        else:
            self.FT_energy_file_name_He=' '
            self.GITR_energy_output_path_He=' '
            self.GITR_energy_output_file_He=[' ',' ']

        if self.ftridyn['inEnergyW'] < 0:
            self.FT_energy_file_name_W = self.FT_ENERGY_INPUT_FILE_W #"W_W0001.ED1"
            self.GITR_energy_output_path_W=self.GITR_OUTPUT_DIR_W #where all the energy distribution files are located
            self.GITR_energy_output_file_W=['dist','.dat']
        else:
            self.FT_energy_file_name_W=' '
            self.GITR_energy_output_path_W=' '
            self.GITR_energy_output_file_W=[' ',' ']

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
                #Xolotl only outputs He_W0001.LAY; but it's same substrate composition for running W->W
                shutil.copyfile(self.FTX_LAY_FILE_He,self.FTX_LAY_FILE_W)
                
                if (self.ftridyn['totalDepth']==0.0):
                    print '\t Totaldepth from last_TRIDYN.dat'
                    self.ftridyn['nTT']=10*numpy.max(numpy.loadtxt('last_TRIDYN.dat')[:,0])
                else:
                    print '\t totalDepth fixed to ', self.ftridyn['totalDepth']
                    self.ftridyn['nTT']=self.ftridyn['totalDepth']


            self.services.update_plasma_state()

            # B) RUN FTRIDYN FOR He-> W

            #component/method calls now include arguments (variables)
            #keep EnergyIn AngleIn WeightAngle explicitely as they are component (He, W) specific
            #eventually think about dictionary with entries for He and W --> loop over, and have it all written once (not doubled)
            self.services.call(ftridyn, 'init', timeStamp, dTime=time, fPrj='He', fTg='W', ftParameters=self.ftridyn , fEnergyIn=self.ftridyn['inEnergyHe'], fAngleIn=self.angleInHe, fWeightAngle=self.weightAngleHe, ft_folder=self.FT_OUTPUT_FOLDER, input_file=self.FT_INPUT_FILE_He, otherInFiles=self.other_ft_input_files_He, ft_energy_file_name=self.FT_energy_file_name_He, orig_energy_files_path=self.GITR_energy_output_path_He, orig_energy_files_pattern=self.GITR_energy_output_file_He)

            self.services.call(ftridyn, 'step', timeStamp, fPrj='He', fTg='W', fEnergyIn=self.ftridyn['inEnergyHe'], fAngleIn=self.angleInHe, fWeightAngle=self.weightAngleHe,)

            self.services.stage_plasma_state()

            # C) POSTPROCESSING OF He -> W

            # 1) uncompress zipped file (do not remane yet)
            #FTRIDYN folder is compressed to be added to plasma state; regardless of 'zipOutput' parameter's value

            zippedFile=self.FT_OUTPUT_FOLDER+'.zip'
            unzip_output='unzipOutputHe.txt'
            unzipString='unzip %s -d %s >> %s' %(zippedFile,cwd,unzip_output)
            subprocess.call([unzipString], shell=True)


            #2) #get maximum projectile range to ensure bins are added correctly in 'translate_ftridyn_to_xolotl'

            angleFolder=self.FT_OUTPUT_FOLDER+'/ANGLE'
            maxDepth=[]
            for j in range(len(self.angleInHe)):
                if (self.weightAngleHe[j] > 0.0):
                    filePrjHe=angleFolder+str(self.angleInHe[j])+'/'+self.FT_OUTPUT_PRJ_FILE_He
                    depth, bla=numpy.loadtxt(filePrjHe, usecols = (2,3) , unpack=True)
                    maxDepth.append(max(depth))
                    
            maxRange=max(maxDepth)
            #print '\t maximum projectile range for He is ', maxRange, ' [A]'

            #3) get the sputtering yield (or use fixed value)                   

            #script always needed to reformat output for xolotl
            spYieldCalc=translate_ftridyn_to_xolotl.ftridyn_to_xolotl(ftridynOnePrjOutput=self.FT_OUTPUT_PRJ_FILE_He, ftridynOneOutOutput=self.FT_OUTPUT_FILE_He, ftridynFolder=angleFolder, fNImpacts=self.ftridyn['nImpacts'], gAngleDistrib=self.angleDistrFileHe, angle=self.angleInHe, prjRange=max(maxDepth), nBins=self.xp.parameters['grid'][0])

            #overwrite spY value if mode is 'calculate'
            if self.ftridyn['spYieldModeHe']=='calculate':
                self.ftridyn['spYieldHe']=spYieldCalc
            
            #4) save tridyn.dat

            #append output to allTridynNN.dat for each species (and save to what folder?)        
            tempfile = open(self.FT_OUTPUT_PROFILE_TEMP,"r")
            f = open(self.FT_OUTPUT_PROFILE_FINAL_He, "a")
            maxRangeXolotlHe=maxRange/10.0 #range in nm for Xolotl
            f.write('%s %s \n' %(tempfile.read().rstrip('\n'),maxRangeXolotlHe))
            f.close()
            tempfile.close()
            
            #keep copies of tridyn.dat
            shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP,self.FT_OUTPUT_PROFILE_TEMP_He)
            shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, timeFolder+'/'+self.FT_OUTPUT_PROFILE_TEMP_He)
            

            #5) COPY FOLDERS WITH TIME-STAMP & RENAME FOR (Tg,Prj) SPECIES  

            shutil.move(self.FT_OUTPUT_FOLDER,timeFolder+'/'+self.FT_OUTPUT_FOLDER+'_HeW')

            self.services.update_plasma_state()


            # D) RUN FTRIDYN FOR  W -> W
            
            #component/method calls now include arguments (variables)              
            if self.wRedepFluxFraction>0.0:

                self.services.call(ftridyn, 'init', timeStamp, dTime=time, fPrj='He', fTg='W', ftParameters=self.ftridyn , fEnergyIn=self.ftridyn['inEnergyW'], fAngleIn=self.angleInW, fWeightAngle=self.weightAngleW, ft_folder=self.FT_OUTPUT_FOLDER, input_file=self.FT_INPUT_FILE_W , otherInFiles=self.other_ft_input_files_W, ft_energy_file_name=self.FT_energy_file_name_W, orig_energy_files_path=self.GITR_energy_output_path_W, orig_energy_files_pattern=self.GITR_energy_output_file_W)

                self.services.call(ftridyn, 'step', timeStamp, fPrj='W', fTg='W', fEnergyIn=self.ftridyn['inEnergyW'], fAngleIn=self.angleInW, fWeightAngle=self.weightAngleW,)

                self.services.stage_plasma_state()


                # E) POSTPROCESSING FOR W -> W                                 

                # 1) uncompress zipped file (do not remane yet)

                zippedFile=self.FT_OUTPUT_FOLDER+'.zip'
                unzip_output='unzipOutputW.txt'
                unzipString='unzip %s -d %s >> %s' %(zippedFile,cwd,unzip_output)
                subprocess.call([unzipString], shell=True)


                #2) #get maximum projectile range to ensure bins are added correctly in 'translate_ftridyn_to_xolotl'

                angleFolder=self.FT_OUTPUT_FOLDER+'/ANGLE'
                maxDepth=[]
                for j in range(len(self.angleInW)):
                    if (self.weightAngleW[j] > 0.0):
                        filePrjW=angleFolder+str(self.angleInW[j])+'/'+self.FT_OUTPUT_PRJ_FILE_W
                        depth, bla=numpy.loadtxt(filePrjW, usecols = (2,3) , unpack=True)
                        maxDepth.append(max(depth))

                maxRange=max(maxDepth)
                print '\t maximum projectile range for W is ', maxRange, ' [A]'


                #3) get the sputtering yield (or use fixed value)

                #script always needed to reformat output for xolotl
                spYieldCalc=translate_ftridyn_to_xolotl.ftridyn_to_xolotl(ftridynOnePrjOutput=self.FT_OUTPUT_PRJ_FILE_W, ftridynOneOutOutput=self.FT_OUTPUT_FILE_W, ftridynFolder=angleFolder, fNImpacts=self.ftridyn['nImpacts'], gAngleDistrib=self.angleDistrFileW, angle=self.angleInW, prjRange=max(maxDepth), nBins=self.xp.parameters['grid'][0])

                #overwrite spY value if mode is 'calculate'                                                            
                if self.ftridyn['spYieldModeW']=='calculate':
                    self.ftridyn['spYieldW']=spYieldCalc


                #4) save tridyn.dat                                                  

                #append output to allTridynNN.dat for each species (and save to what folder?)
                tempfile = open(self.FT_OUTPUT_PROFILE_TEMP,"r")
                f = open(self.FT_OUTPUT_PROFILE_FINAL_W, "a")
                maxRangeXolotlW=maxRange/10.0 #range in nm for Xolotl
                f.write('%s %s \n' %(tempfile.read().rstrip('\n'),maxRangeXolotlW))
                f.close()
                tempfile.close()

                #keep copies of tridyn.dat
                shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, self.FT_OUTPUT_PROFILE_TEMP_W)
                shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, timeFolder+'/'+self.FT_OUTPUT_PROFILE_TEMP_W)
                

                #5) COPY FOLDERS WITH TIME-STAMP & RENAME FOR (Tg,Prj) SPECIES
                shutil.move(self.FT_OUTPUT_FOLDER,timeFolder+'/'+self.FT_OUTPUT_FOLDER+'_WW')

                print '\n'

            #if W fraction == 0.0
            else:
                print 'Skip running FTridyn for W, as fraction of W in plasma is', self.wRedepFluxFraction , '\n'
                self.ftridyn['spYieldW']=0.0
                maxRangeXolotlW=0.0
                outputFTFileW=open(self.FT_OUTPUT_PROFILE_TEMP, "w")
                outputFTFileW.write(" 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0")
                outputFTFileW.close()

                #keep copies of tridyn.dat                                                  
                shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, self.FT_OUTPUT_PROFILE_TEMP_W)
                shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, timeFolder+'/'+self.FT_OUTPUT_PROFILE_TEMP_W)


            ######species independent ############

            #6) write sputtering yields to file so they can be used by Xolotl
            print 'sputtering Yield due to He is ', self.ftridyn['spYieldHe']
            print 'sputtering Yield due to W redeposition is ', self.ftridyn['spYieldW']
            
            #write sp Yields to file (temp) and append to spYield output (final)
            spTempfile = open(self.FTX_SPUT_YIELDS_FILE_TEMP,"w+")
            spTempfile.write("%s %s %s \n" %(time,self.ftridyn['spYieldHe'],self.ftridyn['spYieldW']))
            spFile = open(self.FTX_SPUT_YIELDS_FILE_FINAL, "a+")
            spFile.write("%s %s %s \n" %(time,self.ftridyn['spYieldHe'],self.ftridyn['spYieldW']))

            spFile.close()
            spTempfile.close()
            
            shutil.copyfile(self.FTX_SPUT_YIELDS_FILE_TEMP,timeFolder+'/'+self.FTX_SPUT_YIELDS_FILE_TEMP)
            shutil.copyfile(self.FTX_SPUT_YIELDS_FILE_FINAL,timeFolder+'/'+self.FTX_SPUT_YIELDS_FILE_FINAL) #perhaps unnecessary
            
            #7) write format tridyn.dat to include W redep in Xolotl
            
            if os.path.exists(self.FT_OUTPUT_PROFILE_TEMP):
                os.remove(self.FT_OUTPUT_PROFILE_TEMP)
            combinedFile = open(self.FT_OUTPUT_PROFILE_TEMP, 'a')
            profileHe=open(self.FT_OUTPUT_PROFILE_TEMP_He, "r")
            combinedFile.write('%s%s \n' %(profileHe.read().rstrip('\n'),maxRangeXolotlHe))
            combinedFile.write("%s \n" %(str(self.wRedepFluxFraction)))
            profileW=open(self.FT_OUTPUT_PROFILE_TEMP_W, "r")
            combinedFile.write('%s%s ' %(profileW.read().rstrip('\n'),maxRangeXolotlW))
            
            combinedFile.close()
            profileHe.close()
            profileW.close()
            

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
            totalSpYield=self.ftridyn['spYieldHe']+float(self.wRedepFluxFraction)*self.ftridyn['spYieldW']
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
                print 'this is a test: he_conc printed in this (and every) loop'
                self.petsc_heConc=True
            elif self.driver['XOLOTL_HE_CONC']=='False':
                print 'this is a test: he_conc not printed in this (or any other) loop'
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

                    #dStartMode=self.DRIVER_START_MODE, dMode=self.driverMode, xHe_conc=self.XOLOTL_HE_CONC
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

            #re-bin last_TRIDYN file                                     
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
