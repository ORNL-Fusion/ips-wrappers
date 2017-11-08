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

        #copy the Xolotl paramter template file
        #xolotlTemplateFile='paramXolotlTemplate.txt'
        #print 'copy xolotl template file from ',self.XOLOTL_PARAM_TEMPLATE, ' to ', xolotlTemplateFile 
        #shutil.copyfile(self.XOLOTL_PARAM_TEMPLATE,xolotlTemplateFile)
        
        print 'using the parameter template file', self.XOLOTL_PARAM_TEMPLATE

        self.services.update_plasma_state()
        self.services.stage_plasma_state()


        #### DRIVER PARAMETERS #####

        #start first loop from the beginning (INIT) or from a previous run (RESTART)
        #RESTART mode requires providing a list of input files:
        #for FTridyn: last_TRIDYN.dat; for Xolotl: params.txt (of the last run), networkfile (networkRestart.h5)
        #place them in the 'restart_files' folder. The mode is changed to NEUTRAL after the 1st loop

        self.startMode = 'INIT' # 'RESTART' or 'INIT'
        self.driverMode=self.startMode

        #True=compress output: time folder for FTridyn's output, TRIDYN and heConc for Xolotl's output
        self.zipOutput=True

        self.initTime=0.0
        self.endTime=0.0002
        self.timeStep=0.0001

        print 'running IPS from t = %f to t=%f, in steps of dt=%f' % (self.initTime, self.endTime, self.timeStep)

        #### XOLOTL PARAMETERS ##### 

        #write every parameter that will be used as argumennts in write_xolotl_paramfile function(s) 
        #i.e., anything different from default values (those set to reproduce email-coupling of FTridyn-Xolotl)

        self.xDimensions=1
        if self.xDimensions==1:
            self.fieldsplit_1_pc_type='redundant'
        elif self.xDimensions==2:
            self.fieldsplit_1_pc_type='gamg -fieldsplit_1_ksp_type gmres -ksp_type fgmres -fieldsplit_1_pc_gamg_threshold -1.0'

        print 'running Xolotl in ' , self.xDimensions, 'D'

        self.xolotlCoupling=True
        self.xolotlStartStop=True
        self.xolotlPhaseCut='true'
        self.xolotlMaxVSize=250
        self.xolotlFlux=4.0e4 #ion/nm2
        self.initialV=0.0 #V/nm3 ; e.g., 3.15e-4 V/nm3 = 5ppm
        self.nxGrid=200
        self.dxGrid=0.5
        self.nyGrid=' ' #string, so that it can be replaced by an empty space for 1D paramter file
        self.dyGrid=' '
        self.voidPortion=40 #[%]; std=40%
        if self.startMode=='INIT':
            self.xolotlNetworkFile='notInUse'
        elif self.startMode=='RESTART':
            self.xolotlNetworkFile='networkRestart.h5'
        #True: print at every loop ; False: don't print ; 'Last' (string): print only during the last loop
        self.heConc='Last' 
        self.process='reaction advec modifiedTM diff movingSurface attenuation'

        #bursting: then grouping should happen, phaseCut set to 'false' and smaller maxVSize, e.g., ~50
        self.xolotlBursting=True
        self.xolotlGrouping=False
        self.xolotlGroupHeV=31
        self.xolotlgroupHe=4
        self.xolotlgroupV=4
        
        if self.xolotlBursting:
            self.xolotlGrouping=True
            self.xolotlPhaseCut='false'
            self.xolotlMaxVSize=50

        #CHANGE TO GET FROM FILE 
        self.gFluxFractionW=0.00034 #relative flux of W/He [from GITR!] 
        


        #### FTRIDYN PARAMETERS ##### 

        #TotalDepth: total substrate depth in [A]; set to 0.0 to use what Xolotl passes to ftridyn (as deep as He exists)
        #InitialTotalDepth: if TotalDepth=0.0, choose an appropriate depth for the irradiation energy in the 1st loop
        #     use TotalDepth=0.0 if startMode is RESTART (not understood why, but a fixed totalDepth doesn't work on the 1st loop)
        #NImpacts: number of impacts (NH in generateInput) ;  InEnergy: impact energy (energy in generateInput, [eV]); initialize SpYield
        #if spYield < 0 -> use calculated value; else, use fixed value, usually [0,1) 
            
        self.ftridynTotalDepth=0.0
        self.ftridynInitialTotalDepth=300.0
        self.ftridynNImpacts=1.0e3

        #E or A < 0 -> use distribution(s)
        self.ftridynInEnergyHe=250.0
        self.ftridynInAngleHe=0.0 #wrt surface normal
        self.ftridynInEnergyW=-1
        self.ftridynInAngleW=-1  #wrt surface normal  

        #just have one spYeld to control mode and initialize others to zero
        self.ftridynSpYieldW=-1.0
        self.ftridynSpYieldHe=-1.0

        if self.ftridynInAngleHe < 0 :
            self.angleDistrFileHe = self.GITR_OUTPUT_DIR_He +'/'+self.ANGLE_DISTRIB_FILE
            print '\t angle distribution file for He found; ', self.angleDistrFileHe #test angles are assigned correctly  
            self.angleInHe, self.weightAngleHe = numpy.loadtxt(self.angleDistrFileHe, usecols = (0,1) , unpack=True)
        else:
            self.angleInHe=[self.ftridynInAngleHe] 
            self.weightAngleHe = [1.0]
            self.angleDistrFileHe=' '
            print '\t He angle value as defined by user' #test angles are assigned correctly  

        if self.ftridynInAngleW < 0 :
            self.angleDistrFileW = self.GITR_OUTPUT_DIR_W +'/'+self.ANGLE_DISTRIB_FILE
            print '\t angle distribution file for W found; ', self.angleDistrFileW #test angles are assigned correctly
            self.angleInW, self.weightAngleW = numpy.loadtxt(self.angleDistrFileW, usecols = (0,1) , unpack=True)
        else:
            self.angleInW=[self.ftridynInAngleW]
            self.weightAngleW = [1.0]
            self.angleDistrFileW=' '
            print '\t W angle value as defined by user' #test angles are assigned correctly

        #AND MAYBE SOMETHING SIMILAR WITH ENERGIES?

        if self.ftridynSpYieldHe<0:
            self.ftridynSpYieldModeHe='calculate'
        else:
            self.ftridynSpYieldModeHe='fixed'

        if self.ftridynSpYieldW<0:
            self.ftridynSpYieldModeW='calculate'
        else:
            self.ftridynSpYieldModeW='fixed'

        #FTRIDYN FILES
        #prepare input files; i.e., those transferred from FT init (generateInput) to FT step (run code)
        #leave 'others' empty for a pure FT run
        self.other_ft_input_files_He=[self.SURFACE_FILE, self.LAY_FILE_He]#self.FT_OTHER_INPUT_FILES_He
        self.other_ft_input_files_W=[self.SURFACE_FILE, self.LAY_FILE_W]#self.FT_OTHER_INPUT_FILES_W

        if self.ftridynInEnergyHe < 0:
            self.FT_energy_file_name_He = self.FT_ENERGY_INPUT_FILE_He #"He_W0001.ED1"                              
            self.GITR_energy_output_path_He=self.GITR_OUTPUT_DIR_He #where all the energy distribution files are located
            self.GITR_energy_output_file_He=['dist','.dat']
        else:
            self.FT_energy_file_name_He=' '
            self.GITR_energy_output_path_He=' '
            self.GITR_energy_output_file_He=[' ',' ']


        if self.ftridynInEnergyW < 0:
            self.FT_energy_file_name_W = self.FT_ENERGY_INPUT_FILE_W #"W_W0001.ED1"
            self.GITR_energy_output_path_W=self.GITR_OUTPUT_DIR_W #where all the energy distribution files are located
            self.GITR_energy_output_file_W=['dist','.dat']
        else:
            self.FT_energy_file_name_W=' '
            self.GITR_energy_output_path_W=' '
            self.GITR_energy_output_file_W=[' ',' ']



        #MAYBE THIS CAN ALSO BE WRITTEN MORE ELEGANTLY
        if (self.startMode=='RESTART'):
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

        for time in numpy.arange(self.initTime,self.endTime,self.timeStep):

            self.services.stage_plasma_state()
            print 'driver time (in loop)  %f' %(time)
            self.services.update_plasma_state()

            #keep all files to be saved (not plasma state) in folder with time stamp
            timeFolder=cwd+'/t'+str(time)
            if not os.path.exists(timeFolder):
                os.makedirs(timeFolder)
            print '\n output of this time-loop will be saved in ', timeFolder

            ###################################### 
            ############## run FTridyn ############
            #### for each (Tg,Prj) combination ####
            ###################################### 



            # A) GET INPUT THAT MIGHT CHANGE EVERT LOOP READY

            #determine parameters related to init/restart
            if (self.driverMode == 'INIT'):
                print('\n init mode yes')
                iQ0=0
                nDataPts = 100 #same as default value in generate_ftridyn_input
                if (self.ftridynTotalDepth==0.0):
                    nTT=self.ftridynInitialTotalDepth
                else:
                    nTT=self.ftridynTotalDepth

            else:
                print('\n init mode no')
                iQ0=-1
                nDataPts = translate_xolotl_to_ftridyn.xolotlToLay(totalDepth=self.ftridynTotalDepth)
                #Xolotl only outputs He_W0001.LAY; but it's same substrate composition for running W->W
                shutil.copyfile(self.LAY_FILE_He,self.LAY_FILE_W)
                
                if (self.ftridynTotalDepth==0.0):
                    print '\t Totaldepth from last_TRIDYN.dat'
                    nTT=10*numpy.max(numpy.loadtxt('last_TRIDYN.dat')[:,0])
                else:
                    print '\t totalDepth fixed to ', self.ftridynTotalDepth
                    nTT=self.ftridynTotalDepth


            self.services.update_plasma_state()

            # B) RUN FTRIDYN FOR He-> W

            #component/method calls now include arguments (variables)
            self.services.call(ftridyn, 'init', timeStamp, dTime=time, fPrj='He', fTg='W', fNTT=nTT, fNDataPts=nDataPts, fIQ0=iQ0, fNImpacts=self.ftridynNImpacts, fEnergyIn=self.ftridynInEnergyHe, fAngleIn=self.angleInHe, fWeightAngle=self.weightAngleHe, ft_folder=self.FT_OUTPUT_FOLDER, input_file=self.FT_INPUT_FILE_He, otherInFiles=self.other_ft_input_files_He, ft_energy_file_name=self.FT_energy_file_name_He, orig_energy_files_path=self.GITR_energy_output_path_He, orig_energy_files_pattern=self.GITR_energy_output_file_He)

            self.services.call(ftridyn, 'step', timeStamp, fPrj='He', fTg='W', fEnergyIn=self.ftridynInEnergyHe, fAngleIn=self.angleInHe, fWeightAngle=self.weightAngleHe)

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

            print '\t maximum projectile range for He is ', max(maxDepth), ' [A]'

            #3) get the sputtering yield (or use fixed value)                   

            #script always needed to reformat output for xolotl
            spYieldCalc=translate_ftridyn_to_xolotl.ftridyn_to_xolotl(ftridynOnePrjOutput=self.FT_OUTPUT_PRJ_FILE_He, ftridynOneOutOutput=self.FT_OUTPUT_FILE_He, ftridynFolder=angleFolder, fNImpacts=self.ftridynNImpacts, gAngleDistrib=self.angleDistrFileHe, angle=self.angleInHe, prjRange=max(maxDepth), nBins=self.nxGrid)

            #overwrite spY value if mode is 'calculate'
            if self.ftridynSpYieldModeHe=='calculate':
                self.ftridynSpYieldHe=spYieldCalc
            
            #4) save tridyn.dat

            #append output to allTridynNN.dat for each species (and save to what folder?)        
            tempfile = open(self.FT_OUTPUT_PROFILE_TEMP,"r")
            f = open(self.FT_OUTPUT_PROFILE_FINAL_He, "a")
            f.write(tempfile.read())
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
            if self.gFluxFractionW>0.0:

                self.services.call(ftridyn, 'init', timeStamp, dTime=time, fPrj='W', fTg='W', fNTT=nTT, fNDataPts=nDataPts, fIQ0=iQ0, fNImpacts=self.ftridynNImpacts, fEnergyIn=self.ftridynInEnergyW, fAngleIn=self.angleInW, fWeightAngle=self.weightAngleW, ft_folder=self.FT_OUTPUT_FOLDER, input_file=self.FT_INPUT_FILE_W , otherInFiles=self.other_ft_input_files_W, ft_energy_file_name=self.FT_energy_file_name_W, orig_energy_files_path=self.GITR_energy_output_path_W, orig_energy_files_pattern=self.GITR_energy_output_file_W)
                
                self.services.call(ftridyn, 'step', timeStamp, fPrj='W', fTg='W', fEnergyIn=self.ftridynInEnergyW, fAngleIn=self.angleInW, fWeightAngle=self.weightAngleW)

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

                print '\t maximum projectile range for W is ', max(maxDepth), ' [A]'


                #3) get the sputtering yield (or use fixed value)

                #script always needed to reformat output for xolotl
                spYieldCalc=translate_ftridyn_to_xolotl.ftridyn_to_xolotl(ftridynOnePrjOutput=self.FT_OUTPUT_PRJ_FILE_W, ftridynOneOutOutput=self.FT_OUTPUT_FILE_W, ftridynFolder=angleFolder, fNImpacts=self.ftridynNImpacts, gAngleDistrib=self.angleDistrFileW, angle=self.angleInW, prjRange=max(maxDepth), nBins=self.nxGrid)

                #overwrite spY value if mode is 'calculate'                                                            
                if self.ftridynSpYieldModeW=='calculate':
                    self.ftridynSpYieldW=spYieldCalc


                #4) save tridyn.dat                                                  

                #append output to allTridynNN.dat for each species (and save to what folder?)
                tempfile = open(self.FT_OUTPUT_PROFILE_TEMP,"r")
                f = open(self.FT_OUTPUT_PROFILE_FINAL_W, "a")
                f.write(tempfile.read())
                f.close()
                tempfile.close()

                #keep copies of tridyn.dat
                shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, self.FT_OUTPUT_PROFILE_TEMP_W)
                shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, timeFolder+'/'+self.FT_OUTPUT_PROFILE_TEMP_W)
                

                #5) COPY FOLDERS WITH TIME-STAMP & RENAME FOR (Tg,Prj) SPECIES
                shutil.move(self.FT_OUTPUT_FOLDER,timeFolder+'/'+self.FT_OUTPUT_FOLDER+'_WW')


            #if W fraction == 0.0
            else:
                print 'Skip running FTridyn for W, as fraction of W in plasma is', self.gFluxFractionW
                self.ftridynSpYieldW=0.0
                outputFTFileW=open(self.FT_OUTPUT_PROFILE_TEMP, "w")
                outputFTFileW.write(" 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0")
                outputFTFileW.close()

                #keep copies of tridyn.dat                                                  
                shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, self.FT_OUTPUT_PROFILE_TEMP_W)
                shutil.copyfile(self.FT_OUTPUT_PROFILE_TEMP, timeFolder+'/'+self.FT_OUTPUT_PROFILE_TEMP_W)


            ######species independent ############

            #6) write sputtering yields to file so they can be used by Xolotl
            print 'sputtering Yield due to He is ', self.ftridynSpYieldHe
            print 'sputtering Yield due to W redeposition is ', self.ftridynSpYieldW
            
            #write sp Yields to file (temp) and append to spYield output (final)
            spTempfile = open(self.SPUT_YIELDS_FILE_TEMP,"w+")
            spTempfile.write("%s %s %s \n" %(time,self.ftridynSpYieldHe,self.ftridynSpYieldW))
            spFile = open(self.SPUT_YIELDS_FILE_FINAL, "a+")
            spFile.write("%s %s %s \n" %(time,self.ftridynSpYieldHe,self.ftridynSpYieldW))

            spFile.close()
            spTempfile.close()
            
            shutil.copyfile(self.SPUT_YIELDS_FILE_TEMP,timeFolder+'/'+self.SPUT_YIELDS_FILE_TEMP)
            shutil.copyfile(self.SPUT_YIELDS_FILE_FINAL,timeFolder+'/'+self.SPUT_YIELDS_FILE_FINAL) #perhaps unnecessary
            
            #7) write format tridyn.dat to include W redep in Xolotl
            
            if os.path.exists(self.FT_OUTPUT_PROFILE_TEMP):
                os.remove(self.FT_OUTPUT_PROFILE_TEMP)
            combinedFile = open(self.FT_OUTPUT_PROFILE_TEMP, 'a')
            profileHe=open(self.FT_OUTPUT_PROFILE_TEMP_He, "r")
            combinedFile.write(profileHe.read())
            combinedFile.write("%s \n" %(str(self.gFluxFractionW)))
            profileW=open(self.FT_OUTPUT_PROFILE_TEMP_W, "r")
            combinedFile.write(profileW.read())
            
            combinedFile.close()
            profileHe.close()
            profileW.close()
            

            #compress output
            if self.zipOutput:
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


            if self.heConc=='Last':
                if time+1.5*self.timeStep>self.endTime:  #*1.5, to give marging of error 
                    self.petsc_heConc=True
                    print 'printing He concentrations in the last loop'
                elif time<(self.endTime-self.timeStep):
                    self.petsc_heConc=False


            #calculate effective sputtering yield; i.e., weighted by relative flux of W-to-He
            totalSpYield=self.ftridynSpYieldHe+self.gFluxFractionW*self.ftridynSpYieldW

            self.services.call(xolotl, 'init', timeStamp, dStartMode=self.startMode, dMode=self.driverMode, dTime=time, dTimeStep=self.timeStep, xFtCoupling=self.xolotlCoupling, dZipOutput=self.zipOutput, xParamTemplate=self.XOLOTL_PARAM_TEMPLATE, xNetworkFile=self.xolotlNetworkFile, xDimensions=self.xDimensions, xFieldsplit_1_pc_type=self.fieldsplit_1_pc_type, xStartStop=self.xolotlStartStop, xPhaseCut=self.xolotlPhaseCut, xMaxVSize=self.xolotlMaxVSize,  xFlux=self.xolotlFlux, xInitialV=self.initialV, xNxGrid=self.nxGrid, xNyGrid=self.nyGrid, xDxGrid=self.dxGrid, xDyGrid=self.dyGrid, xBursting=self.xolotlBursting, xGrouping=self.xolotlGrouping, xGroupHeV=self.xolotlGroupHeV, xGroupHe=self.xolotlgroupHe, xGroupV=self.xolotlgroupV, fNImpacts=self.ftridynNImpacts, gFractionW=self.gFluxFractionW, xHe_conc=self.petsc_heConc, xProcess=self.process, xVoidPortion=self.voidPortion, weightedSpYield=totalSpYield)

            self.services.call(xolotl, 'step', timeStamp, dTime=time)

            self.services.stage_plasma_state()


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
            tempfileRet = open(self.RETENTION_XOLOTL_TEMP,"r")
            fRet = open(self.RETENTION_XOLOTL_FINAL, "a")
            fRet.write(tempfileRet.read())
            fRet.close()
            tempfileRet.close()
            
            #surface
            tempfileSurf = open(self.SURFACE_XOLOTL_TEMP,"r")
            fSurf = open(self.SURFACE_XOLOTL_FINAL, "a")
            fSurf.write(tempfileSurf.read())
            fSurf.close()
            tempfileSurf.close()

            #save network file with a different name to use in the next time step
            currentXolotlNetworkFile='xolotlStop_%f.h5' %time
            shutil.copyfile('xolotlStop.h5',currentXolotlNetworkFile)

            #update driver mode after the 1st loop, from INIT to RESTART
            if self.driverMode == 'INIT':
                self.driverMode = 'RESTART'
                print 'switched driverMode to ', self.driverMode

            if self.startMode != 'NEUTRAL':
                self.startMode = 'NEUTRAL'
                print 'switched startMode to ',  self.startMode 

            self.xolotlNetworkFile='xolotlStop_%f.h5' %(time)

            self.services.update_plasma_state()

    def finalize(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: finalize')

        #can we add compressing output here? e.g., last_TRIDYN, xolotlStop...
        #and remove large output files? e.g., FTRIDYN.zip

        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')
        self.services.call(ftridyn, 'finalize', timeStamp)
        self.services.call(xolotl, 'finalize', timeStamp)
