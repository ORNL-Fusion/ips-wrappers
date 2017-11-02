#! /usr/bin/env python

from  component import Component
import os
import shutil
import glob
import sys
import translate_xolotl_to_ftridyn
import translate_ftridyn_to_xolotl
import generate_ftridyn_input #new script to generate FTridyn input
import numpy as np
import subprocess
import re

class ftridynWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        #get name of FTridyn executable from config file
        #This is actually a shell script that calls FTridyn and pipe the input to the executable
        self.ftridyn_exe = self.FTRIDYN_EXE

    def init(self, timeStamp=0.0,**keywords):
        print('\n fridyn_worker: init')
        #stage plasma state files for use on execution of FTridyn
        self.services.stage_plasma_state()

        print 'check that all arguments are read well by ftridyn-init'
        for (k, v) in keywords.iteritems():
            print '\t', k, " = ", v

        #asign a local variable to arguments used multiple times
        initialTotalDepth=keywords['fInitialTotalDepth']
        totalDepth=keywords['fTotalDepth']
        self.driverTime=keywords['dTime']
        self.angleInHe=keywords['fAngleInHe']
        self.angleInW = keywords['fAngleInW']
        self.weightAngleHe=keywords['fWeightAngleHe']
        self.weightAngleW=keywords['fWeightAngleW']
        self.energyInHe=keywords['fEnergyInHe']
        self.energyInW = keywords['fEnergyInW']
        self.gitrOuputFolder_He=keywords['gOutputFolderHe']
        self.gitrOuputFolder_W=keywords['gOutputFolderW']
        self.nImpacts=keywords['fNImpacts']
        self.fluxFractionW=keywords['gFractionW']

        cwd = self.services.get_working_dir()

        #prepare and save (for each loop) ftridyn input
        if (keywords['dMode'] == 'INIT'):
            print('\n init mode yes')
            if (totalDepth==0.0):
                nTT=initialTotalDepth
            else:
                nTT=totalDepth
            
            #generateInput for He -> W
            print '\t calling script to generate FTridyn input with depth=%d, NH=%d and Ein=%f' % (nTT,self.nImpacts,self.energyInHe)
            for j in range(len(self.angleInHe)):
                if (self.weightAngleHe[j] == 0.0):
                    print '\t weight of angle ',  self.angleInHe[j], ' is ' , self.weightAngleHe[j] , ', so skip generateInput'
                elif (self.weightAngleHe[j] > 0.0):
                    generate_ftridyn_input.Prj_Tg_xolotl(depth=nTT,number_histories=self.nImpacts,incident_energy=self.energyInHe,incident_angle=self.angleInHe[j],projectile_name='He',target_name='W')
                    pathFolderHe=cwd+'/'+self.FOLDER_He+str(self.angleInHe[j])# + "_"+str(energyInW[j])
                    if not os.path.exists(pathFolderHe):
                        os.makedirs(pathFolderHe)

                    #copy the input file Prj_Tg0001.IN to its folder
                    shutil.copyfile(cwd+'/'+self.INPUT_FILE_He, pathFolderHe+"/"+self.INPUT_FILE_He)
                    shutil.copyfile(self.SURFACE_FILE, pathFolderHe+"/"+self.SURFACE_FILE)
                    if self.energyInHe < 0:
                        #copy energy distribution file                            
                        shutil.copyfile(self.gitrOuputFolder_He+"/dist"+str(int(j))+".dat", pathFolderHe+"/"+self.ENERGY_INPUT_FILE_He)#"He_W0001.ED1"   

            print '\t DONE with GenerateInput (init) He->W \n'
                    
            
            #run generateInput for W->W: 
            #for W redep, W energy distr given by GITR -> energyInW=-1 -> no loop in energy needed
            if self.fluxFractionW>0.0 :
                for j in range(len(self.angleInW)):
                    if (self.weightAngleW[j] == 0.0):
                        print '\t weight of angle ',  self.angleInW[j], ' is ' , self.weightAngleW[j] , ', so skip generateInput'
                    elif (self.weightAngleW[j] > 0.0):
                        generate_ftridyn_input.Prj_Tg_xolotl(depth=nTT,number_histories=self.nImpacts,incident_energy=self.energyInW,incident_angle=self.angleInW[j],projectile_name='W',target_name='W')
                        pathFolderW = cwd+'/'+self.FOLDER_W+str(self.angleInW[j])# + "_"+str(energyInW[j])
                        if not os.path.exists(pathFolderW):
                            os.makedirs(pathFolderW)
                
                        #copy the input file Prj_Tg0001.IN to its folder
                        shutil.copyfile(cwd+'/'+self.INPUT_FILE_W, pathFolderW+"/"+self.INPUT_FILE_W) 
                        shutil.copyfile(self.SURFACE_FILE, pathFolderW+"/"+self.SURFACE_FILE)
                        if self.energyInW < 0:
                            #copy energy distribution file
                            shutil.copyfile(self.gitrOuputFolder_W+"/dist"+str(int(j))+".dat", pathFolderW+"/"+self.ENERGY_INPUT_FILE_W)#"W_W0001.ED1"

                print '\t DONE with GenerateInput (init) W->W \n'

            else:
                print 'Skip running generateInput for W; as W fraction in plasma is ',  self.fluxFractionW

        else:
            print('\n init mode no')
            nDataPts = translate_xolotl_to_ftridyn.xolotlToLay(totalDepth=totalDepth)
            #newestLay = max(glob.iglob('*.LAY'), key=os.path.getctime)
            #currentFtridynLayFileHe='%s_%f' %(self.layFileHe,self.driverTime)
            #currentFtridynLayFileW='%s_%f' %(self.layFileW,self.driverTime)
            #shutil.copyfile(self.layFileHe, currentFtridynLayFileHe)
            #shutil.copyfile(self.layFileW, currentFtridynLayFileW)

            if (totalDepth==0.0):
                print '\t Totaldepth from last_TRIDYN.dat'
                nTT=10*np.max(np.loadtxt('last_TRIDYN.dat')[:,0]) 
            else:
                print '\t totalDepth fixed to ', totalDepth 
                nTT=totalDepth

            #new script to generate FTridyn input for He->W  
            print '\t calling script to generate FTridyn input with IQO=%d, number_layers=%d, depth=%d, NH=%d and Ein=%f' %(-1,nDataPts, nTT,self.nImpacts,self.energyInHe)
            for j in range(len(self.angleInHe)):
                if (self.weightAngleHe[j] == 0.0):
                    print '\t weight of angle ',  self.angleInHe[j], ' is ' , self.weightAngleHe[j] , ', so skip generateInput'
                elif (self.weightAngleHe[j] > 0.0):
                    generate_ftridyn_input.Prj_Tg_xolotl(IQ0=-1,number_layers=nDataPts,depth=nTT,number_histories=self.nImpacts,incident_energy=self.energyInHe,incident_angle=self.angleInW[j],projectile_name='He',target_name='W')
                    pathFolderHe = cwd+'/'+self.FOLDER_He+str(self.angleInHe[j])# + "_"+str(energyInHe[j])
                    if not os.path.exists(pathFolderHe):
                        os.makedirs(pathFolderHe)


                    #copy the input file Prj_Tg0001.IN to its folder
                    shutil.copyfile(cwd+'/'+self.INPUT_FILE_He, pathFolderHe+"/"+self.INPUT_FILE_He)
                    shutil.copyfile(self.SURFACE_FILE, pathFolderHe+"/"+self.SURFACE_FILE)
                    shutil.copyfile(self.LAY_FILE_He, pathFolderHe+"/"+self.LAY_FILE_He)
                    #shutil.copyfile('surface.surf', pathFolderHe+"/"+'surface.surf')
                    if self.energyInHe < 0:
                        #copy energy distribution file 
                        shutil.copyfile(self.gitrOuputFolder_He+"/dist"+str(int(j))+".dat", pathFolderHe+"/"+self.ENERGY_INPUT_FILE_He)#"He_W0001.ED1"                

            print '\t DONE with GenerateInput (restart) He->W \n'

            #run generateInput for W->W
            #for W redep, W energy distr given by GITR -> energyInW=-1 -> no loop in energy needed
            if self.fluxFractionW>0.0:            
                for j in range(len(self.angleInW)):
                    if (self.weightAngleW[j] == 0.0):
                        print '\t weight of angle ',  self.angleInW[j], ' is ' , self.weightAngleW[j] , ', so skip generateInput'
                    elif (self.weightAngleW[j] > 0.0):
                        generate_ftridyn_input.Prj_Tg_xolotl(IQ0=-1,number_layers=nDataPts,depth=nTT,number_histories=self.nImpacts,incident_energy=self.energyInW,incident_angle=self.angleInW[j],projectile_name='W',target_name='W')
                        pathFolderW = cwd+'/'+self.FOLDER_W+str(self.angleInW[j])# + "_"+str(energyInW[j])
                        if not os.path.exists(pathFolderW):
                            os.makedirs(pathFolderW)
                
                        #copy the input file Prj_Tg0001.IN to its folder
                        shutil.copyfile(cwd+'/'+self.INPUT_FILE_W, pathFolderW+"/"+self.INPUT_FILE_W)
                        #shutil.copyfile('surface.surf', pathFolderW+"/"+'surface.surf')
                        shutil.copyfile(self.SURFACE_FILE, pathFolderW+"/"+self.SURFACE_FILE)
                        shutil.copyfile(self.LAY_FILE_He, pathFolderW+"/"+self.LAY_FILE_W)

                        if self.energyInW < 0:
                            #copy energy distribution file              
                            shutil.copyfile(self.gitrOuputFolder_W+"/dist"+str(int(j))+".dat", pathFolderW+"/"+self.ENERGY_INPUT_FILE_W)#"W_W0001.ED1"
                    
                print '\t DONE with GenerateInput (restart) W->W \n'
            else:
                print 'Skip running generateInput for W; as W fraction in plasma is ',  self.fluxFractionW

        ##TO DO: DO WE NEED THIS?!? IT IS FOR UPDATING THE PLASMA STATE -> MOVE IT TO TOP OR TO DRIVER?
        ##get name of FTridyn input file from config file to copy newly generated files to           
        #current_ftridyn_namelist = self.services.get_config_param('FTRIDYN_INPUT_FILE')
        ##this may be more than one file, not sure yet - need to learn more about FTridyn I/O
        #from_file_list = self.COPY_FILES.split()
        #file_list = current_ftridyn_namelist.split()

        ##copy newly generated files to names specified in config file
        #for index in range(0,1): #range(len(file_list)):
        #    print('copying ', from_file_list[index], ' to ', file_list[index])
        #    shutil.copyfile(from_file_list[index], file_list[index])

        self.services.update_plasma_state()


    def step(self, timeStamp=0.0,**keywords):
        print('ftridyn_worker: step')
        self.services.stage_plasma_state()

        print 'check that all arguments are read well by ftridyn-step'
        for (k, v) in keywords.iteritems():
            print '\t', k, " = ", v

        angleInputFile=keywords['gAngleInFile']
        spYieldsTemp=keywords['spYieldsFile_temp']
        spYieldsFinal=keywords['spYieldsFile_final']
       
        spYieldModeHe=keywords['fSpYieldModeHe']
        spYieldModeW=keywords['fSpYieldModeW']

        #RUN FTRIDYN for He->W
        #call shell script that runs FTridyn and pipes input file
        #task_id = self.services.launch_task(self.NPROC,
        #                                    self.services.get_working_dir(),
        #                                    self.FTRIDYN_EXE_He)

        #monitor task until complete
        #if (self.services.wait_task(task_id)):
        #    self.services.error('ftridyn_worker: step failed.')


        cwd = self.services.get_working_dir()
        ftOutputTempFile=cwd+'/'+self.OUTPUT_FTRIDYN_TEMP

        #keep all files to be saved (not plasma state) in folder with time stamp
        timeFolder=cwd+'/t'+str(self.driverTime)
        if not os.path.exists(timeFolder):
            os.makedirs(timeFolder)
        print '\n saving the output of this time-loop in ', timeFolder

        poolHe = self.services.create_task_pool('poolHe')
        #for i in range(len(energy)):
        for j in range(len(self.angleInHe)):
            if (self.weightAngleHe[j] == 0.0):
                    print '\t weight of angle ',  self.angleInHe[j], ' is ' , self.weightAngleHe[j] , ', so skip running F-Tridyn'
            elif (self.weightAngleHe[j] > 0.0):
                pathFolderHe = cwd+'/'+self.FOLDER_He+str(self.angleInHe[j])# + "_"+str(energyInW[j])
                poolInputFile=self.INPUT_FILE_He+'_angle'+str(self.angleInHe[j])
                shutil.copyfile(self.INPUT_FILE_He, pathFolderHe+'/'+'FTridyn.IN')
                self.services.add_task('poolHe', 'task'+str(self.angleInHe[j]), 1, pathFolderHe, self.FTRIDYN_EXE, poolInputFile,logfile=pathFolderHe+'/task.log' )
            

        ret_val = self.services.submit_tasks('poolHe')
        print '\t ret_val = ', ret_val
        exit_status = self.services.get_finished_tasks('poolHe')
        print exit_status
        self.services.remove_task_pool('poolHe')


        #post-processing: re-format output (fit function to profile) for Xolotl
        #how to get W sputtering yield from FT output using He_WSPYL.dat or He_WOUT.dat is described in 'note_calculateSpY'
        folderHe=cwd+'/'+self.FOLDER_He
        pathAngleFile_He=self.gitrOuputFolder_He+'/'+angleInputFile
        maxDepthHe=[]
        for j in range(len(self.angleInHe)):
            if (self.weightAngleHe[j] > 0.0):
                filePrjHe=cwd+'/'+self.FOLDER_He+str(self.angleInHe[j])+'/'+self.OUTPUT_PRJ_FILE_He
#                print 'loading file ', filePrjHe, 'to find maximum projectile depth'
                depth, bla=np.loadtxt(filePrjHe, usecols = (2,3) , unpack=True)
                #get maximum projectile range to ensure bins are added correctly in 'translate_ftridyn_to_xolotl'
                maxDepthHe.append(max(depth))
                
        print '\t maximum projectile range for He is ', max(maxDepthHe) , ' [A]'

        spYieldHe=translate_ftridyn_to_xolotl.ftridyn_to_xolotl(ftridynOnePrjOutput=self.OUTPUT_PRJ_FILE_He, ftridynOneOutOutput=self.OUTPUT_FILE_He, ftridynFolder=folderHe, fNImpacts=self.nImpacts, gAngleDistrib=pathAngleFile_He, angle=self.angleInHe, prjRange=max(maxDepthHe),nBins=keywords['xNGrid'])
        #overwrite spY value is mode is fixed
        if spYieldModeHe=='fixed':
            spYieldHe=keywords['fSpYieldHe']


        #append output and save to folder 
        tempfile = open(ftOutputTempFile,"r")
        f = open(self.OUTPUT_FTRIDYN_FINAL_He, "a")
        f.write(tempfile.read())
        f.close()
        tempfile.close()


        #keep copies of tridyn.dat
        shutil.copyfile(ftOutputTempFile, cwd+'/'+self.OUTPUT_FTRIDYN_TEMP_He)
        shutil.copyfile(ftOutputTempFile, timeFolder+'/'+self.OUTPUT_FTRIDYN_TEMP_He)
        #print "\t done with copying and appending tridyn.dat"

        #COPY FOLDERS TO TIME-STAMPED FOLDER 
        for j in range(len(self.angleInHe)):
            if (self.weightAngleHe[j] > 0.0):
                pathFolderHe = folderHe+str(self.angleInHe[j])# + "_"+str(energyInHe[j])
                shutil.copytree(pathFolderHe,timeFolder+'/'+self.FOLDER_He+str(self.angleInHe[j]))


        print '\t DONE with FTridyn  He->W \n'


        #RUN FTRIDYN for W->W

        if self.fluxFractionW>0.0:
            poolW = self.services.create_task_pool('poolW')
            for j in range(len(self.angleInW)): #for i in range(len(energy)):
                if (self.weightAngleW[j] == 0.0):
                    print '\t weight of angle ',  self.angleInW[j], ' is ' , self.weightAngleW[j] , ', so skip running F-Tridyn'
                elif (self.weightAngleW[j] > 0.0):
                    pathFolderW = cwd+'/'+self.FOLDER_W+str(self.angleInW[j])# + "_"+str(energyInW[j])
                    poolInputFile=self.INPUT_FILE_W+'_angle'+str(self.angleInW[j])
                    shutil.copyfile(self.INPUT_FILE_W, pathFolderW+'/'+'FTridyn.IN')
                    self.services.add_task('poolW', 'task'+str(self.angleInW[j]), 1, pathFolderW, self.FTRIDYN_EXE,poolInputFile,logfile=pathFolderW+'/task.log' )
            
            ret_val = self.services.submit_tasks('poolW')
            print '\t ret_val = ', ret_val
            exit_status = self.services.get_finished_tasks('poolW')
            print exit_status
            self.services.remove_task_pool('poolW')
        

            #post-processing: re-format output (fit function to profile) for Xolotl 
            #ftridynFolderPath='task'#+str(energy[i]) + "_"
            folderW=cwd+'/'+self.FOLDER_W
            pathAngleFile_W=self.gitrOuputFolder_W+'/'+angleInputFile
            maxDepthW=[]
            for j in range(len(self.angleInW)):
                if (self.weightAngleW[j] > 0.0):
                    filePrjW=cwd+'/'+self.FOLDER_W+str(self.angleInW[j])+'/'+self.OUTPUT_PRJ_FILE_W
                    #                print 'loading file ', filePrjW, 'to find maximum projectile depth'
                    depth, bla=np.loadtxt(filePrjW, usecols = (2,3) , unpack=True)
                    #get maximum projectile range to ensure bins are added correctly in 'translate_ftridyn_to_xolotl'
                    maxDepthW.append(max(depth))

            print '\t maximum projectile range for W is ', max(maxDepthW), ' [A]'

            spYieldW=translate_ftridyn_to_xolotl.ftridyn_to_xolotl(ftridynOnePrjOutput=self.OUTPUT_PRJ_FILE_W, ftridynOneOutOutput=self.OUTPUT_FILE_W, ftridynFolder=folderW, fNImpacts=self.nImpacts, gAngleDistrib=pathAngleFile_W, angle=self.angleInW, prjRange=max(maxDepthW), nBins=keywords['xNGrid'])
            #overwrite spY value is mode is fixed
            if spYieldModeW=='fixed':
                spYieldW=keywords['fSpYieldW']


            #append output (and save to what folder?)
            tempfile = open(ftOutputTempFile,"r")
            f = open(self.OUTPUT_FTRIDYN_FINAL_W, "a")
            f.write(tempfile.read())
            f.close()
            tempfile.close()

            #keep copies of tridyn.dat
            shutil.copyfile(ftOutputTempFile, cwd+'/'+self.OUTPUT_FTRIDYN_TEMP_W)
            shutil.copyfile(ftOutputTempFile, timeFolder+'/'+self.OUTPUT_FTRIDYN_TEMP_W)
        

            #COPY FOLDERS WITH TIME-STAMP
            #shutil.copyfile(pathFolderHe,pathFolderHe+timeStamp)
            for j in range(len(self.angleInW)):
                if (self.weightAngleW[j] > 0.0):
                    pathFolderW = folderW+str(self.angleInW[j])# + "_"+str(energyInW[j])
                    shutil.copytree(pathFolderW,timeFolder+'/'+self.FOLDER_W+str(self.angleInW[j]))

            print '\t DONE with FTridyn  W->W \n'
            
        else:
            print 'Skip running FTridyn for W, as fraction of W in plasma is', self.fluxFractionW
            spYieldW=0.0
            outputFTFileW=open(cwd+'/'+self.OUTPUT_FTRIDYN_TEMP_W, "w")
            outputFTFileW.write(" 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0")
            outputFTFileW.close()
        
        #write sputtering yields to file so they can be used by Xolotl
        print 'sputtering Yield due to He is ', spYieldHe
        print 'sputtering Yield due to W redeposition is ', spYieldW        

        #write sp Yields to file (temp) and append to spYield output (final)
        spTempfile = open(cwd+'/'+spYieldsTemp,"w+")
        spTempfile.write("%s %s %s \n" %(self.driverTime,spYieldHe,spYieldW))
        spFile = open(cwd+'/'+spYieldsFinal, "a+")
        spFile.write("%s %s %s \n" %(self.driverTime,spYieldHe,spYieldW))
        #spFile.write(spTempfile.read())
        #print 'spYieldsTemp: ', cwd+'/'+spYieldsTemp, ' reads: (THIS IS A TEST)'
        #for line in spTempfile:
        #    print line
        #    spFile.write(line)
        #print 'spYieldsFinal: ', cwd+'/'+spYieldsFinal ,' reads: (THIS IS A TEST)'
        #for line in spFile:
        #    print line
        spFile.close()
        spTempfile.close()

        shutil.copyfile(cwd+'/'+spYieldsTemp,timeFolder+'/'+spYieldsTemp)
        shutil.copyfile(cwd+'/'+spYieldsFinal,timeFolder+'/'+spYieldsFinal) #perhaps unnecessary

        #write format tridyn.dat to include W redep in Xolotl
        
        if os.path.exists(ftOutputTempFile):
            os.remove(ftOutputTempFile)
        combinedFile = open(ftOutputTempFile, 'a')
        outputFTFileHe=open(cwd+'/'+self.OUTPUT_FTRIDYN_TEMP_He, "r")
        combinedFile.write(outputFTFileHe.read())
        combinedFile.write("%s \n" %(str(self.fluxFractionW)))
        outputFTFileW=open(cwd+'/'+self.OUTPUT_FTRIDYN_TEMP_W, "r")
        combinedFile.write(outputFTFileW.read())

        combinedFile.close()
        outputFTFileHe.close()
        outputFTFileW.close()
        
        shutil.copyfile(ftOutputTempFile,timeFolder+'/'+self.OUTPUT_FTRIDYN_TEMP)

        #updates plasma state FTridyn output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
