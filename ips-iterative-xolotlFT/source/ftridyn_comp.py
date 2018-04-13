#! /usr/bin/env python

from  component import Component
import os
import shutil
import glob
import sys
#import translate_xolotl_to_ftridyn
#import translate_ftridyn_to_xolotl
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
        print '\n'
        print('fridyn_worker: init')
        #stage plasma state files for use on execution of FTridyn
        self.services.stage_plasma_state()

        print 'check that all arguments are read well by ftridyn-init'
        for (k, v) in keywords.iteritems():
            print '\t', k, " = ", v

        #asign a local variable to arguments used multiple times
        driverTime=keywords['dTime']

        ftridyn=keywords['ftParameters']

        energyIn=keywords['fEnergyIn']
        angleIn = keywords['fAngleIn']
        weightAngle=keywords['fWeightAngle']
        prj=keywords['fPrj']
        tg=keywords['fTargetList']

        self.ft_folder=keywords['ft_folder']
        fInFile=keywords['input_file']

        otherInFiles=keywords['otherInFiles']

        energyFileName=keywords['energy_file_name']
        origEnergyFilePath=keywords['orig_energy_files_path']
        origEnergyFilePattern=keywords['orig_energy_files_pattern']

            
        #run new script to generate FTridyn input for He->W  
    
        if os.path.exists(self.ft_folder):
            shutil.rmtree(self.ft_folder)

        os.makedirs(self.ft_folder)

        for j in range(len(angleIn)):
            if (weightAngle[j] == 0.0):
                print '\t weight of angle ',  angleIn[j], ' is ' , weightAngle[j] , ', so skip generateInput'
            elif (weightAngle[j] > 0.0):
                generate_ftridyn_input.Prj_Tg_xolotl(IQ0=ftridyn['iQ0'],number_layers=ftridyn['nDataPts'],depth=ftridyn['nTT'],number_histories=ftridyn['nImpacts'],incident_energy=energyIn,incident_angle=angleIn[j],projectile_name=str(prj),target1_name=str(tg[0]), target2_name=str(tg[1]),target3_name=str(tg[2]),target4_name=str(tg[3]))
            pathFolder = self.ft_folder+'/ANGLE'+str(angleIn[j])# + "_"+str(energyInW[j])
            #if not os.path.exists(pathFolder):
            os.makedirs(pathFolder)
            
            shutil.copyfile(fInFile, pathFolder+'/'+'FTridyn.IN')

            if otherInFiles:
                for f in otherInFiles:
                    shutil.copyfile(f, pathFolder+"/"+f)

            if energyIn < 0:
                shutil.copyfile(origEnergyFilePath+'/'+origEnergyFilePattern[0]+str(int(j))+origEnergyFilePattern[1],pathFolder+"/"+energyFileName)


        print '\t DONE with GenerateInput \n'

        self.services.update_plasma_state()


    def step(self, timeStamp=0.0,**keywords):
        print('ftridyn_worker: step')
        self.services.stage_plasma_state()

        cwd = self.services.get_working_dir()

        print 'check that all arguments are read well by ftridyn-step'
        for (k, v) in keywords.iteritems():
            print '\t', k, " = ", v

        energyIn=keywords['fEnergyIn']
        angleIn = keywords['fAngleIn']
        weightAngle=keywords['fWeightAngle']
        #tg=keywords['fTg']
        #prj=keywords['fPrj']

        #call shell script that runs FTridyn and pipes input file
        #task_id = self.services.launch_task(self.NPROC,
        #                                    self.services.get_working_dir(),
        #                                    self.FTRIDYN_EXE_He)

        #monitor task until complete
        #if (self.services.wait_task(task_id)):
        #    self.services.error('ftridyn_worker: step failed.')


        #RUN FTRIDYN


        pool_ftx = self.services.create_task_pool('pool_ftx')
        for j in range(len(angleIn)): #for i in range(len(energy)):
            if (weightAngle[j] == 0.0):
                print '\t weight of angle ',  angleIn[j], ' is ' , weightAngle[j] , ', so skip running F-Tridyn'
            elif (weightAngle[j] > 0.0):
                pathFolder = self.ft_folder+'/ANGLE'+str(angleIn[j])# + "_"+str(energyIn[j])
                poolInput_ftx='FTridyn_angle'+str(angleIn[j])
                self.services.add_task('pool_ftx', 'task'+str(angleIn[j]), 1, pathFolder, self.FTRIDYN_EXE,poolInput_ftx,logfile=pathFolder+'/task.log' )
            
        ret_val = self.services.submit_tasks('pool_ftx')
        print '\t ret_val = ', ret_val
        exit_status = self.services.get_finished_tasks('pool_ftx')
        print exit_status
        self.services.remove_task_pool('pool_ftx')
        
        #write the path to the current working directory, so the driver can access the data
        outputPathString="outputPath="+cwd
        outputPathFile = open(self.PWD_PATH, "w")
        outputPathFile.write(outputPathString)
        outputPathFile.write('\n')
        outputPathFile.close()

        print 'path of FTRIDYNs output: ', cwd
        print '\t written to file: ', self.PWD_PATH , '\n'

        #post-processing: COMPRESS ALL OUTPUT INTO A GENERIC, SINGLE ZIP FILE AND ADD TO PLASMA STATE 
        #zippedFile=str(self.ft_folder)+'.zip'
        #if os.path.isfile(zippedFile):
        #    os.remove(zippedFile)

        #zip_output='zipOutput.txt'
        #zipString='zip -r %s %s >> %s' %(zippedFile,self.ft_folder,zip_output)
        #subprocess.call([zipString], shell=True)

        #shutil.rmtree(self.ft_folder)

        #updates plasma state FTridyn output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
