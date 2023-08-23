#! /usr/bin/env Python

from ipsframework import Component
import os
import shutil
import glob
import sys
#import translate_xolotl_to_ftridyn
#import translate_ftridyn_to_xolotl
from ips_xolotlFT.python_scripts_for_coupling import generate_ftridyn_input #new script to generate FTridyn input
import numpy as np
import subprocess
import re
import pickle
import math

class ftridynWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))
        #get name of FTridyn executable from config file
        #This is actually a shell script that calls FTridyn and pipe the input to the executable
        self.ftridyn_exe = self.FTRIDYN_EXE

    def init(self, timeStamp=0.0,**keywords):

        cwd = self.services.get_working_dir()

        #stage plasma state files for use on execution of FTridyn
        self.services.stage_state()

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

        print_test=keywords['print_test']
        
        if 'output_file' in keywords:
            outFile=keywords['output_file']
            if outFile  is not None:
                print('\t redirect F-TRIDYN:init output')
                print('\t \t of ', cwd )
                print('\t \t to:', outFile)
                outF=open(outFile , 'a')
                sys.stdout = outF
            else:
                print ('\t not log file defined in keywords or config file')
                print ('\t print output of F-TRIDYN:init to default sys.stdout')
        sys.stdout.flush()
        print(' ')
        print('fridyn_worker: init : for ', str(prj))
        print(' ')

        #print lines if in test mode:
        if print_test:
            print('\t check that all arguments are read well by ftridyn-init \n')
            for (k, v) in keywords.items():
                print(('\t \t {0} = {1}'.format(k,v)))

        if os.path.exists(self.ft_folder):
            shutil.rmtree(self.ft_folder)
        os.makedirs(self.ft_folder)

        #run new script to generate FTridyn input for He->W 
        for j in range(len(angleIn)):

            if (weightAngle[j] == 0.0):
                print(('\t weight of angle {0} is {1}, so skip generateInput'.format(angleIn[j],weightAngle[j]))) 

            elif (weightAngle[j] > 0.0):
                opts = {}
                for opt in ["ED_He", "EF_He", "E0_He", "ALPHA0_He", "ED_W", "EF_W", "SBV_W"]:
                    if opt in ftridyn.keys():
                        opts[opt] = ftridyn[opt]
                if print_test:
                    print('The generate_ftridyn_input path in ftridyn_comp is:')
                    print(os.path.abspath(generate_ftridyn_input.__file__)) 
                generate_ftridyn_input.Prj_Tg_xolotl(IQ0=ftridyn['iQ0'],number_layers=ftridyn['nDataPts'],depth=ftridyn['nTT'],number_histories=ftridyn['nImpacts'],incident_energy=energyIn,incident_angle=angleIn[j],projectile_name=str(prj),target1_name=str(tg[0]), target2_name=str(tg[1]),target3_name=str(tg[2]),target4_name=str(tg[3]), **opts)
                sys.stdout.flush()
                
            pathFolder = self.ft_folder+'/ANGLE'+str(angleIn[j])# + "_"+str(energyInW[j])
            #if not os.path.exists(pathFolder):
            os.makedirs(pathFolder)            
            shutil.copyfile(fInFile, pathFolder+'/'+'FTridyn.IN')
            
            if otherInFiles:
                for f in otherInFiles:
                    shutil.copyfile(f, pathFolder+"/"+f)
                    
            if energyIn < 0:
                shutil.copyfile(origEnergyFilePath+'/'+origEnergyFilePattern[0]+str(int(j))+origEnergyFilePattern[1],pathFolder+"/"+energyFileName)
                
        print('\t DONE with GenerateInput \n')

        sys.stdout.flush()
        self.services.update_state()

        
    def step(self, timeStamp=0.0,**keywords):

        cwd = self.services.get_working_dir()
        self.services.stage_state()

        energyIn=keywords['fEnergyIn']
        angleIn = keywords['fAngleIn']
        weightAngle=keywords['fWeightAngle']
        ftridyn=keywords['ftParameters']
        prj=keywords['fPrj']
        #tg=keywords['fTg']
        print_test=keywords['print_test']
        
        if 'output_file' in keywords:
            outFile=keywords['output_file']
            if outFile  is not None:
                print('\t redirect F-TRIDYN:step output')
                print('\t \t of ', cwd )
                print('\t \t to:', outFile)
                outF=open(outFile , 'a')
                sys.stdout = outF
            else:
                print ('\t no log file defined in keywords or config file')
                print ('\t print output of F-TRIDYN:step to default sys.stdout')

        print(' ')
        print('ftridyn_worker: step : for ', str(prj))
        print(' ')

        #if test:
        if print_test:
            print('\t TEST: check that all arguments are read well by ftridyn-step \n')
            for (k, v) in keywords.items():
                print(('\t \t {0} = {1} '.format(k,v)))

        #call shell script that runs FTridyn and pipes input file
        #task_id = self.services.launch_task(self.NPROC,
        #                                    self.services.get_working_dir(),
        #                                    self.FTRIDYN_EXE_He)

        #monitor task until complete
        #if (self.services.wait_task(task_id)):
        #    self.services.error('ftridyn_worker: step failed.')


        #RUN FTRIDYN
        os.environ['OMP_NUM_THREADS']=self.THREADS_PER_TASK
        file_namelist = []
        #pool_ftx = self.services.create_task_pool('pool_ftx')
        for j in range(len(angleIn)): #for i in range(len(energy)):

            if (weightAngle[j] == 0.0):
                print(('\t weight of angle {0} is {1}, so skip running F-Tridyn'.format(angleIn[j],weightAngle[j])))

            elif (weightAngle[j] > 0.0):
                pathFolder = self.ft_folder+'/ANGLE'+str(angleIn[j])# + "_"+str(energyIn[j])
                file_namelist.append(pathFolder)
                #poolInput_ftx='FTridyn_angle'+str(angleIn[j])
                #self.services.add_task('pool_ftx', 'task'+str(angleIn[j]), 1, pathFolder, self.FTRIDYN_EXE,poolInput_ftx,logfile=pathFolder+'/task.log' )
                 
        nFTruns = len(file_namelist)
        nFTrunsPerNode = int(self.FTMPI_PPN)
        nFTrunNodes = int(self.FTMPI_NODES)
        if(nFTruns <= nFTrunsPerNode*nFTrunNodes): 
            nFTrunsPerNode = int(math.ceil(1.0*nFTruns/nFTrunNodes))
            nFTpoolTasks = nFTrunNodes
        else:
            nFTpoolTasks = int(math.ceil(1.0*nFTruns/nFTrunsPerNode))

        with open('ftridyn_file_namelist.pkl', 'wb') as f:
            pickle.dump(list(reversed(file_namelist)), f)

        pool_ftx = self.services.create_task_pool('pool_ftx')
        for i in range(nFTpoolTasks):
            task_pool_log='task_pool'+str(i)+'.log'
            #poolInput_ftx='FTridyn_angle'+str(angleIn[j])
            self.services.add_task('pool_ftx', 'task'+str(i), nFTrunsPerNode, cwd, 'python',self.FTMPI_EXEC,str(i),str(nFTrunsPerNode),str(self.FTRIDYN_EXE),"FTridyn.IN",task_ppn= nFTrunsPerNode,logfile=task_pool_log )
        ret_val = self.services.submit_tasks('pool_ftx')
        print(' ')
        print(('\t ret_val = {}'.format(ret_val)))
        exit_status = self.services.get_finished_tasks('pool_ftx')
        print(('\t exit status {} \n'.format(exit_status)))
        self.services.remove_task_pool('pool_ftx')

        for i in range(nFTpoolTasks):
            task_pool_log='task_pool'+str(i)+'.log'
            print('\t keep a copy of ',task_pool_log,' in every loop by moving to ',self.ft_folder)
            shutil.move(task_pool_log,self.ft_folder+'/'+task_pool_log)

            
        #write the path to the current working directory, so the driver can access the data
        outputPathString="outputPath="+cwd
        outputPathFile = open(self.PWD_PATH, "w")
        outputPathFile.write(outputPathString)
        outputPathFile.write('\n')
        outputPathFile.close()
        
        print('\t path of FTRIDYNs output:')
        print('\t \t {} '.format(cwd))
        print('\t written to file: ')
        print('\t \t{} '.format(self.PWD_PATH))
        print(' ')

        
        #post-processing: COMPRESS ALL OUTPUT INTO A GENERIC, SINGLE ZIP FILE AND ADD TO PLASMA STATE 
        #zippedFile=str(self.ft_folder)+'.zip'
        #if os.path.isfile(zippedFile):
        #    os.remove(zippedFile)

        #zip_output='zipOutput.txt'
        #zipString='zip -r %s %s >> %s' %(zippedFile,self.ft_folder,zip_output)
        #subprocess.call([zipString], shell=True)

        #shutil.rmtree(self.ft_folder)

        #updates plasma state FTridyn output files
        sys.stdout.flush()
        self.services.update_state()

    def finalize(self, timeStamp=0.0):
        return
