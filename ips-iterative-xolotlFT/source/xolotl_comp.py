#! /usr/bin/env python

from  component import Component
import os
import shutil
import subprocess
import glob
import param_handler #write_xolotl_paramfile
import sys
import numpy as np

class xolotlWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        

    def init(self, timeStamp=0.0, **keywords):
        print('xolotl_worker: init')
        self.services.stage_plasma_state()

        print 'check that all arguments are read well by xolotl-init and write Xolotl input file (from dictionary)' 
        for (k, v) in keywords.iteritems():
            print '\t', k, " = ", v

        #asign a local variable to arguments used multiple times 
        self.driverTime=keywords['dTime']
        self.coupling=keywords['xFtCoupling']

        cwd = self.services.get_working_dir()

        xp = param_handler.xolotl_params()
        xp.parameters=keywords['xParameters'] 

        #write and store xolotls parameter for each loop 
        xp.write('params.txt')
        currentXolotlParamFile='params_%f.txt' %self.driverTime
        shutil.copyfile('params.txt',currentXolotlParamFile) 
        
        self.services.update_plasma_state()

    def step(self, timeStamp=0.0,**keywords):
        print('xolotl_worker: step')

        self.services.stage_plasma_state()

        #asign a local variable to arguments used multiple times

        print 'check that all arguments are read well by xolotl-step'
        for (k, v) in keywords.iteritems():
            print '\t', k, " = ", v
            
        zipOutput=keywords['dZipOutput']
        petscHeConc=keywords['xHe_conc']
        xp_parameters=keywords['xParameters']

        #xolotlLogFile='xolotl_t%f.log' %self.driverTime
        #print '\t Xolotl log file ', xolotlLogFile

        #call shell script that runs Xolotl and pipes input file
        #task_id = self.services.launch_task(self.NPROC,
        #                                    self.services.get_working_dir(),
        #                                    self.XOLOTL_EXE, 'params.txt', 
        #                                    logfile=xolotlLogFile)

        #monitor task until complete
        #if (self.services.wait_task(task_id)):
        #    self.services.error('xolotl_worker: step failed.')


        import time
        num_trials = 2 #SET AS DRIVER ARGUMENT
        for i in range(num_trials):
            xolotlLogFile='xolotl_t%f_%d.log' %(self.driverTime,i)
            task_id = self.services.launch_task(self.NPROC,self.services.get_working_dir(),
                                                self.XOLOTL_EXE, 'params.txt',logfile=xolotlLogFile)
            ret_val = self.services.wait_task(task_id)

            if (ret_val == 0):
                break
            else:
                self.services.error('xolotl_worker: step failed in trial %d.' %i)
                time.sleep(5)
                #if it failed, save last networkFile before a new try and set newest network file to use in the
                # next try, so that it starts from the last saved time step, not from the beginning of the loop
                shutil.copyfile(xp_parameters['networkFile'],'networkFile_%f_%d.h5' %(self.driverTime,i))
                shutil.copyfile('xolotlStop.h5',xp_parameters['networkFile'])

        else:
            self.services.error('xolotl_worker: Aborting after %d num_trials)trials' %num_trials)
            raise Exception("Aborting simulation after %d failed xolotl runs" % num_trials)

        #ALREADY IN DRIVER save network file with a different name to use in the next time step
        #currentXolotlNetworkFile='xolotlStop_%f.h5' %self.driverTime
        #shutil.copyfile('xolotlStop.h5',currentXolotlNetworkFile)
        #shutil.copyfile('xolotlStop.h5',xp.parameters['networkFile'])

        newest = max(glob.iglob('TRIDYN_*.dat'), key=os.path.getctime)
        print '\t newest file ' , newest
        shutil.copyfile(newest, 'last_TRIDYN.dat')


        #save TRIDYN_*.dat files, zipped
        TRIDYNFiles='TRIDYN_*.dat'

        if self.coupling=='True' and zipOutput=='True':
            TRIDYNZipped='allTRIDYN_t%f.zip' %self.driverTime
            zip_ouput='zipTridynDatOuput.txt'

            print '\t save and zip output: ', TRIDYNFiles
            zipString='zip %s %s >> %s ' %(TRIDYNZipped, TRIDYNFiles, zip_ouput)
            subprocess.call([zipString], shell=True)

            rmString='rm '+ TRIDYNFiles
            subprocess.call([rmString], shell=True)

        elif self.coupling=='True':
            print '\t leaving ',  TRIDYNFiles , 'uncompressed'

        else:
            print '\t no ', TRIDYNFiles , ' generated in this simulation'

        #save helium concentration files, zipped
        heConcFiles='heliumConc_*.dat'

        if petscHeConc and zipOutput=='True':
            heConcZipped='allHeliumConc_t%f.zip' %self.driverTime
            zip_ouput='zipHeConcOuput.txt'
            
            print '\t save and zip output: ', heConcFiles
            zipString='zip %s %s >> %s ' %(heConcZipped, heConcFiles, zip_ouput)
            subprocess.call([zipString], shell=True)

            rmString='rm '+heConcFiles
            subprocess.call([rmString], shell=True)
            
        elif petscHeConc:
            print '\t leaving ', heConcFiles ,'uncompressed'

        else:
            print '\t no ', heConcFiles , ' in this loops output'

        statusFile=open(self.EXIT_STATUS, "r")
        exitStatus=statusFile.read().rstrip('\n')
        
        print '\n'

        if exitStatus=='collapsed':
            print '\t simulation exited loop with status collapse'
            print '\t rename output files as _collapsed before trying again'

            currentXolotlNetworkFile='xolotlStop_%f.h5' %self.driverTime
            networkFile_unfinished='xolotlStop_%f_collapsed.h5' %self.driverTime
            os.rename(currentXolotlNetworkFile,networkFile_unfinished)

            retentionFile = self.RET_FILE
            rententionUnfinished = 'retention_t%f_collapsed.out' %self.driverTime
            shutil.copyfile(retentionFile,rententionUnfinished)
            
            surfaceFile=self.SURFACE_FILE
            surfaceUnfinished='surface_t%f_collapsed.txt' %self.driverTime
            shutil.copyfile(surfaceFile,surfaceUnfinished)

            if zipOutput=='True':
                TRIDYNZipped='allTRIDYN_t%f.zip' %self.driverTime
                TRIDYNUnfinished='allTRIDYN_t%f_collapsed.zip' %self.driverTime
                os.rename(TRIDYNZipped,TRIDYNUnfinished)
                if petscHeConc:
                    heConcZipped='allHeliumConc_t%f.zip' %self.driverTime
                    heConcUnfinished='allHeliumConc_t%f_collapsed.zip' %self.driverTime
                    os.rename(heConcZipped,heConcUnfinished)
            
        #updates plasma state Xolotl output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
