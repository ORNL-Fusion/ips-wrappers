#! /usr/bin/env python

from  component import Component
import os
import shutil
import subprocess
import glob
import translate_xolotl_to_ftridyn
import write_xolotl_paramfile
import sys

class xolotlWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        

    def init(self, timeStamp=0.0):
        print('xolotl_worker: init')
        self.services.stage_plasma_state()

        sys.path.append(os.getcwd())
        import parameterConfig
        reload(parameterConfig)

        print('from parameterConfig file, mode is %s', parameterConfig.mode)
        print('from parameterConfig file, driverTime is %s', parameterConfig.driverTime)
        print('from parameterConfig file, driverTimeStep is %s', parameterConfig.driverTimeStep)
        runEndTime=parameterConfig.driverTime+parameterConfig.driverTimeStep
        print('from parameterConfig file, this runs end time is  %f', runEndTime)

        if (parameterConfig.mode == 'INIT'):
            print('run xolotl preprocessor')
            #run prepocessor and copy params.txt input file to plasma state
            os.system('java -Djava.library.path=/project/projectdirs/atom/atom-install-edison/xolotl/xolotl-trunk-source/gov.ornl.xolotl.preprocessor/deps -cp .:/project/projectdirs/atom/atom-install-edison/xolotl/xolotl-trunk-source/gov.ornl.xolotl.preprocessor/deps/*:/project/projectdirs/atom/atom-install-edison/xolotl/xolotl-trunk-build/gov.ornl.xolotl.preprocessor/preprocessor/CMakeFiles/xolotlPreprocessor.dir/ gov.ornl.xolotl.preprocessor.Main --nxGrid 160 --maxVSize 250 --phaseCut')        
            write_xolotl_paramfile.writeXolotlParameterFile_fromPreprocessor(start_stop=runEndTime,ts_final_time=runEndTime)
            shutil.copyfile('params.txt','paramsInit.txt')  #store file to look into network issue; line to be removed
        else:
            write_xolotl_paramfile.writeXolotlParameterFile_fromTemplate(start_stop=runEndTime,ts_final_time=runEndTime,networkFile="xolotlStop.h5")
            shutil.copyfile('params.txt','paramsRestart.txt') #store file to look into network issue; line to be removed
        self.services.update_plasma_state()

    def step(self, timeStamp=0.0):
        print('xolotl_worker: step')
        self.services.stage_plasma_state()
        #call shell script that runs Xolotl and pipes input file
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.XOLOTL_EXE, 'params.txt')
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('xolotl_worker: step failed.')

        newest = max(glob.iglob('TRIDYN_*.dat'), key=os.path.getctime)
        print('newest file ' , newest)
        shutil.copyfile(newest, 'TRIDYN_last.dat')

        #append output
        tempfile = open(self.OUTPUT_XOLOTL_TEMP,"r")
        f = open(self.OUTPUT_XOLOTL_FINAL, "a")
        f.write(tempfile.read())
        f.close()
        tempfile.close()

        #updates plasma state Xolotl output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
