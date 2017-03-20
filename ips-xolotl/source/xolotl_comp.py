#! /usr/bin/env python

from  component import Component
import os
import shutil

class xolotlWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        #get name of FTridyn executable from config file
        #This is actually a shell script that calls FTridyn and pipes
        #the input to the executable
        #self.ftridyn_exe = self.FTRIDYN_EXE

    def init(self, timeStamp=0.0):
        print('run xolotl preprocessor')
        #run prepocessor and copy params.txt input file to plasma state
        os.system('java -Djava.library.path=/project/projectdirs/atom/users/tyounkin/xolotl-trunk-source/gov.ornl.xolotl.preprocessor/deps -cp .:/project/projectdirs/atom/users/tyounkin/xolotl-trunk-source/gov.ornl.xolotl.preprocessor/deps/*:/project/projectdirs/atom/users/tyounkin/xolotl-trunk-build/gov.ornl.xolotl.preprocessor/preprocessor/CMakeFiles/xolotlPreprocessor.dir/ gov.ornl.xolotl.preprocessor.Main --nxGrid 160 --xStepSize 0.5 --dimensions 1')
        #self.services.stage_plasma_state()
        self.services.update_plasma_state()

    def step(self, timeStamp=0.0):
        print('ftridyn_worker: step')
        self.services.stage_plasma_state()
        #call shell script that runs FTridyn and pipes input file
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.XOLOTL_EXE, 'params.txt')
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('ftridyn_worker: step failed.')
        
        #updates plasma state FTridyn output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
