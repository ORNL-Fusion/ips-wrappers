#! /usr/bin/env python

from ipsframework import Component
import os
import shutil
import sys
import translate_xolotl_to_lay

class ftridynWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        #get name of FTridyn executable from config file
        #This is actually a shell script that calls FTridyn and pipes
        #the input to the executable
        self.ftridyn_exe = self.FTRIDYN_EXE

    def init(self, timeStamp=0.0):
        #stage plasma state files for use on execution of FTridyn
        self.services.stage_plasma_state()
        sys.path.append(os.getcwd())
        import parameterConfig
        reload(parameterConfig)

        if (parameterConfig.mode == 'RESTART'):
            nDataPts = translate_xolotl_to_lay.xolotlToLay()
            self.services.update_plasma_state()
    def step(self, timeStamp=0.0):
        print('ftridyn_worker: step')
        self.services.stage_plasma_state()
        #call shell script that runs FTridyn and pipes input file
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.FTRIDYN_EXE)
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('ftridyn_worker: step failed.')

        os.system(' '.join(['python', self.POSTPROCESSING_SCRIPT]))

        #append output
        tempfile = open(self.OUTPUT_FTRIDYN_TEMP,"r")
        f = open(self.OUTPUT_FTRIDYN_FINAL, "a")
        f.write(tempfile.read())
        f.close()
        tempfile.close()

        #updates plasma state FTridyn output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
