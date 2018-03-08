#! /usr/bin/env python

from  component import Component
import gitr
import os
import pickle
import shutil

class gitr_comp(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        #Set up input deck
        print('input dir and cwd ', self.INPUT_DIR, ' ', os.getcwd())
        gitr.copy_folder(self.GITR_INPUT_DIR,os.getcwd())
        gitr.modifyInputParam(nT=int(self.NT))
        return

    def step(self, timeStamp=0.0):
        print 'Hello from gitr_comp'
        self.services.stage_plasma_state()
        shutil.copyfile('ftridyn.nc','input/ftridyn.nc')
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.GITR_EXE)
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('gitr_comp: step failed.')
        
        gitr.piscesProcessing()
        self.services.update_plasma_state()         
        
        return
    
    def finalize(self, timeStamp=0.0):
        return
    
