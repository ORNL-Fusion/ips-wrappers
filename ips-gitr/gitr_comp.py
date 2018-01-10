#! /usr/bin/env python

from  component import Component
import gitr
import os
import pickle

class gitr_comp(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        #get plasma state files
        self.services.update_plasma_state()         
        #Set up input deck
        gitr.copy_folder(self.INPUT_DIR,os.getcwd())
        return

    def step(self, timeStamp=0.0):
        print 'Hello from gitr_comp'
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.GITR_EXE)
        #monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('gitr_comp: step failed.')
        return
    
    def finalize(self, timeStamp=0.0):
        return
    
