#! /usr/bin/env python

from  component import Component
import os
import numpy as np

class ftridynDriver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        print('ftridyn_driver: init')
        current_ftridyn_namelist = self.services.get_config_param('PLASMA_STATE_FILES')
        ##split filenames into a list of strings
        file_list = current_ftridyn_namelist.split()
        ##loop over file names and create dummy files in ftridynInit work area
        for index in range(len(file_list)):
            if os.path.isfile(file_list[index]) == 0 :
                open(file_list[index], 'a').close()


        #update plasma state from relevant files in ftridynInit work area
        self.services.update_plasma_state()        

        #Get info from config file about port 'WORKER'
        self.ftridyn_comp = self.services.get_port('WORKER')

    def step(self, timeStamp=0.0):
        print('ftridyn_driver: step')
        
        self.services.call(self.ftridyn_comp, 'step', timeStamp)
    
    def finalize(self, timeStamp=0.0):
        print('ftridyn_driver: finalize')
        self.services.call(self.ftridyn_comp, 'finalize', timeStamp)
