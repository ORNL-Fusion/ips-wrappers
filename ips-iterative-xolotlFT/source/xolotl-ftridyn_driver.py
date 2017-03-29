#! /usr/bin/env python

from  component import Component
import sys
import os


class xolotlFtridynDriver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: init')
        
        plasma_state_file = self.services.get_config_param('PLASMA_STATE_FILES')
        plasma_state_list = plasma_state_file.split()
        for index in range(len(plasma_state_list)):
            open(plasma_state_list[index], 'a').close()

        self.services.update_plasma_state()

    def step(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: step')

        ftridyn = self.services.get_port('WORKER')
        xolotl = self.services.get_port('XWORKER')

        self.services.stage_plasma_state() 

        #get config file and write initial state
        xolotl_ftridyn_config_file = self.services.get_config_param('PARAMETER_CONFIG_FILE')
        fid = open(xolotl_ftridyn_config_file, 'w')
        fid.write("mode='INIT'")
        fid.close()       

        self.services.update_plasma_state()
        sys.path.append(os.getcwd())

        for i in range(1,3):
            print 'loop %d' %(i)
            self.services.call(ftridyn, 'init', timeStamp)
            self.services.call(ftridyn, 'step', timeStamp)
            
            self.services.call(xolotl, 'init', timeStamp)
            self.services.call(xolotl, 'step', timeStamp)    
            
            self.services.stage_plasma_state() 
            
            import parameterConfig
            reload(parameterConfig)
            if (parameterConfig.mode == 'INIT'):
                fid = open(xolotl_ftridyn_config_file, 'w')
                fid.write("mode='RESTART'")
                fid.close()      

            self.services.update_plasma_state()

    def finalize(self, timeStamp=0.0):
        print('xolotl-ftridyn_driver: finalize')
        self.services.call(ftridyn, 'finalize', timeStamp)
