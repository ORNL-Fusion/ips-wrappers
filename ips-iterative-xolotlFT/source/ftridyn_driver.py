#! /usr/bin/env python

from  component import Component

class ftridynDriver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        print('ftridyn_driver: init')
        
        #Get info from config file about port 'CODE_INIT'
        self.ftridyn_init = self.services.get_port('CODE_INIT')

        #call init method of CODE_INIT ('ftridyn_init.py')
        self.services.call(self.ftridyn_init, 'init', timeStamp)
        
        #Get info from config file about port 'WORKER'
        self.ftridyn_comp = self.services.get_port('WORKER')

        #Get info from config file about 'XWORKER'
        self.xolotl_comp = self.services.get_port('XWORKER')

    def step(self, timeStamp=0.0):
        print('ftridyn_driver: step')
        self.services.stage_plasma_state() 
        #get config file and write initial state
        parameter_config_file = self.services.get_config_param('PARAMETER_CONFIG_FILE')
        fid = open(parameter_config_file, 'w')
        fid.write("mode='INIT'")
        fid.close()       
        self.services.update_plasma_state()
        
        self.services.call(self.xolotl_comp, 'init', timeStamp)
        ##call step method of CODE_INIT ('ftridyn_init.py')
        self.services.call(self.ftridyn_init, 'step', timeStamp)

        #call init method of WORKER ('ftridyn_worker.py')
        self.services.call(self.ftridyn_comp, 'init', timeStamp)
        ##call step method of WORKER ('ftridyn_worker.py')
        self.services.call(self.ftridyn_comp, 'step', timeStamp)

        self.services.call(self.xolotl_comp, 'step', timeStamp)    
        
        self.services.stage_plasma_state() 
        fid = open(parameter_config_file, 'w')
        fid.write("mode='RESTART'")
        fid.close()      
        self.services.update_plasma_state()
        self.services.call(self.ftridyn_comp, 'init', timeStamp)
        self.services.call(self.ftridyn_init, 'step', timeStamp) 
        self.services.call(self.ftridyn_comp, 'step', timeStamp)
        self.services.call(self.xolotl_comp, 'init', timeStamp)
        self.services.call(self.xolotl_comp, 'step', timeStamp)    
    def finalize(self, timeStamp=0.0):
        print('ftridyn_driver: finalize')
        self.services.call(self.ftridyn_comp, 'finalize', timeStamp)
