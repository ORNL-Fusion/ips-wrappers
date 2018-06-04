#! /usr/bin/env python

from  component import Component
import gitr

class gitr_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        try:
            self.worker_comp = self.services.get_port('WORKER')
        except Exception:
            self.services.exception('Error accessing worker component')
            raise
    
        #Get output file names from config file
        outfile_param = self.services.get_config_param('PLASMA_STATE_FILES')
        #split filenames into a list of strings
        outfile_list = outfile_param.split()
        #loop over file names and create dummy files in ftridynInit work area
        for index in range(len(outfile_list)):
            open(outfile_list[index], 'a').close()

        self.services.update_plasma_state()
        self.services.call(self.worker_comp, 'init', 0.0)
    
        return

    def step(self, timeStamp=0.0):
        print 'gitr_driver: beginning step call' 
        self.services.call(self.worker_comp, 'step', 0.0)
        print 'GITRDriver: finished worker call' 
        return
    def finalize(self, timeStamp=0.0):
        print 'GITR Driver Finalized' 
        return

