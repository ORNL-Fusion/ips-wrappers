#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper to bridge the DAKOTA init component. This wapper only initalizes
#  empty dakota file.
#
#-------------------------------------------------------------------------------

from component import Component
import shutil

#-------------------------------------------------------------------------------
#
#  DAKOTA init Component Constructor
#
#-------------------------------------------------------------------------------
class dakota_init(Component):
    def __init__(self, services, config):
        print('dakota_init: Construct')
        Component.__init__(self, services, config)
    
#-------------------------------------------------------------------------------
#
#  DAKOTA init Component init method. This method reads the
#  parameter file from datoka and adds to an existing namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('dakota_init: init')
        
        if (self.INPUT_FILES != ''):
            self.services.stage_input_files(self.INPUT_FILES)
        
            current_dakota_param_file = self.services.get_config_param('CURRENT_DAKOTA_PARAM_FILE')
            shutil.copyfile(self.INPUT_FILES, current_dakota_param_file)

#  Create a dummy result file so the plasma state has something to update to.
            open(self.services.get_config_param('CURRENT_DAKOTA_COST_FILE'), 'a').close()
        
            self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  DAKOTA init Component step method. Not used.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('dakota_init: step')
    
#-------------------------------------------------------------------------------
#
#  DAKOTA init Component finalize method. This cleans up afterwards. Not
#  used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('dakota_init: finalize')
