#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC init component. This wapper only takes a VMEC input file
#  and runs VMEC.
#
#-------------------------------------------------------------------------------

from component import Component
import shutil

#-------------------------------------------------------------------------------
#
#  VMEC init Component Constructor
#
#-------------------------------------------------------------------------------
class vmec_init(Component):
    def __init__(self, services, config):
        print('vmec_init: Construct')
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  VMEC init Component init method. This method prepairs the namelist input
#  file and creates a dummy out put file. This allows staging the plasma state
#  files.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('vmec_init: init')

        self.services.stage_input_files(self.INPUT_FILES)
        
        current_vmec_namelist = self.services.get_config_param('CURRENT_VMEC_NAMELIST')
        shutil.copyfile(self.INPUT_FILES, current_vmec_namelist)
        
#  Create a dummy wout file so the plasma state has something to update to.
        open(self.services.get_config_param('CURRENT_VMEC_WOUT_FILE'), 'a').close()
        
        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  VMEC init Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('vmec_init: step')

#-------------------------------------------------------------------------------
#
#  VMEC init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('vmec_init: finalize')
