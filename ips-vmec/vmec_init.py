#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC init component. This wapper only takes a VMEC input file
#  and runs VMEC.
#
#-------------------------------------------------------------------------------

from component import Component
import os

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
#  VMEC init Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('vmec_init: init')
    
        self.services.stage_input_files(self.INPUT_FILES)
        
        current_vmec_namelist = self.services.get_config_param('CURRENT_VMEC_NAMELIST')
        
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            'mv',
                                            self.INPUT_FILES,
                                            current_vmec_namelist,
                                            logfile = 'vmec_init.log')
        if (self.services.wait_task(task_id)):
            self.services.error('vmec_init: init failed.')

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
