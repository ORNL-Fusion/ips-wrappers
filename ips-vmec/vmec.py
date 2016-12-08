#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC component. This wapper only takes a VMEC input file and
#  runs VMEC.
#
#-------------------------------------------------------------------------------

from component import Component
import os

#-------------------------------------------------------------------------------
#
#  VMEC Component Constructor
#
#-------------------------------------------------------------------------------
class vmec(Component):
    def __init__(self, services, config):
        print('vmec: Construct')
        Component.__init__(self, services, config)
        
        self.vmec_exe = self.VMEC_EXE

#-------------------------------------------------------------------------------
#
#  VMEC Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('vmec: init')
        self.current_vmec_namelist = self.services.get_config_param('CURRENT_VMEC_NAMELIST')
        self.services.stage_plasma_state()

#-------------------------------------------------------------------------------
#
#  VMEC Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('vmec: step')

        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.VMEC_EXE,
                                            self.current_vmec_namelist,
                                            logfile = 'vmec.log')
    
        if (self.services.wait_task(task_id)):
            self.services.error('vmec: step failed.')

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  VMEC Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('vmec: finalize')
