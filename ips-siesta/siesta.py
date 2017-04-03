#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for SIESTA component. This wapper only takes a SIESTA input file
#  and runs SIESTA.
#
#-------------------------------------------------------------------------------

from component import Component
import os

#-------------------------------------------------------------------------------
#
#  SIESTA Component Constructor
#
#-------------------------------------------------------------------------------
class siesta(Component):
    def __init__(self, services, config):
        print('siesta: Construct')
        Component.__init__(self, services, config)
        
        self.siesta_exe = self.SIESTA_EXE

#-------------------------------------------------------------------------------
#
#  SIESTA Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('siesta: init')
        self.services.stage_plasma_state()

#-------------------------------------------------------------------------------
#
#  SIESTA Component step method. This runs siesta.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('siesta: step')

        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.SIESTA_EXE,
                                            logfile = 'siesta.log')
    
        if (self.services.wait_task(task_id)):
            self.services.error('siesta: step failed.')

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  SIESTA Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('siesta: finalize')
