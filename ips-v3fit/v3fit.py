#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for V3FIT component. This wapper only takes a V3FIT input file
#  and runs V3FIT.
#
#-------------------------------------------------------------------------------

from component import Component
import fileinput

#-------------------------------------------------------------------------------
#
#  V3FIT Component Constructor
#
#-------------------------------------------------------------------------------
class v3fit(Component):
    def __init__(self, services, config):
        print('v3fit: Construct')
        Component.__init__(self, services, config)
        
        self.v3fit_exe = self.V3FIT_EXE

#-------------------------------------------------------------------------------
#
#  V3FIT Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('v3fit: init')
        
        self.services.stage_input_files(self.INPUT_FILES)
        
        self.current_v3fit_namelist = self.services.get_config_param('CURRENT_V3FIT_NAMELIST')
        self.services.stage_plasma_state()
    
#-------------------------------------------------------------------------------
#
#  V3FIT Component step method. This runs V3FIT.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('v3fit: step')

        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.V3FIT_EXE,
                                            self.current_v3fit_namelist,
                                            logfile = 'v3fit.log')
    
        if (self.services.wait_task(task_id)):
            self.services.error('v3fit: step failed.')

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  V3FIT Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('v3fit: finalize')
