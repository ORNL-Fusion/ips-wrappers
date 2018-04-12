#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC init component. This wapper only takes a VMEC input file
#  and runs VMEC.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ZipState
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
#  VMEC init Component init method. This method prepairs the plasma state. Input
#  files can either be a new namelist input, a new plasma state, or both.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('vmec_init: init')

#  Get config filenames.
        self.current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        current_wout_file = 'wout_{}.nc'.format(self.current_vmec_namelist.replace('input.','',1))
        self.current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        
#  Stage input files. Remove any old files that might still exist.
        if os.path.exists(self.current_vmec_namelist):
            os.remove(self.current_vmec_namelist)
        if os.path.exists(self.current_vmec_state):
            os.remove(self.current_vmec_state)
        self.services.stage_input_files(self.INPUT_FILES)
        
#  Create plasma state from files. Input files can either be a new plasma state,
#  namelist input file or both. If both file were staged, replace the namelist
#  input file. If the namelist file is present remove the wout file and assume
#  that woutfile should be regenerated.
        with ZipState.ZipState(self.current_vmec_state, 'a') as zip_ref:
            if os.path.exists(self.current_vmec_namelist):
                zip_ref.write(self.current_vmec_namelist)
                zip_ref.remove(current_wout_file)

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
