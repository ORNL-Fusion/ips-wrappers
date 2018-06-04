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
        current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        
#  Stage input files. Remove old namelist input if it exists.
        if os.path.exists(current_vmec_namelist):
            os.remove(current_vmec_namelist)
        if os.path.exists(current_vmec_namelist):
            os.remove(current_vmec_state)
        
        self.services.stage_input_files(self.INPUT_FILES)
        
#  Create plasma state from files. Input files can either be a new plasma state,
#  namelist input file or both. If both file were staged, replace the namelist
#  input file. If the namelist file is present flag the plasma state as needing
#  to be updated.
        with ZipState.ZipState(current_vmec_state, 'a') as zip_ref:
            if os.path.exists(current_vmec_namelist):
                zip_ref.write(current_vmec_namelist)
                zip_ref.set_state(state='needs_update')
                    
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
