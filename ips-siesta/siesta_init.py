#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC init component. This wapper only takes a SIESTA input
#  file and runs SIESTA.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ZipState
import os

#-------------------------------------------------------------------------------
#
#  SIESTA init Component Constructor
#
#-------------------------------------------------------------------------------
class siesta_init(Component):
    def __init__(self, services, config):
        print('siesta_init: Construct')
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  SIESTA init Component init method. This method prepairs the namelist input
#  file and creates a dummy out put file. This allows staging the plasma state
#  files.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('siesta_init: init')

#  Get config filenames.
        current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        self.current_siesta_namelist = self.services.get_config_param('SIESTA_NAMELIST_INPUT')
        self.current_siesta_state = self.services.get_config_param('CURRENT_SIESTA_STATE')

#  Stage input files. Remove an old namelist input if it exists.
        if os.path.exists(self.current_siesta_namelist):
            os.remove(self.current_siesta_namelist)
        if os.path.exists(current_vmec_namelist):
            os.remove(current_vmec_namelist)
        self.services.stage_input_files(self.INPUT_FILES)
    
#  Create plasma state from files. Input files can either be a new plasma state,
#  namelist input file or both. If both file were staged, replace the namelist
#  input file. If the namelist file is present flag the plasma state as needing
#  to be updated.
        with ZipState.ZipState(self.current_siesta_state, 'a') as zip_ref:
            if os.path.exists(self.current_siesta_namelist):
                zip_ref.write(self.current_siesta_namelist)
                zip_ref.set_state(state='needs_update')

            with ZipState.ZipState(current_vmec_state, 'a') as zip_vmec_ref:
                if os.path.exists(current_vmec_namelist):
                    zip_vmec_ref.write(current_vmec_namelist)
                    zip_vmec_ref.set_state(state='needs_update')
            zip_ref.write(current_vmec_state)
            
        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  SIESTA init Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('siesta_init: step')

#-------------------------------------------------------------------------------
#
#  SIESTA init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('siesta_init: finalize')
