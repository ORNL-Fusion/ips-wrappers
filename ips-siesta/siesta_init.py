#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC init component. This wapper only takes a SIESTA input
#  file and runs SIESTA.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ZipState
from utilities import ScreenWriter
import os

#-------------------------------------------------------------------------------
#
#  SIESTA init Component Constructor
#
#-------------------------------------------------------------------------------
class siesta_init(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  SIESTA init Component init method. This method prepairs the namelist input
#  file and creates a dummy out put file. This allows staging the state files.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'siesta_init: init')

#  Get config filenames.
        current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        current_siesta_namelist = self.services.get_config_param('SIESTA_NAMELIST_INPUT')
        current_siesta_state = self.services.get_config_param('CURRENT_SIESTA_STATE')

#  Remove old inputs. Stage input files.
        for file in os.listdir('.'):
            os.remove(file)
                        
        self.services.stage_input_files(self.INPUT_FILES)

#  Create a vmec state. If the vmec namelist file exists add the namelist input
#  file.
        with ZipState.ZipState(current_vmec_state, 'a') as zip_ref:
            if os.path.exists(current_vmec_namelist):
                zip_ref.write(current_vmec_namelist)
                zip_ref.set_state(state='needs_update')

#  Create state from files. Input files can either be a new state, namelist
#  input file or both. If both files were staged, replace the namelist input
#  file. If the namelist file is present flag the state as needing to be
#  updated. Sub states will automatically merge.
        with ZipState.ZipState(current_siesta_state, 'a') as zip_ref:
            if os.path.exists(current_siesta_namelist):
                zip_ref.write(current_siesta_namelist)
                zip_ref.set_state(state='needs_update')

#  The vmec state will be merged with any existing vmec state in the siesta
#  state.
            zip_ref.write(current_vmec_state)
            
        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  SIESTA init Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'siesta_init: step')

#-------------------------------------------------------------------------------
#
#  SIESTA init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'siesta_init: finalize')
