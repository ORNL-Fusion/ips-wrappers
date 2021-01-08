#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC init component. This wapper only takes a VMEC input file
#  and runs VMEC.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from utilities import ZipState
from utilities import ScreenWriter
import os

#-------------------------------------------------------------------------------
#
#  VMEC init Component Constructor
#
#-------------------------------------------------------------------------------
class vmec_init(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  VMEC init Component init method. This method prepairs the state. Input files
#  can either be a new namelist input, a new state, or both.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'vmec_init: init')

#  Get config filenames.
        current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        
#  Remove old inputs. Stage input files.
        for file in os.listdir('.'):
            os.remove(file)
        
        self.services.stage_input_files(self.INPUT_FILES)
        
#  Create state from files. Input files can either be a new state, namelist
#  input file or both. If both file were staged, replace the namelist input
#  file. If the namelist file is present flag the state as needing to be
#  updated.
        with ZipState.ZipState(current_vmec_state, 'a') as zip_ref:
            if os.path.exists(current_vmec_namelist):
                zip_ref.write(current_vmec_namelist)
                zip_ref.set_state(state='needs_update')

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  VMEC init Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'vmec_init: step')

#-------------------------------------------------------------------------------
#
#  VMEC init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'vmec_init: finalize')
