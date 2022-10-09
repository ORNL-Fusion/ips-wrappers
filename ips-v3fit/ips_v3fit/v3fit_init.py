#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for V3FIT init component. This wapper only takes a V3FIT input
#  file and runs V3FIT.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from ips_component_utilities import ZipState
from ips_component_utilities import ScreenWriter
import os

#-------------------------------------------------------------------------------
#
#  V3FIT init Component Constructor
#
#-------------------------------------------------------------------------------
class v3fit_init(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  V3FIT init Component init method. This method prepairs the namelist input
#  file and creates a dummy out put file. This allows staging the state files.
#  In the v3fit namelist input file configure the v3fit namelist input with the
#  task, internal vmec input name and optional name of the wout file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'v3fit_init: init')

#  Get config filenames.
        current_model_config = self.services.get_config_param('MODEL_CONFIG')
        currnet_model_state = self.services.get_config_param('MODEL_STATE')

        current_v3fit_namelist = self.services.get_config_param('V3FIT_NAMELIST_INPUT')
        current_v3fit_state = self.services.get_config_param('CURRENT_V3FIT_STATE')

#  Remove old inputs. Stage input files.
        for file in os.listdir('.'):
            os.remove(file)
        
        self.services.stage_input_files(self.INPUT_FILES)

#  Create state from files. Input files can either be a new state, namelist
#  input file or both. If both files were staged, replace the namelist input
#  file. If the namelist file is present flag the state as needing to be
#  updated.
        with ZipState.ZipState(current_v3fit_state, 'a') as zip_ref:
            if os.path.exists(current_v3fit_namelist):
                zip_ref.write(current_v3fit_namelist)
                zip_ref.write(current_model_config)
                zip_ref.write(currnet_model_state)
                zip_ref.set_state(state='needs_update')

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  V3FIT init Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'v3fit_init: step')

#-------------------------------------------------------------------------------
#
#  V3FIT init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'v3fit_init: finalize')
