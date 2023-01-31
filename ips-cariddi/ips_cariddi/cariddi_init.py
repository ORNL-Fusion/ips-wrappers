#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for CARIDDI init component. This wapper only takes a CARIDDI
#  inputs and creates the inital state file.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from ips_component_utilities import ZipState
from ips_component_utilities import ScreenWriter
import os

#-------------------------------------------------------------------------------
#
#  CARIDDI init Component Constructor
#
#-------------------------------------------------------------------------------
class cariddi_init(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  CARIDDI init Component init method. This method prepairs the state. Input
#  files can either be a new input, a new state, or both.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi_init: init')

#  Get config filenames.
        current_cariddi_state = self.services.get_config_param('CURRENT_CARIDDI_STATE')
        current_v3fit_config = self.services.get_config_param('V3FIT_CONFIG')
        current_v3fit_state = self.services.get_config_param('CURRENT_V3FIT_STATE')

#  Remove old inputs. Stage input files.
        for file in os.listdir('.'):
            os.remove(file)

        self.services.stage_input_files(self.INPUT_FILES)

#  Create state from files. Input files can either be a new state, input file or
#  both. If both files were staged, replace the input file. If the input file is
#  present flag the state as needing to be updated.
        with ZipState.ZipState(current_cariddi_state, 'a') as zip_ref:
            zip_ref.write(current_v3fit_state)
            zip_ref.write(current_v3fit_config)

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  CARIDDI init Component step method. Not used.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi_init: step')

#-------------------------------------------------------------------------------
#
#  CARIDDI init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi_init: finalize')
