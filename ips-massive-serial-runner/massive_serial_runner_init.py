#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive Serial Runner init component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
import os

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component Constructor
#
#-------------------------------------------------------------------------------
class massive_serial_runner_init(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component init method. This method prepairs the
#  state.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_init: init')

#  Get config filenames.
        current_batch = self.services.get_config_param('CURRENT_MSR_BATCH')
        current_state = self.services.get_config_param('CURRENT_MSR_STATE')

#  Remove old inputs.
        for file in os.listdir('.'):
            os.remove(file)

#  Stage input files and setup inital state.
        self.services.stage_input_files(self.INPUT_FILES)

#  Create plasma state from files. Input files can either be a new plasma state,
#  training data file or both. If both file were staged, replace the training
#  data input file. If the training data file is present flag the plasma state
#  as needing to be updated.
        with ZipState.ZipState(current_state, 'a') as zip_ref:
            if os.path.exists(current_batch):
                zip_ref.write(current_batch)
                zip_ref.set_state(state='needs_update')

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component step method. Not used.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_init: step')

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component finalize method. This cleans up
#  afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_init: finalize')
