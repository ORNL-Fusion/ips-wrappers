#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for QUASI-NEWTON init component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ZipState
from utilities import ScreenWriter
import os

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Init Constructor
#
#-------------------------------------------------------------------------------
class quasi_newton_init(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Init init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_init: init')
    
#  Get config filenames.
        current_model_state = self.services.get_config_param('MODEL_INPUT')
        quasi_newton_config = self.services.get_config_param('QUASI_NEWTON_CONFIG')
        current_quasi_newton_state = self.services.get_config_param('CURRENT_QUASI_NEWTON_STATE')

#  State input files. Remove old files if they exist.
        if os.path.exists(current_model_state):
            os.remove(current_model_state)
        if os.path.exists(quasi_newton_config):
            os.remove(quasi_newton_config)
        if os.path.exists(current_quasi_newton_state):
            os.remove(current_quasi_newton_state)
        self.services.stage_input_files(self.INPUT_FILES)

#  Create plasma state from files.
        with ZipState.ZipState(current_quasi_newton_state, 'a') as zip_ref:
            zip_ref.write(quasi_newton_config)
            zip_ref.write(current_model_state)

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Init step method. This runs the vmec component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_init: step')

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Init finalize method. This cleans up afterwards.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_init: finalize')
