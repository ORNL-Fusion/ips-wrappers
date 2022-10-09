#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for QUASI-NEWTON init component.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from ips_component_utilities import ZipState
from ips_component_utilities import ScreenWriter
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

#  Remove old inputs. Stage input files.
        for file in os.listdir('.'):
            os.remove(file)

        self.services.stage_input_files(self.INPUT_FILES)

#  Create state from files.
        with ZipState.ZipState(current_quasi_newton_state, 'a') as zip_ref:
            zip_ref.write(quasi_newton_config)
            zip_ref.write(current_model_state)

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Init step method.
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
