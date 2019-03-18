#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for ML Train component. This driver only runs the ML Train
#  component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
import os

#-------------------------------------------------------------------------------
#
#  ML Train Driver Constructor
#
#-------------------------------------------------------------------------------
class ml_train_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  ML Train Driver init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train_driver: init')

#  Initialize ml_train.
        self.ml_train_port = self.services.get_port('ML_TRAIN')
        self.wait = self.services.call_nonblocking(self.ml_train_port, 'init',
                                                   timeStamp, **keywords)

#-------------------------------------------------------------------------------
#
#  ML Train Driver step method. This runs the ml_train component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train_driver: step')

#  Run ml_train.
        self.services.wait_call(self.wait, True)
        self.services.call(self.ml_train_port, 'step', timeStamp)

#  Prepare the output files for a super work flow. Need to remove any old output
#  files first before staging the plasma state.
        if os.path.exists(self.OUTPUT_FILES):
            os.remove(self.OUTPUT_FILES)
        self.services.stage_plasma_state()

#  The super flow may need to rename the output file. Check is the current state
#  matches if output file. If it does not rename the plasma state so it can be
#  staged.
        if not os.path.exists(self.OUTPUT_FILES):
            os.rename(self.services.get_config_param('CURRENT_ML_TRAIN_STATE'),
                      self.OUTPUT_FILES)

#-------------------------------------------------------------------------------
#
#  ML Train Driver finalize method. This cleans up afterwards.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train_driver: finalize')
        self.services.call(self.ml_train_port, 'finalize', timeStamp)
