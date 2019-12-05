#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for ML Train component. This wapper runs the actual training.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState

#-------------------------------------------------------------------------------
#
#  ML Train Constructor
#
#-------------------------------------------------------------------------------
class ml_train(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  ML Train init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train: init')

        current_ml_train_state = self.services.get_config_param('CURRENT_ML_TRAIN_STATE')

#  Stage state.
        self.services.stage_state()

#  Unzip files from the current state. Use mode a so files can be read and
#  written to.
        self.zip_ref = ZipState.ZipState(current_ml_train_state, 'a')
        self.zip_ref.extract('training_data.json')

#-------------------------------------------------------------------------------
#
#  ML Train step method. This runs the ml_train component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train: step')

        flags = self.zip_ref.get_state()

        if 'state' in flags and flags['state'] == 'needs_update':
            task_wait = self.services.launch_task(self.NPROC,
                                                  self.services.get_working_dir(),
                                                  'python',
                                                  self.ML_TRAIN_EXE,
                                                  logfile = 'ml_train.log')

#  Update flags.
            self.zip_ref.set_state(state='updated')

#  Wait for Training to finish. FIXME: Need to check that the outputs exist to
#  check errors
            if (self.services.wait_task(task_wait) and False):
                self.services.error('ml_train: step failed.')

#  Add outputs to state.
        else:
#  Update flags.
            self.zip_ref.set_state(state='unchanged')

        self.zip_ref.close()

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  ML Train finalize method. This cleans up afterwards.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train: finalize')
