#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for ML Train init component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
import os

#-------------------------------------------------------------------------------
#
#  ML Train init Component Constructor
#
#-------------------------------------------------------------------------------
class ml_train_init(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  ML Train init Component init method. This method prepairs the plasma state.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train_init: init')

#  Get config filenames.
        if timeStamp == 0.0:
            self.current_ml_train_state = self.services.get_config_param('CURRENT_ML_TRAIN_STATE')

            self.data_gen_config = self.services.get_config_param('DATA_GEN_CONFIG')
            self.data_gen_state = self.services.get_config_param('DATA_GEN_STATE')

            self.training_data = self.services.get_config_param('TRAINING_DATA')
            self.new_data = self.services.get_config_param('NEW_DATA')
            self.prediction_data = self.services.get_config_param('PREDICTION_DATA')

            self.nn_model_config = self.services.get_config_param('NN_MODEL_CONFIG')
            self.nn_model_matrix = self.services.get_config_param('NN_MODEL_MATRIX')
            self.nn_model = self.services.get_config_param('NN_MODEL')
            self.ml_train_args = self.services.get_config_param('ML_TRAIN_ARGS')

#  Remove old inputs. Stage input files.
        for file in os.listdir('.'):
            os.remove(file)

#  State input files and setup the inital state.
        self.services.stage_input_files(self.INPUT_FILES)

#  Create plasma state from files. Input files can either be a new plasma state,
#  training data file or both. If both file were staged, replace the training
#  data input file. If the training data file is present flag the plasma state
#  as needing to be updated.
        with ZipState.ZipState(self.current_ml_train_state, 'a') as zip_ref:
            zip_ref.write_or_check(self.data_gen_config)
            zip_ref.write_or_check(self.data_gen_state)

            zip_ref.write_optional(self.training_data)
            zip_ref.write_optional(self.new_data)
            zip_ref.write_optional(self.prediction_data)

            zip_ref.write_or_check(self.nn_model_config)
            zip_ref.write_optional(self.nn_model_matrix)
            zip_ref.write_optional(self.nn_model)
            zip_ref.write_or_check(self.nn_model_config)

            zip_ref.set_state(state='needs_update')

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  ML Train init Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train_init: step')
    
#-------------------------------------------------------------------------------
#
#  ML Train init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train_init: finalize')
