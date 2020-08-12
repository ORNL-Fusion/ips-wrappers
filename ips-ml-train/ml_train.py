#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for ML Train component. This wapper runs the actual training.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
import json
import os

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

        if timeStamp == 0.0:
            self.current_ml_train_state = self.services.get_config_param('CURRENT_ML_TRAIN_STATE')
            self.training_data = self.services.get_config_param('TRAINING_DATA')
            self.new_data = self.services.get_config_param('NEW_DATA')
            self.prediction_data = self.services.get_config_param('PREDICTION_DATA')

            self.nn_model_config = self.services.get_config_param('NN_MODEL_CONFIG')
            self.nn_model_matrix = self.services.get_config_param('NN_MODEL_MATRIX')
            self.nn_model = self.services.get_config_param('NN_MODEL')
            self.batch_size = self.services.get_config_param('BATCH_SIZE')

            self.constraint_path = self.services.get_config_param('MODULE_PATH')
            self.constraint_name = self.services.get_config_param('MODULE_NAME')

#  Stage state.
        self.services.stage_state()

#  Unzip files from the current state. Use mode a so files can be read and
#  written to.
        self.zip_ref = ZipState.ZipState(self.current_ml_train_state, 'a')

        if timeStamp == 0.0:
            ml_train_args = self.services.get_config_param('ML_TRAIN_ARGS')
            self.zip_ref.extract_or_check(ml_train_args)
            self.zip_ref.extract_or_check(self.nn_model_config)

            with open(ml_train_args, 'r') as args_ref:
                self.args = json.load(args_ref)


        self.zip_ref.extract_or_check(self.training_data)
        self.zip_ref.extract_optional('{}.zip'.format(self.nn_model))

        if os.path.exists('{}.zip'.format(self.nn_model)):
            with ZipState.ZipState('{}.zip'.format(self.nn_model), 'r') as nn_ref:
                nn_ref.extractall()

#-------------------------------------------------------------------------------
#
#  ML Train step method. This runs the ml_train component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train: step')

        flags = self.zip_ref.get_state()

        if 'state' in flags and flags['state'] == 'needs_update':
            if os.path.exists(self.nn_model):
                task_wait = self.services.launch_task(self.NPROC,
                                                      self.services.get_working_dir(),
                                                      self.ML_TRAIN_EXE,
                                                      '--config={}'.format(self.nn_model_config),
                                                      '--model={}'.format(self.nn_model),
                                                      '--training_data={}'.format(self.training_data),
                                                      '--supplemental_data={}'.format(self.new_data),
                                                      '--prediction_data={}'.format(self.prediction_data),
                                                      '--batch_size={}'.format(self.batch_size),
                                                      '--iterations={}'.format(self.args['--iterations']),
                                                      '--epochs={}'.format(self.args['--epochs']),
                                                      '--param_covar_matrix={}'.format(self.nn_model_matrix),
                                                      '--validation_split={}'.format(self.args['--validation_split']),
                                                      '--module_path={}'.format(self.constraint_path),
                                                      '--module={}'.format(self.constraint_name),
                                                      logfile = 'ml_train_{}.log'.format(timeStamp))
            else:
                task_wait = self.services.launch_task(self.NPROC,
                                                      self.services.get_working_dir(),
                                                      self.ML_TRAIN_EXE,
                                                      '--config={}'.format(self.nn_model_config),
                                                      '--model={}'.format(self.nn_model),
                                                      '--activation={}'.format(self.args['--activation']),
                                                      '--training_data={}'.format(self.training_data),
                                                      '--supplemental_data={}'.format(self.new_data),
                                                      '--prediction_data={}'.format(self.prediction_data),
                                                      '--batch_size={}'.format(self.batch_size),
                                                      '--iterations={}'.format(self.args['--iterations']),
                                                      '--epochs={}'.format(self.args['--epochs']),
                                                      '--num_layers={}'.format(self.args['--num_layers']),
                                                      '--layer_width={}'.format(self.args['--layer_width']),
                                                      '--param_covar_matrix={}'.format(self.nn_model_matrix),
                                                      '--l1_factor={}'.format(self.args['--l1_factor']),
                                                      '--l2_factor={}'.format(self.args['--l2_factor']),
                                                      '--validation_split={}'.format(self.args['--validation_split']),
                                                      '--module_path={}'.format(self.constraint_path),
                                                      '--module={}'.format(self.constraint_name),
                                                      logfile = 'ml_train_{}.log'.format(timeStamp))



#  Update flags.
            self.zip_ref.set_state(state='updated')

#  Wait for Training to finish. FIXME: Need to check that the outputs exist to
#  check errors
            if (self.services.wait_task(task_wait) and False):
                self.services.error('ml_train: step failed.')

#  NN models may be a directory. Zip them first before adding them to the state.
            with ZipState.ZipState('{}.zip'.format(self.nn_model), 'w') as nn_ref:
                nn_ref.write(self.nn_model)
            self.zip_ref.write('{}.zip'.format(self.nn_model))

            self.zip_ref.write(self.new_data)
            self.zip_ref.write(self.prediction_data)
            self.zip_ref.write(self.nn_model_matrix)

#  Add outputs to state.
        else:
#  Update flags.
            self.zip_ref.set_state(state='unchanged')

        self.zip_ref.close()

        self.services.update_state()
        self.services.stage_output_files(timeStamp, self.current_ml_train_state)

#-------------------------------------------------------------------------------
#
#  ML Train finalize method. This cleans up afterwards.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train: finalize')
