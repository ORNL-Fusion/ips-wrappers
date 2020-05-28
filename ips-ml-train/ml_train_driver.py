#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for ML Train component. This driver only runs the ML Train
#  component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
import os
import shutil

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

        if timeStamp == 0.0:
            self.current_ml_train_state = self.services.get_config_param('CURRENT_ML_TRAIN_STATE')
            self.new_data = self.services.get_config_param('NEW_DATA')
            self.training_data = self.services.get_config_param('TRAINING_DATA')
            self.data_gen_state = self.services.get_config_param('DATA_GEN_STATE')

            self.ml_train_port = self.services.get_port('ML_TRAIN')

        self.services.stage_state()

#  Extract files needed to set up the data gen model.
        with ZipState.ZipState(self.current_ml_train_state, 'a') as zip_ref:
            zip_ref.extract_optional(self.new_data)
            zip_ref.extract_or_check(self.data_gen_state)

            with ZipState.ZipState(self.data_gen_state, 'a') as model_state_ref:
                model_state_ref.write_optional(self.new_data)

#  Setup the data generation subworkflow.
            if timeStamp == 0.0:
                data_gen_config = self.services.get_config_param('DATA_GEN_CONFIG')
                zip_ref.extract_or_check(data_gen_config)

#  Get keys for the sub workflow.
                keys = {
                    'pwd'           : self.services.get_config_param('PWD'),
                    'SIM_NAME'      : '{}_gen_data'.format(self.services.get_config_param('SIM_NAME')),
                    'LOG_FILE'      : 'log.gen_data.warning',
                    'OUTPUT_LEVEL'  : self.services.get_config_param('OUTPUT_LEVEL'),
                    'CURRENT_BATCH' : self.new_data
                }

                if os.path.exists('data_gen_input_dir'):
                    shutil.rmtree('data_gen_input_dir')
                os.mkdir('data_gen_input_dir')

                self.data_gen = {
                    'sim_name' : None,
                    'init'     : None,
                    'driver'   : None
                }
                (self.data_gen['sim_name'],
                 self.data_gen['init'],
                 self.data_gen['driver']) = self.services.create_sub_workflow('data_gen',
                                                                              data_gen_config,
                                                                              keys,
                                                                              'data_gen_input_dir')

            shutil.copy2(self.data_gen_state, 'data_gen_input_dir')

#  If new data exists or a training data does not exist. Generate a data batch.
            if os.path.exists(self.new_data) or self.training_data not in zip_ref:
                self.services.call(self.data_gen['init'], 'init', timeStamp)
                self.services.call(self.data_gen['driver'], 'init', timeStamp)
                self.services.call(self.data_gen['driver'], 'driver', timeStamp)

                self.services.stage_subflow_output_files()

                with ZipState.ZipState(self.data_gen_state, 'r') as model_state_ref:
                    model_state_ref.extract_or_check(self.new_data)

                if self.training_data not in zip_ref:
                    os.rename(self.new_data, self.training_data)
                else:
                    zip_ref.extract(self.training_data)
                    self.append_data()

                zip_ref.write(self.training_data)
                zip_ref.write(self.data_gen_state)
                zip_ref.set_state(state='needs_update')

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  ML Train Driver step method. This runs the ml_train component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train_driver: step')

        with ZipState.ZipState(self.current_ml_train_state, 'r') as zip_ref:
            flags = zip_ref.get_state()

#  Adaptive training loop. Hard code the number of iterations for now.
        for i in range(10):

#  Train the NN model.
            self.services.call(self.ml_train_port, 'init', timeStamp)
            self.services.call(self.ml_train_port, 'step', timeStamp)

#  Get the new data batch and add it to the gen data state.
            self.services.stage_state()
            with ZipState.ZipState(self.current_ml_train_state, 'a') as zip_ref:
                zip_ref.extract_or_check(self.new_data)
                zip_ref.extract_or_check(self.training_data)

                with ZipState.ZipState(self.data_gen_state, 'a') as model_state_ref:
                    model_state_ref.write(self.new_data)

#  Generate new training data.
                shutil.copy2(self.data_gen_state, 'data_gen_input_dir')
                self.services.call(self.data_gen['init'], 'init', timeStamp)
                self.services.call(self.data_gen['driver'], 'init', timeStamp)
                self.services.call(self.data_gen['driver'], 'driver', timeStamp)

                self.services.stage_subflow_output_files()

#  Append the new data to he training data.
                with ZipState.ZipState(self.data_gen_state, 'r') as model_state_ref:
                    model_state_ref.extract_or_check(self.new_data)
                self.append_data()

                zip_ref.set_state(state='needs_update')
                flags = zip_ref.get_state()

            self.services.update_state()
            timeStamp = timeStamp + 1.0

#-------------------------------------------------------------------------------
#
#  ML Train Driver finalize method. This cleans up afterwards.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_train_driver: finalize')
        self.services.call(self.ml_train_port, 'finalize', timeStamp)

#-------------------------------------------------------------------------------
#
#  ML Train Driver append data. Appends new data to the training data.
#
#-------------------------------------------------------------------------------
    def append_data(self)
        with open(self.training_data, 'a') as training_data_ref:
            train = json.load(training_data_ref)

            with open(self.new_data, 'r') as new_data_ref:
                new = json.load(new_data_ref)

            for k, v in train.items():
                train[k] = v + new[k]

            json.dump(train, training_data_ref)
