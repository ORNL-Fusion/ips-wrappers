#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for ML Gen Data component. This wapper launches a massive serial
#  job to generate more training data.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
import os
import numpy
import random
import subprocess

#-------------------------------------------------------------------------------
#
#  ML Gen Data Constructor
#
#-------------------------------------------------------------------------------
class ml_gen_data(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  ML Gen Data Init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'ml_gen_data: init')

#  Perform the inital setup.
        current_ml_model_state = ''
        if timeStamp == 0.0:
            os.environ['PWD'] = os.getcwd()
            self.current_ml_train_state = self.services.get_config_param('CURRENT_ML_TRAIN_STATE')
            self.current_ml_train_data = self.services.get_config_param('CURRENT_ML_TRAIN_DATA')
            self.current_ml_train_new_data = self.services.get_config_param('CURRENT_ML_TRAIN_NEW_DATA')

            current_ml_model_state = self.services.get_config_param('CURRENT_ML_MODEL_STATE')

            self.model_sim_config = self.services.get_config_param('MODEL_SIM_CONFIG')

            if not os.path.exists('inputs')
                os.makedirs('inputs')

#  Stage state.
        self.services.stage_state()

#  Unzip files from the current state. Use mode a so files can be read and
#  written to.
        self.zip_ref = ZipState.ZipState(self.current_ml_train_state, 'a')
        
        if self.current_ml_train_new_data in self.zip_ref:
            self.zip_ref.set_state(state='needs_update')

            self.zip_ref.extract(self.current_ml_train_new_data)

            with open(self.current_ml_train_new_data, 'r') as new_data_ref:
                new_data = json.load(new_data_ref)

                with open('massive_input') as massive_ref:
                    for param in keywords['input_params']:
                        massive_ref.write('{prefix}:{name}:{type} '.format(param))
                    massive_ref.write('\n')

                    permutation_size = min(keywords['batch_size'],
                                           len(new_data[keywords['input_params'][0]['name']]))
                    permutation = numpy.random.permutation(permutation_size)
                    for i in permutation:
                        for param in keywords['input_params']:
                            massive_ref.write('{} '.format(new_data[param['name']][i]))
                        massive_ref.write('\n')

        elif self.current_ml_train_data in self.zip_ref:
            self.zip_ref.set_state(state='updated')

        else:
            self.zip_ref.set_state(state='needs_update')

            with open('massive_input') as massive_ref:
                for param in keywords['input_params']:
                    massive_ref.write('{prefix}:{name}:{type} '.format(param))
                massive_ref.write('\n')

                for i in batch_size:
                    for param in keywords['input_params']:
                        massive_ref.write('{} '.format(random.uniform(param['lower_range'], param['upper_range'])))
                    massive_ref.write('\n')

        self.zip_ref.extract(current_ml_model_state)
        if keywords['extract'] and timeStamp == 0.0:
            with ZipState.ZipState(current_ml_model_state, 'a') as model_ref:
                model_ref.extractall('inputs')
        elif timeStamp == 0.0:
            os.rename(current_ml_model_state, 'inputs/{}'.format(current_ml_model_state))

#-------------------------------------------------------------------------------
#
#  ML Gen Data Step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self. 'verbose', 'ml_gen_data: step')

        flags = self.zip_ref.get_state()

        if 'state' in flags and flags['state'] == 'needs_update':
            process = subprocess.Popen([, massive_input, self.model_sim_config, , 1])
            process.wait();
        
        

        else:
#  Update flags.
            self.zip_ref.set_state(state='unchanged')

        zip_ref.close()

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  ML Gen Data Finalize method.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'ml_gen_data: finalize')
