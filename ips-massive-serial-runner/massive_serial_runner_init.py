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
import json

#-------------------------------------------------------------------------------
#
#  Checks if the file needs to read, extracted, or rasies an exception if it
#  couldn't be found.
#
#-------------------------------------------------------------------------------
def write_or_extract(zip_ref, file):
    if os.path.exists(file):
        zip_ref.write(file)
    elif file in zip_ref:
        zip_ref.extract(file)
    else:
        raise Exception('Missing {} file.'.format(file))

#-------------------------------------------------------------------------------
#
#  Create an inscan file. The first line contains the headers starting with the
#  TIME_ID as a string type. The remaining headers are formated as the
#  prefix:name:float. To data is filled in by looping through the size of the
#  batch_size and pulling the correct index from the batch.
#
#-------------------------------------------------------------------------------
def create_inscan(batch, inscan_config):
    input = ':TIME_ID:str'

    for k, v in inscan_config.items():
        input = '{} {}'.format(input, '{}:{}:float'.format(v[0], v[1]))

    batch_size = len(batch[batch.keys()[0]])
    for i in range(batch_size):
        input = '{}\n{:05d}'.format(input, i)

        for k in inscan_config:
            input = '{} {}'.format(input, batch[k][i])

    with open('inscan', 'w') as inscan_ref:
        inscan_ref.write(input)

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
#  state and generates an inscan file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_init: init')

#  Get config filenames.
        current_batch = self.services.get_config_param('CURRENT_MSR_BATCH')
        current_state = self.services.get_config_param('CURRENT_MSR_STATE')
        inscan_config = self.services.get_config_param('INSCAN_CONFIG')
        model_config = self.services.get_config_param('MODEL_CONFIG')
        database_config = self.services.get_config_param('DATABASE_CONFIG')

#  Remove old inputs.
        for file in os.listdir('.'):
            os.remove(file)

#  Stage input files and setup inital state.
        self.services.stage_input_files(self.INPUT_FILES)

#  Load or create a masive serial runner zip state. And over write the model
#  current batch, inscan config and model_config is they were staged as input
#  files. Otherwise extract them from the zipstate.
        with ZipState.ZipState(current_state, 'a') as zip_ref:
            write_or_extract(zip_ref, inscan_config)
            write_or_extract(zip_ref, model_config)
            write_or_extract(zip_ref, database_config)

#  Load the inscan config file and model config file.
            with open(current_inscan_config, 'r') as inscan_ref:
                inscan_config = json.load(inscan_ref)
            with open(model_config, 'r') as model_ref:
                model_config = json.load(model_ref)

#  Check that model inputs are valid and not constant.
            for value in model_config['inputs']:
                if value['name'] not in inscan_config['io']:
                    raise Exception('Model input {} is not a valid parameter.'.format(value['name']))
                if value['name'] in inscan_config['const']:
                    raise Exception('Model input {} is not a valid scanable.'.format(value['name']))

#  Check that the model inputs match the inscan inputs.
            for key in inscan_config['scan'].keys()
                if key not in model_config['inputs']:
                    raise Exception('Model input {} is not a valid')

#  Check if a new batch of data exists. If it does create the new inscan file
#  and write into the
            if os.path.exists(current_batch):
                with open(current_batch, 'r') as json_ref:
                    create_inscan(json.load(json_ref), inscan_config)
                zip_ref.write('inscan')
            elif 'inscan' not in zip_ref:
                task_wait = self.services.launch_task(self.NPROC,
                                                      self.services.get_working_dir(),
                                                      self.services.SAMPLE_EXE,
                                                      '--input={}'format(inscan_config),
                                                      '--output=inscan',
                                                      '--nscan=32')
                zip_ref.write('inscan')

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
