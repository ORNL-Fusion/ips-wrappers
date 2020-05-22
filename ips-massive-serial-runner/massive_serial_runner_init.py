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
#  Checks if the file needs to written, extracted, or rasies an exception if it
#  couldn't be found. This should be used when input files have been staged but
#  still being used by this component.
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
#  If the file can be written into the zip state file. If it cannot be check if
#  a file already exists. This should be used when input files have been staged
#  but not being used by this component.
#
#-------------------------------------------------------------------------------
def write_or_check(zip_ref, file):
    if os.path.exists(file):
        zip_ref.write(file):
    elif file not in zip_ref:
        raise Exception('Missing {} file.'.format(file))

#-------------------------------------------------------------------------------
#
#  Check if the file already exists. If it doesn't try to extract it from the
#  zip state file.
#
#-------------------------------------------------------------------------------
def extract_if_needed(zip_ref, file)
    if not os.path.exists(file) and file in zip_ref:
        zip_ref.extract(file)

#-------------------------------------------------------------------------------
#
#  Create an inscan file. The first line contains the headers starting with the
#  TIME_ID as a string type. The remaining headers are formated as the
#  prefix:name:float. To data is filled in by looping through the size of the
#  batch_size and pulling the correct index from the batch.
#
#-------------------------------------------------------------------------------
def create_inscan(current_batch, inscan_config):
    with open(current_batch, 'r') as json_ref:
        batch = json.load(json_ref)

    input = ':TIME_ID:str'

    for k, v in inscan_config.items():
        if k not in batch:
            raise Exception('Missing {} in {} file'.format(k, current_batch))
        input = '{} {}'.format(input, '{}:{}:float'.format(v[0], v[1]))

    batch_size = len(batch[batch.keys()[0]])
    for i in range(batch_size):
        input = '{}\n{:05d}'.format(input, i)

        for k in inscan_config:
            input = '{} {}'.format(input, batch[k][i])

    with open('inscan', 'w') as inscan_ref:
        inscan_ref.write(input)

    return batch_size

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
        if timeStamp = 0.0:
            self.current_state = self.services.get_config_param('CURRENT_MSR_STATE')
            msr_model_config = self.services.get_config_param('MSR_MODEL_CONFIG')
            self.current_batch = self.services.get_config_param('CURRENT_BATCH')
            database_config = self.services.get_config_param('DATABASE_CONFIG')
            self.inscan_config_file = self.services.get_config_param('INSCAN_CONFIG')

#  Remove old inputs.
        for file in os.listdir('.'):
            os.remove(file)

#  Stage input files and setup inital state.
        self.services.stage_input_files(self.INPUT_FILES)

#  Load or create a masive serial runner zip state.
        with ZipState.ZipState(current_state, 'a') as zip_ref:

#  Overwrite the msr_model_config and database_config file if they were staged
#  as input files. Over the write inscan_config if it was staged. otherwise
#  extract it. These files are not expected to change so we only need todo this
#  once.
            if timeStamp = 0.0:
                write_or_extract(zip_ref, self.inscan_config)
                write_or_check(zip_ref, database_config)
                write_or_check(zip_ref, msr_model_config)

                #  Load the inscan config file once.
                with open(self.inscan_config_file, 'r') as inscan_ref:
                    self.inscan_config = json.load(inscan_ref)

#  Batch files are optional. If a batch file was not staged as an input, extract
#  if from the plasma state if one exists inside it.
            extract_if_needed(self.current_batch)

#  Check if a new batch of data exists. If it does create the new inscan file
#  and write into the
            if os.path.exists(self.current_batch):
                zip_ref.set_state(batch_size=create_inscan(self.current_batch, self.inscan_config))
                zip_ref.write('inscan')

#  There maybe an existing inscan file in the state file. If one doesn't exist,
#  create a new one.
            elif 'inscan' not in zip_ref:
                task_wait = self.services.launch_task(self.NPROC,
                                                      self.services.get_working_dir(),
                                                      self.services.SAMPLE_EXE,
                                                      '--input={}'format(self.inscan_config_file),
                                                      '--output=inscan',
                                                      '--nscan=32',
                                                      logfile='sample_{}.log'.format(timeStamp))
                zip_ref.write('inscan')
                zip_ref.set_state(batch_size=32)

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
