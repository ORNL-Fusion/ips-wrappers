#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive Serial Runner init component.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from ips_component_utilities import ZipState
from ips_component_utilities import ScreenWriter
import os
import json

#-------------------------------------------------------------------------------
#
#  Check if the file already exists. If it doesn't try to extract it from the
#  zip state file.
#
#-------------------------------------------------------------------------------
def extract_if_needed(zip_ref, file):
    if not os.path.exists(file) and file in zip_ref:
        zip_ref.extract(file)

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
        if timeStamp == 0.0:
            self.current_state = self.services.get_config_param('CURRENT_MSR_STATE')
            self.msr_config = self.services.get_config_param('MSR_CONFIG')
            self.msr_global_config = self.services.get_config_param('MSR_GLOBAL_CONFIG')
            self.msr_model_config = self.services.get_config_param('MSR_MODEL_CONFIG')
            self.msr_platform_conf = self.services.get_config_param('MSR_PLATFORM_FILE')
            self.current_batch = self.services.get_config_param('CURRENT_BATCH')
            self.database_config = self.services.get_config_param('DATABASE_CONFIG')
            self.inscan_config_file = self.services.get_config_param('INSCAN_CONFIG')
            self.batch_size = self.services.get_config_param('BATCH_SIZE')

#  Remove old inputs.
        for file in os.listdir('.'):
            os.remove(file)

#  Stage input files and setup inital state.
        self.services.stage_input_files(self.INPUT_FILES)

#  Load or create a masive serial runner zip state.
        with ZipState.ZipState(self.current_state, 'a') as zip_ref:

#  Overwrite the msr_model_config and database_config file if they were staged
#  as input files. Over the write inscan_config if it was staged. otherwise
#  extract it. These files are not expected to change so we only need todo this
#  once.
            zip_ref.write_or_extract(self.inscan_config_file)
            zip_ref.write_or_check(self.database_config)
            zip_ref.write_or_check(self.msr_config)
            zip_ref.write_or_check(self.msr_global_config)
            zip_ref.write_or_check(self.msr_model_config)
            zip_ref.write_or_check(self.msr_platform_conf)

#  This is the inputs to fastran. Should only need todo this once.
            if timeStamp == 0.0:
                zip_ref.write('input.zip')

#  Batch files are optional. If a batch file was not staged as an input, extract
#  if from the plasma state if one exists inside it.
            extract_if_needed(zip_ref, self.current_batch)

#  Check if a new batch of data exists. If it does create the new inscan file
#  and write into the
            if os.path.exists(self.current_batch):
                task_wait = self.services.launch_task(self.NPROC,
                                                      self.services.get_working_dir(),
                                                      self.SAMPLE_EXE,
                                                      '--input={}'.format(self.inscan_config_file),
                                                      '--output=inscan',
                                                      '--nscan={}'.format(self.batch_size),
                                                      '--new={}'.format(self.current_batch),
                                                      logfile='sample_{}.log'.format(timeStamp))

                if self.services.wait_task(task_wait):
                    self.services.error('massive_serial_runner_init: failed to generate inscan sample')

#  There maybe an existing inscan file in the state file. If one doesn't exist,
#  create a new one.
            elif 'inscan' not in zip_ref:
                task_wait = self.services.launch_task(self.NPROC,
                                                      self.services.get_working_dir(),
                                                      self.SAMPLE_EXE,
                                                      '--input={}'.format(self.inscan_config_file),
                                                      '--output=inscan',
                                                      '--nscan={}'.format(self.batch_size),
                                                      logfile='sample_{}.log'.format(timeStamp))

                if self.services.wait_task(task_wait):
                    self.services.error('massive_serial_runner_init: failed to generate inscan sample')

            else:
                raise Expection('Expected inscan or {} file not found.'.format(self.current_batch))

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
