#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive Serial Runner component.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from ips_component_utilities import ZipState
from ips_component_utilities import ScreenWriter
import os
import json
import subprocess
import math

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component Constructor
#
#-------------------------------------------------------------------------------
class massive_serial_runner(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component init method. This method prepairs the state.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: init')

#  Get config filenames.
        if timeStamp == 0.0:
            self.current_state = self.services.get_config_param('CURRENT_MSR_STATE')
            self.database_config = self.services.get_config_param('DATABASE_CONFIG')
            self.current_batch = self.services.get_config_param('CURRENT_BATCH')

            self.constraint_path = self.services.get_config_param('MODULE_PATH')
            self.constraint_name = self.services.get_config_param('MODULE_NAME')

#  IPS framework config parameters.
            os.environ['PWD'] = os.getcwd()

            self.msr_config = self.services.get_config_param('MSR_CONFIG')
            self.msr_global_config = self.services.get_config_param('MSR_GLOBAL_CONFIG')
            self.msr_model_config = self.services.get_config_param('MSR_MODEL_CONFIG')
            self.msr_platform_conf = self.services.get_config_param('MSR_PLATFORM_FILE')
            self.msr_node_conf = self.services.get_config_param('MSR_NODE_FILE')
            os.environ['PLATFORM'] = self.msr_platform_conf
            os.environ['MSR_CONFIG'] = self.services.get_config_param('MSR_CONFIG')
            os.environ['MSR_MODEL_CONFIG'] = self.msr_model_config
            os.environ['NNODES'] = self.services.get_config_param('NNODES')

#  Stage state.
        self.services.stage_state()

#  Unzip files from the state. Use mode a so files an be read and written to.
        self.zip_ref = ZipState.ZipState(self.current_state, 'a')
        self.zip_ref.extract('inscan')

#  These files should never change so only extract them once.
        if timeStamp == 0.0:
            self.zip_ref.extract(self.database_config)
            self.zip_ref.extract(self.msr_config)
            self.zip_ref.extract(self.msr_global_config)
            self.zip_ref.extract(self.msr_model_config)
            self.zip_ref.extract(self.msr_platform_conf)
            self.zip_ref.extract(self.msr_node_conf)

#  We need the input directory to exist in a directory called input. So we must
#  make that directory first than extract the files. Remember to change back to
#  the orginal working directory after extraction.
            self.zip_ref.extract('input.zip')
            os.makedirs('input')
            os.chdir('input')
            with ZipState.ZipState('../input.zip', 'r') as input_ref:
                input_ref.extractall()
            os.chdir('../')

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component step method. Not used.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: step')

        flags = self.zip_ref.get_state()

#  Run the massive serial workflow.
        if 'state' in flags and flags['state'] == 'needs_update':
            process = subprocess.Popen([self.MASSIVE_SERIAL_EXE,
                                        '--platform={}'.format(self.msr_platform_conf),
                                        '--simulation={}'.format(self.msr_global_config),
                                        '--log=massive_serial_{}.log'.format(timeStamp)],
                                       env=os.environ)

            database = 'db_{}.dat'.format(timeStamp)

#  Collect results of the workflow into the database file.
            if process.wait():
                self.services.error('massive_serial_runner: step failed to run massive serial')

            print('here1')
            task_wait = self.services.launch_task(1, self.services.get_working_dir(),
                                                  self.MAKE_DATABASE_EXE,
                                                  '--rdir=output',
                                                  '--input={}'.format(self.database_config),
                                                  '--output={}'.format(database),
                                                  logfile='make_db_{}.log'.format(timeStamp))
            if self.services.wait_task(task_wait):
                self.services.error('massive_serial_runner: step failed to make database')

#  Save the new database entries.
            self.services.stage_output_files(timeStamp, database)

#  Convert the database file to json format.
            print('here2')
            task_wait = self.services.launch_task(1, self.services.get_working_dir(),
                                                  self.TO_JSON_EXE,
                                                  '--input_file={}'.format(database),
                                                  '--output_file={}'.format(self.current_batch),
                                                  '--module_path={}'.format(self.constraint_path),
                                                  '--module={}'.format(self.constraint_name),
                                                  logfile='to_json_{}.log'.format(timeStamp))

            if self.services.wait_task(task_wait):
                self.services.error('massive_serial_runner: step failed to make json')

            self.zip_ref.write(self.current_batch)
            self.zip_ref.set_state(state='updated')

        else:
            self.zip_ref.write(state='unchanged')

        self.zip_ref.close()

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component finalize method. This cleans up
#  afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: finalize')
