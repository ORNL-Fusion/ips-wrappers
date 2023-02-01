#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive Parallel Runner component.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from ips_component_utilities import ZipState
from ips_component_utilities import ScreenWriter
import os
import json
import subprocess
import math
import shutil
from configobj import ConfigObj

#-------------------------------------------------------------------------------
#
#  Massive Parallel Runner Component Constructor
#
#-------------------------------------------------------------------------------
class massive_parallel_runner(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  Massive Parallel Runner Component init method. This method prepairs the state.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'massive_parallel_runner: init')

#  Stage state.
        self.services.stage_state()

#  Unzip files from the state. Use mode a so files an be read and written to.
        if timeStamp == 0.0:
            self.current_state = self.services.get_config_param('CURRENT_MSR_STATE')
            self.save_logs = int(self.services.get_config_param('SAVE_LOGS'))

        self.zip_ref = ZipState.ZipState(self.current_state, 'a')

#  Get config filenames.
        if timeStamp == 0.0:
            self.database_config = self.services.get_config_param('DATABASE_CONFIG')
            self.zip_ref.extract(self.database_config)

            self.current_batch = self.services.get_config_param('CURRENT_BATCH')

            self.constraint_path = self.services.get_config_param('MODULE_PATH')
            self.constraint_name = self.services.get_config_param('MODULE_NAME')

#  Keys for the massiver parallel subworkflow.
            keys = {
                'PWD'            : self.services.get_config_param('PWD'),
                'SIM_NAME'       : 'massive_parallel_runner_sub',
                'LOG_FILE'       : 'log.massive_parallel_runner',
                'NNODES'         : self.services.get_config_param('MSR_NNODES'),
                'INPUT_DIR_SIM'  : 'massive_parallel_runner_input_dir',
                'OUTPUT_DIR_SIM' : '{}/massive_parallel_runner_output_dir'.format(os.getcwd())
            }

            if os.path.exists('massive_parallel_runner_input_dir'):
                shutil.rmtree('massive_parallel_runner_input_dir')
            os.mkdir('massive_parallel_runner_input_dir')

            self.massive_parallel_worker = {
                'sim_name' : None,
                'init'     : None,
                'driver'   : None
            }

            msr_global = self.services.get_config_param('MSR_GLOBAL_CONFIG')
            self.zip_ref.extract(msr_global)

            (self.massive_parallel_worker['sim_name'],
             self.massive_parallel_worker['init'],
             self.massive_parallel_worker['driver']) = self.services.create_sub_workflow('massive_parallel',
                                                                                         msr_global,
                                                                                         keys,
                                                                                         'massive_parallel_runner_input_dir')

        self.zip_ref.extract('inscan')
        shutil.copy2('inscan', 'massive_parallel_runner_input_dir')

#  These files should never change so only extract them once.
        if timeStamp == 0.0:
            ms_state = self.services.get_config_param('MSR_SERIAL_STATE')

            self.zip_ref.extract(ms_state)
            shutil.copy2(ms_state, 'massive_parallel_runner_input_dir')

            os.chdir('massive_parallel_runner_input_dir')

#  We need the input directory to exist in a directory called input. So we must
#  make that directory first than extract the files. Remember to change back to
#  the orginal working directory after extraction.
            with ZipState.ZipState(ms_state, 'r') as zip_ref:
                zip_ref.extractall()
            with ZipState.ZipState('input.zip', 'r') as input_ref:
                input_ref.extractall()

            override = ConfigObj(infile=self.services.get_config_param('MSR_SERIAL_NODE_CONFIG'), interpolation='template', file_error=True)
            override['INPUT_DIR_SIM'] = os.getcwd()
            override.write()

            override2 = ConfigObj(infile=self.services.get_config_param('MSR_MODEL_CONFIG'), interpolation='template', file_error=True)
            override2['INPUT_DIR_SIM'] = os.getcwd()
            override2.write()

            os.chdir('../')

#-------------------------------------------------------------------------------
#
#  Massive Parallel Runner Component step method. Not used.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_parallel_runner: step')

        flags = self.zip_ref.get_state()

#  Run the massive parallel workflow.
        if 'state' in flags and flags['state'] == 'needs_update':
            self.services.call(self.massive_parallel_worker['init'], 'init', timeStamp)
            self.services.call(self.massive_parallel_worker['driver'], 'init', timeStamp)
            wait = self.services.call_nonblocking(self.massive_parallel_worker['driver'], 'step', timeStamp)

            database = 'db_{}.dat'.format(timeStamp)

#  Collect results of the workflow into the database file.
            if self.services.wait_call(wait, True):
                self.services.error('massive_parallel_runner: step failed to run massive parallel')

            if self.save_logs:
                logfile='make_db_{}.log'.format(timeStamp)
            else:
                logfile = 'make_db.log'

            task_wait = self.services.launch_task(1, self.services.get_working_dir(),
                                                  self.MAKE_DATABASE_EXE,
                                                  '--rdir=massive_parallel_runner_output_dir',
                                                  '--input={}'.format(self.database_config),
                                                  '--output={}'.format(database),
                                                  '--ndir=0', #  FIXME: This command option works around a bug in makedb which shouldn't get called.
                                                  logfile=logfile)

            if self.services.wait_task(task_wait):
                self.services.error('massive_parallel_runner: step failed to make database')

#  Save the new database entries.
            self.services.stage_output_files(timeStamp, database)

#  Convert the database file to json format.
            if self.save_logs:
                logfile='to_json_{}.log'.format(timeStamp)
            else:
                logfile = 'to_json.log'

            task_wait = self.services.launch_task(1, self.services.get_working_dir(),
                                                  self.TO_JSON_EXE,
                                                  '--input_file={}'.format(database),
                                                  '--output_file={}'.format(self.current_batch),
                                                  '--module_path={}'.format(self.constraint_path),
                                                  '--module={}'.format(self.constraint_name),
                                                  logfile=logfile)

            if self.services.wait_task(task_wait):
                self.services.error('massive_parallel_runner: step failed to make json')

            self.zip_ref.write(self.current_batch)
            self.zip_ref.set_state(state='updated')

        else:
            self.zip_ref.write(state='unchanged')

        self.zip_ref.close()

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  Massive Parallel Runner Component finalize method. This cleans up
#  afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_parallel_runner: finalize')
