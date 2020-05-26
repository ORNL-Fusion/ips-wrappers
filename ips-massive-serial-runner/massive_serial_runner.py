#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive Serial Runner component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
import os
import json

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
            self.msr_config = self.services.get_config_param('MSR_MODEL_CONFIG')

#  Stage state.
        self.services.stage_state()

#  Unzip files from the state. Use mode a so files an be read and written to.
        self.zip_ref = ZipState.ZipState(self.current_state, 'a')
        self.zip_ref.extract('inscan')

#  These file should never change so only extract them once.
        if timeStamp == 0.0:
            self.zip_ref.extract(self.database_config)
            self.zip_ref.extract(self.msr_config)

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
            task_wait = self.services.launch_task(self.NPROC,
                                                  self.services.get_working_dir(),
                                                  self.MASSIVE_SERIAL_EXE,
                                                  'inscan',
                                                  self.msr_config,
                                                  '{}'.format(max(32,
                                                              int(flags['batch_size']))),
                                                  '0',
                                                  logfile='massive_serial_{}.log'.format(timeStamp),
                                                  whole_nodes=True)

            database = 'db_{}.dat'.format(timeStamp)

#  Collect results of the workflow into the database file.
            if self.services.wait_task(task_wait):
                self.services.error('massive_serial_runner: step failed to run massive serial')

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
            task_wait = self.services.launch_task(1, self.services.get_working_dir(),
                                                  self.TO_JSON_EXE,
                                                  '--input_file={}'.format(database),
                                                  '--output_file={}'.format(self.current_batch),
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
