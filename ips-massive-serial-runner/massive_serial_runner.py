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
import numpy
import netCDF4
from omfit.classes.omfit_eqdsk import OMFITeqdsk

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component init method. This method prepairs the state.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: init')

#  Get config filenames.
        self.current_state = self.services.get_config_param('CURRENT_MSR_STATE')
        self.database_config = self.services.get_config_param('DATABASE_CONFIG')
        self.current_batch = self.services.get_config_param('CURRENT_MSR_BATCH')

#  Stage state.
        self.services.stage_state()

#  Unzip files from the state. Use mode a so files an be read and written to.
        self.zip_ref = ZipState.ZipState(self.current_state, 'a')
        self.zip_ref.extract('inscan')

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
            with open(current_batch, 'r') as json_ref:
                keys = json.load(json_ref).keys()

            task_wait = self.services.launch_task(self.NPROC,
                                                  self.services.get_working_dir(),
                                                  self.MASSIVE_SERIAL_EXE,
                                                  'inscan',
                                                  self.fastran_config,
                                                  MAX(self.NPROC, len(keys)),
                                                  0)

            database = 'db_{}.dat'.format(timeStamp)

#  Collect results of the workflow.
            if (self.services.wait_task(task_wait)):
                self.services.error('massive_serial_runner: step failed to run massive serial')

            task_wait = self.services.launch_task(1, self.services.get_working_dir(),
                                                  self.MAKE_DATABASE_EXE,
                                                  '--rdir=output',
                                                  '--input={}'.format(self.database_config),
                                                  '--output={}'.format(database))
            if (self.services.wait_task(task_wait)):
                self.services.error('massive_serial_runner: step failed to make database')

            task_wait = self.services.launch_task(1, self.services.get_working_dir(),
                                                  self.TO_JSON_EXE,
                                                  '--input_file={}'.format(database),
                                                  '--output_file={}'.format(self.current_batch))

            if (self.services.wait_task(task_wait)):
                self.services.error('massive_serial_runner: step failed to make json')

            self.zip_ref.write(self.current_batch)
            self.zip_ref.set_state(state='updated')

            self.services.stage_output_files(timeStamp, self.current_batch)

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
