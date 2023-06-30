#! /usr/bin/env python

"""
basic_init.py  Batchelor (3-11-2018)

Version 2.0 (Batchelor 10/31/2020)
The purpose is to collect a complete set of initial state files and stage them to the 
/work/state/ directory. It is simplified and adapted from generic_ps_init.py, but 
eliminating reference to many features specific to plasma physics.  The immediate  
application is to simple, example simulations which do not make use of the SWIM Plasma 
State system. The SWIM Plasma State need not be used at all, but can be.

This script first touches all the files listed as STATE_FILES in the config file.  If
config variable INIT_MODE = TOUCH_ONLY that is all that is done, then returns.
If INPUT_FILES are listed in the [init] section of the config file these are then staged
to the working directory thereby overwriting the dummy files generated before.  It may be
convenient to initialize a state file from an input file of a different name, therefore a
service is provided to copy a file to another name before updating state.  This is
controlled by config parameters COPY_FILES and COPIED_FILES_NEW_NAMES

COPY_FILES = list of files to be copied
COPIED_FILES_NEW_NAMES = list of new names for files, must match COPY_FILES

The term CURRENT_STATE refers to a SWIM PLASMA_STATE file if it is being used.  If a 
CURRENT_STATE is specified in the simulation config file, this script adds the Plasma State
variables to it -> run_id and time loop variables tinit, and tfinal.  

If more work needs to be done before the individual components do their own
init, one can specify an INIT_HELPER_CODE (full path) in the config file which
if present will be executed here. If input files are needed for the helper code
they must also be specified in the [init] section of the config file.

"""

import subprocess
import os
# import utils.simple_assignment_file_edit as edit
# import utils.get_IPS_config_parameters as config
import simple_file_editing_functions as edit
import simple_file_editing_functions as config
from ipsframework import Component


class basic_init (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

# ------------------------------------------------------------------------------
#
# init function
#
# Does nothing.
#
# ------------------------------------------------------------------------------

    def init(self, timestamp=0.0):
        print(' ')
        print('basic_init.init() called')
        return

# ------------------------------------------------------------------------------
#
# step function
#
# Calls fortran executable init_empty_state and updates  state
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print(' ')
        print('basic_init.step() called')

        services = self.services

# Check if this is a restart simulation
        simulation_mode = config.get_global_param(
            self, services, 'SIMULATION_MODE')

        if simulation_mode == 'RESTART':
            print('basic_init: RESTART')
        if simulation_mode not in ['RESTART', 'NORMAL']:
            logMsg = 'basic_init: unrecoginzed SIMULATION_MODE: ' + \
                      simulation_mode
            self.services.error(logMsg)
            raise ValueError(logMsg)

# ------------------------------------------------------------------------------
#
# RESTART simulation mode
#
# ------------------------------------------------------------------------------

        if simulation_mode == 'RESTART':
            # Get restart files listed in config file. Here just the  state
            # files.
            restart_root = config.get_global_param(
                self, services, 'RESTART_ROOT')
            restart_time = config.get_global_param(
                self, services, 'RESTART_TIME')
            try:
                services.get_restart_files(
                    restart_root, restart_time, self.RESTART_FILES)
            except BaseException:
                logMsg = 'Error in call to get_restart_files()'
                self.services.exception(logMsg)
                raise

        # Check if there is a config parameter CURRENT_STATE and add data
            cur_state_file = config.get_global_param(
                self, services, 'CURRENT_STATE', optional=True)
            if cur_state_file is not None and len(cur_state_file) > 0:
                timeloop = services.get_time_loop()
                tfinal = timeloop[-1]

                # Put data into current state. For restart only tfinal gets
                # changed
                variable_dict = {'tfinal': tfinal}
                edit.add_variables_to_output_file(
                    variable_dict, cur_state_file)

# ------------------------------------------------------------------------------
#
# NORMAL simulation mode
#
# ------------------------------------------------------------------------------
        else:

            print('basic_init: simulation mode NORMAL')
            state_file_list = config.get_global_param(
                self, services, 'STATE_FILES').split(' ')

        # Generate state files as dummies so framework will have a complete set
            for file in state_file_list:
                print('touching state file = ', file)
                try:
                    subprocess.call(['touch', file])
                except Exception:
                    print('No file ', file)

            init_mode = config.get_component_param(
                self, services, 'INIT_MODE', optional=True)
            if init_mode in ['touch_only', 'TOUCH_ONLY']:
                # Update  state
                try:
                    services.update_state()
                except Exception as e:
                    print('Error in call to updateState()', e)
                    raise
                return

        # Stage input files if any, thereby overwriting the dummy files generated above
            try:
                services.stage_input_files(self.INPUT_FILES)
            except Exception:
                message = 'basic_init: Error in staging input files'
                print(message)
                services.exception(message)
                raise

        # Check if there is a config parameter CURRENT_STATE and add data if
        # so.
            cur_state_file = config.get_global_param(
                self, services, 'CURRENT_STATE', optional=True)
            if cur_state_file is not None and len(cur_state_file) > 0:
                run_id = config.get_global_param(self, services, 'RUN_ID')

                timeloop = services.get_time_loop()
                tinit = timeloop[0]
                tfinal = timeloop[-1]

                # Put data into current state
                variable_dict = {
                    'run_id': run_id,
                    'tinit': tinit,
                    'tfinal': tfinal}
                edit.add_variables_to_output_file(
                    variable_dict, cur_state_file)

       # Execute HELPER_CODES if any
            HELPER_CODES = config.get_component_param(self, services,
                'INIT_HELPER_CODES', optional = True, verbose = True)
            # Check if there are codes to be run
            if HELPER_CODES is not None and len(HELPER_CODES) > 0:
                HELPER_CODES = HELPER_CODES.split(' ')
                for code in HELPER_CODES:
                    cmd = [code]
                    print('Executing ', cmd)
                    services.send_portal_event(event_type='COMPONENT_EVENT',
                                               event_comment=cmd)
                    retcode = subprocess.call(cmd)
                    if (retcode != 0):
                        logMsg = 'Error executing '.join(map(str, cmd))
                        self.services.error(logMsg)
                        raise Exception(logMsg)

        # Copy files to new names if any

            COPY_FILES = config.get_component_param(self, services,
                'COPY_FILES', optional = True)
            COPIED_FILES_NEW_NAMES = config.get_component_param(self, services,
                'COPIED_FILES_NEW_NAMES', optional = True)
            # Check if there are files to be copied
            if COPY_FILES is not None and len(COPY_FILES) > 0:
                COPY_FILES = COPY_FILES.split(' ')
                COPIED_FILES_NEW_NAMES = COPIED_FILES_NEW_NAMES.split(' ')
                # Verify that list lengths are the same
                if len(COPY_FILES) != len(COPIED_FILES_NEW_NAMES):
                    message = ('Error in generic_component init: COPY_FILES and '
                        'COPIED_FILES_NEW_NAMES lists are different lengths')
                    print(message)
                    self.services.error(message)
                    raise Exception(message)

                # Copy the files
                for i in range(len(COPY_FILES)):
                    try:
                        os.system('cp ' + COPY_FILES[i] + ' ' + COPIED_FILES_NEW_NAMES[i])
                    except OSError as xxx_todo_changeme:
                        (errno, strerror) = xxx_todo_changeme.args
                        print('Error copying file %s to %s' % (COPY_FILES[i],
                         COPIED_FILES_NEW_NAMES[i]), strerror)
                        services.error(COPY_FILES[i] + '-> ' + COPIED_FILES_NEW_NAMES[i])
                        raise Exception('Error copying COPY_FILES -> COPIED_FILES_NEW_NAMES')


# Update  state
        try:
            services.update_state()
        except Exception as e:
            print('Error in call to updateState()', e)
            raise

# "Archive" output files in history directory
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        print(' ')

# ------------------------------------------------------------------------------
#
# checkpoint function
#
# Saves  state files to restart directory
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print('basic_init.checkpoint() called')

        services = self.services
        services.stage_state()
        services.save_restart_files(timestamp, self.RESTART_FILES)
        print(' ')

# ------------------------------------------------------------------------------
#
# finalize function
#
# Does nothing
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print('basic_init.finalize() called')
        print(' ')
