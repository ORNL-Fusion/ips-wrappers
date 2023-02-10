#! /usr/bin/env python

# version 0.4 4/27/08 (Batchelor)

# ------------------------------------------------------------------------------
#
# EPA component script to drive model executable.
# ! This version generates a next plasma state (ps_next) containing the source parameters
# ! for the next time step (presently this only consists of ps%power_IC(1) ).  The stored
# ! current plasma state contains the source parameters for the current time step.
#
# !      The executable requires 5 commandline arguments:
# !      1) current state file!
# !      2) Next state file
# !      3) current eqdsk file (for completeness, this version doesn't use eqdsk file)
# !      4) mode = one of "INIT", "STEP", "FINALIZE"
# !      5) timeStamp = initial time for "INIT", or = time at end of time stamp for "STEP"
#
# ------------------------------------------------------------------------------

import sys
import os
import subprocess
import getopt
import shutil
import string
from ipsframework import Component

class model_EPA(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

# ------------------------------------------------------------------------------
#
# init function
#
# model_EPA init function allocates plasma profiles and initializes rf sources
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp):
        print('model_epa.init() called')

        services = self.services

# Copy current and prior state over to working directory
        self.services.stage_state()

        self.services.stage_input_files(self.INPUT_FILES)
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        next_state_file = self.services.get_config_param('NEXT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        change_bin = os.path.join(self.BIN_PATH, 'model_epa')

        print('Executing ', [change_bin, cur_state_file, 'INIT', timeStamp])
        retcode = subprocess.call([change_bin, cur_state_file, next_state_file,
        cur_eqdsk_file, 'INIT', timeStamp])
        if (retcode != 0):
            print('Error executing ', change_bin)
            return 1

# Update (original) plasma state
        services.update_state()

        return

# ------------------------------------------------------------------------------
#
# STEP function
#
# Does nothing out of the ordinary for a component script 'step' function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print('model_epa.step() called')

        if (self.services == None) :
            print('Error in model_epa: step () : init() function not called before step().')
            sys.exit(1)
        services = self.services


# Copy current and prior state over to working directory
        services.stage_state()

        cur_state_file = services.get_config_param('CURRENT_STATE')
        next_state_file = self.services.get_config_param('NEXT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

# Call model_epa
        change_bin = os.path.join(self.BIN_PATH, 'model_epa')
        retcode = subprocess.call([change_bin, cur_state_file, next_state_file,
        cur_eqdsk_file, 'STEP', timeStamp])
        if (retcode != 0):
            print('Error executing command: ', change_bin)
            sys.exit(1)

# Update plasma state
        services.update_state()

# "Archive" output files in history directory
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

# ------------------------------------------------------------------------------
#
# finalize function
#
# Does nothing
# ------------------------------------------------------------------------------



    def finalize(self, timestamp=0.0):
        print('model_epa finalize() called')
