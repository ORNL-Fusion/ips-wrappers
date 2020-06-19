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
from component import Component
from get_IPS_config_parameters import get_global_param, get_component_param

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
        self.services.stage_plasma_state()

        self.services.stage_input_files(self.INPUT_FILES)
        cur_state_file = get_global_param(self, services,'CURRENT_STATE')
        next_state_file = get_global_param(self, services,'NEXT_STATE')
        cur_eqdsk_file = get_global_param(self, services,'CURRENT_EQDSK')
        
        model_epa_bin = os.path.join(self.BIN_PATH, 'model_epa_ps_file_init')

        print('Executing ', [model_epa_bin, cur_state_file, 'INIT', timeStamp])
        
        try:
            retcode = subprocess.call([model_epa_bin, cur_state_file, next_state_file,
            cur_eqdsk_file, 'INIT', timeStamp])
        except Exception:
            services.error(' error executing model_epa_bin')
            raise Exception(' error executing model_epa_bin')

# Update (original) plasma state
        services.update_plasma_state()

# "Archive" output files in history directory
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

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
            services.error('Error in model_epa step (): No self.services')
            raise Exception('Error in model_epa step (): No self.services')
        services = self.services


# Copy current and prior state over to working directory
        services.stage_plasma_state()

        cur_state_file = services.get_global_param('CURRENT_STATE')
        next_state_file = get_global_param(self, services,'NEXT_STATE')
        cur_eqdsk_file = get_global_param(self, services,'CURRENT_EQDSK')

# Call model_epa
        model_epa_bin = os.path.join(self.BIN_PATH, 'model_epa_ps_file_init')
        
        try:
            retcode = subprocess.call([model_epa_bin, cur_state_file, next_state_file,
            cur_eqdsk_file, 'STEP', timeStamp])
        except Exception:
            services.error(' error executing model_epa_bin')
            raise Exception(' error executing model_epa_bin')

# Update plasma state
        services.update_plasma_state()

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
