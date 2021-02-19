#! /usr/bin/env python

# version 0.0 8/7/09 (Batchelor)

# ------------------------------------------------------------------------------
#
# FUS component script to drive model_FUS executable. The executable requires 5
# commandline arguments:
#    1) current state file
#    2) current eqdsk file
#    3) mode = one of "INIT", "STEP", "FINALIZE"
#    4) timeStamp = initial time for "INIT", or = time at end of time stamp for "STEP"
#
# ------------------------------------------------------------------------------

import sys
import os
import subprocess
import getopt
import shutil
import string
from ipsframework import Component

class model_FUS (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp=0):
        print 'py model_FUS.init() called'

        services = self.services


# Copy current state files and input files over to working directory
        self.services.stage_state()
        self.services.stage_input_files(self.INPUT_FILES)
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        prior_state_file = self.services.get_config_param('PRIOR_STATE')
        next_state_file = self.services.get_config_param('NEXT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        FUS_bin = os.path.join(self.BIN_PATH, 'model_FUS')

        print 'Executing ', [FUS_bin, cur_state_file, 'INIT', timeStamp]
        retcode = subprocess.call(['touch','FUS.out'])
        if (retcode != 0):
            print 'Error executing ', FUS_bin
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
        print 'model_FUS.step() called'

        if (self.services == None) :
            print 'Error in model_FUS:;step () : init() function not called before step().'
            return 1
        services = self.services

# Copy current state files and input files over to working directory
        self.services.stage_state()
        self.services.stage_input_files(self.INPUT_FILES)
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        prior_state_file = self.services.get_config_param('PRIOR_STATE')
        next_state_file = self.services.get_config_param('NEXT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

# Call model_FUS
        FUS_bin = os.path.join(self.BIN_PATH, 'model_FUS')

#       print 'Executing ', [FUS_bin, cur_state_file, 'STEP', timeStamp]
#       retcode = subprocess.call([FUS_bin, cur_state_file, cur_eqdsk_file,
#       'STEP', timeStamp])
#       if (retcode != 0):
#           print 'Error executing ', FUS_bin
#           return 1

# Update plasma state
        services.update_state()

# "Archive" output files in history directory
        retcode = subprocess.call(['touch','FUS.out'])
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

# ------------------------------------------------------------------------------
#
# finalize function
#
# Does nothing
# ------------------------------------------------------------------------------



    def finalize(self):
        print 'model_FUS finalize() called'
