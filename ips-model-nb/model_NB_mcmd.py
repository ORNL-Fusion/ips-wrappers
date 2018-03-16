#! /usr/bin/env python

# version 0.0 3/28/08 (Batchelor)

# ------------------------------------------------------------------------------
#
# NB component script to drive model_NB executable. The executable requires 5
# commandline arguments:
#    1) current state file
#    2) current eqdsk file
#    3) current jsdsk file
#    4) mode = one of "INIT", "STEP", "FINALIZE"
#    5) timeStamp = initial time for "INIT", or = time at end of time stamp for "STEP"
#
# ------------------------------------------------------------------------------

import sys
import os
import subprocess
import getopt
import shutil
import string
from component import Component

class model_NB (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp=0):
        print 'model_NB.init() called'

        services = self.services


# Copy current state files and input files over to working directory
        self.services.stage_plasma_state()
        self.services.stage_input_files(self.INPUT_FILES)
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        prior_state_file = self.services.get_config_param('PRIOR_STATE')
        next_state_file = self.services.get_config_param('NEXT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_jsdsk_file = self.services.get_config_param('CURRENT_JSDSK')

        NB_bin = os.path.join(self.BIN_PATH, 'model_NB')

        print 'Executing ', [NB_bin, cur_state_file, 'INIT', timeStamp]
        retcode = subprocess.call([NB_bin, cur_state_file, cur_eqdsk_file,
        cur_jsdsk_file, 'INIT', timeStamp])
        pwd = os.getcwd()
#     task_id_NB = services.launch_task(1, pwd, NB_bin, cur_state_file, cur_eqdsk_file,
#        cur_jsdsk_file, 'INIT', timeStamp)
#     retcode = wait_task(task_id_NB)
        if (retcode != 0):
            print 'Error executing ', NB_bin, 'for nb init'
            return 1

# Update (original) plasma state
        services.update_plasma_state()

        return

# ------------------------------------------------------------------------------
#
# STEP function
#
# Does nothing out of the ordinary for a component script 'step' function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'model_NB.step() called'

        if (self.services == None) :
            print 'Error in model_NB:;step () : init() function not called before step().'
            sys.exit(1)
        services = self.services

# Copy current state files and input files over to working directory
        self.services.stage_plasma_state()
        self.services.stage_input_files(self.INPUT_FILES)
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        prior_state_file = self.services.get_config_param('PRIOR_STATE')
        next_state_file = self.services.get_config_param('NEXT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_jsdsk_file = self.services.get_config_param('CURRENT_JSDSK')


# Call model_NB
        NB_bin = os.path.join(self.BIN_PATH, 'model_NB')

        print 'Executing ', [NB_bin, cur_state_file, 'STEP', timeStamp]
        pwd = os.getcwd()
#     retcode = subprocess.call([NB_bin, cur_state_file, cur_eqdsk_file,
#        cur_jsdsk_file, 'STEP', timeStamp])
        task_id_NB = services.launch_task(1, pwd, NB_bin, cur_state_file, cur_eqdsk_file,
           cur_jsdsk_file, 'STEP', timeStamp)
        retcode = services.wait_task(task_id_NB)
        print 'Launched model NB task'
        if (retcode != 0):
            print 'Error executing ', NB_bin
            return 1

# Update plasma state
        services.update_plasma_state()

# "Archive" output files in history directory
        print 'ERROR_3'
        print type(timeStamp), type(self.OUTPUT_FILES)
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

# ------------------------------------------------------------------------------
#
# finalize function
#
# Does nothing
# ------------------------------------------------------------------------------



    def finalize(self, timestamp=0.0):
        print 'model_NB finalize() called'
