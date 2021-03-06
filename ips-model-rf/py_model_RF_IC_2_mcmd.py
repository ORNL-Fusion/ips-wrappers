#! /usr/bin/env python

# version 0.0 3/28/08 (Batchelor)

# ------------------------------------------------------------------------------
#
# RF_IC component script to drive change_ICRF_profiles executable. The executable requires 3
# commandline arguments:
# 1) current state file
# 2) mode = one of "INIT", "STEP", "FINALIZE"
# 3) timeStamp = initial time for "INIT", or = time at end of time stamp for "STEP"
#
# ------------------------------------------------------------------------------

import sys
import os
import subprocess
import getopt
import shutil
import string
from component import Component

class model_RF_IC_2 (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp=0):
        print 'model_RF_IC.init() called'

        services = self.services


# Copy current state files and input files over to working directory
        self.services.stage_plasma_state()
        self.services.stage_input_files(self.INPUT_FILES)
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_cql_file = self.services.get_config_param('CURRENT_CQL')
        cur_dql_file = self.services.get_config_param('CURRENT_DQL')

        RF_IC_bin = os.path.join(self.BIN_PATH, 'py_model_RF_IC_2_mcmd')

        print 'Executing ', [RF_IC_bin, cur_state_file, 'INIT', timeStamp]
#       retcode = subprocess.call([RF_IC_bin, cur_state_file, cur_eqdsk_file,
#       cur_cql_file, cur_dql_file, 'INIT', timeStamp])
        retcode = subprocess.Popen(['touch','test.file'])
        print 'ret code for first touch ', retcode
        retcode = subprocess.Popen(['touch',cur_state_file])
        print 'ret code for second touch ', retcode
        if (retcode != 0):
            print 'Error executing ', RF_IC_bin
            return 1
# Update (original) plasma state
        services.update_plasma_state()
# "Archive" output files in history directory
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        return 0

# ------------------------------------------------------------------------------
#
# STEP function
#
# Does nothing out of the ordinary for a component script 'step' function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'py model_RF_IC_2_mcmd.step() called'

        if (self.services == None) :
            print 'Error in model_RF_IC:;step () : init() function not called before step().'
            return 1
        services = self.services

# Copy current state files and input files over to working directory
        self.services.stage_plasma_state()
        self.services.stage_input_files(self.INPUT_FILES)
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_cql_file = self.services.get_config_param('CURRENT_CQL')
        cur_dql_file = self.services.get_config_param('CURRENT_DQL')
        cwd = os.getcwd()

# Call model_RF_IC
        RF_IC_bin = os.path.join(self.BIN_PATH, 'py_model_RF_IC_2_mcmd')

        print 'Executing ', [RF_IC_bin, cur_state_file, 'STEP', timeStamp]
        task_id  = services.launch_task(1, cwd, RF_IC_bin)
        retcode = services.wait_task(task_id)
        partial_file = cwd + '/RF_IC_' + cur_state_file
        if (retcode != 0):
            print 'Error executing ', RF_IC_bin
            return 1

# Update plasma state
#       services.merge_current_plasma_state(partial_file)

# "Archive" output files in history directory
        subprocess.Popen(['touch',self.OUTPUT_FILES])
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        return partial_file
# ------------------------------------------------------------------------------
#
# finalize function
#
# Does nothing
# ------------------------------------------------------------------------------



    def finalize(self, timestamp=0.0):
        print 'model_RF_IC finalize() called'
