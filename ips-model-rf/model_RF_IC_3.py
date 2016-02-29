#! /usr/bin/env python

# version 3.0 2/28/2015 (Batchelor)
# Worked on the exception handling

# version 2.0 3/7/11 (Batchelor)
# ------------------------------------------------------------------------------
# Added checkpoint/restart. Made exception handling more uniform across all of the
# model compoinents.  Added reference to config parameters (e.g. NPROC) that are used in
# the real components so the config files will be more uniform.  Changed the executable
# launch in STEP to services.launch_task so the model components can be used concurrently
# and will behave more like real components

# version 0.0 3/28/08 (Batchelor)
# ------------------------------------------------------------------------------
# RF_IC component script to drive change_ICRF_profiles executable. The executable requires 3
# commandline arguments:
#       1) path to the current plasma state file
#       2) path to the current plasma eqdsk file
#       3) path to current cql distribution function file - cur_cql_file
#       4) path to current quasilinear operator file - cur_dql_file
#       5) action mode, i.e. one of: "INIT", "STEP", or "FINALIZE"
#       6) time stamp the time set by the driver component to which the simulation is 
#          supposed to advance.
#
# ------------------------------------------------------------------------------

import sys
import os
import subprocess
import getopt
import shutil
import string
from component import Component

class model_RF_IC_3 (Component):
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

        if (self.services == None) :
            services.error('Error in model_NB init (): No self.services')
            raise Exception('Error in model_NB init (): No self.services')
        services = self.services

    # Get global configuration parameters
        cur_state_file = self.try_get_config_param(services,'CURRENT_STATE')
        cur_eqdsk_file = self.try_get_config_param(services,'CURRENT_EQDSK')
        cur_cql_file = self.try_get_config_param(services,'CURRENT_CQL')
        cur_dql_file = self.try_get_config_param(services,'CURRENT_DQL')

    # Get component-specific configuration parameters. Note: Not all of these are
    # used in 'init' but if any are missing we get an exception now instead of
    # later
        BIN_PATH = self.try_get_config_param(services,'BIN_PATH')
        RESTART_FILES = self.try_get_config_param(services,'RESTART_FILES')
        NPROC = self.try_get_config_param(services,'NPROC')

    # Copy plasma state files over to working directory
        try:
          services.stage_plasma_state()
        except Exception, e:
          print 'Error in call to stage_plasma_state()' , e
          services.error('Error in call to stage_plasma_state()')
          raise
      
    # Get input files  
        try:
          services.stage_input_files(self.INPUT_FILES)
        except Exception, e:
          print 'Error in call to stageInputFiles()' , e
          self.services.error('Error in call to stageInputFiles()')
          raise

        RF_IC_bin = os.path.join(BIN_PATH, 'model_RF_IC_3')

        print 'Executing ', [RF_IC_bin, cur_state_file, 'INIT', timeStamp]
        try:
            retcode = subprocess.call([RF_IC_bin, cur_state_file, cur_eqdsk_file,
                cur_cql_file, cur_dql_file, 'INIT', timeStamp])     
        except Exception: 
            message = "Error executing " +  RF_IC_bin
            self.services.error(message)  
            raise
        else: 
            if (retcode != 0):
                message = "Abnormal termination of " + RF_IC_bin
                print message
                self.services.error(message)
                raise Exception(message)

# Update plasma state files in plasma_state work directory
        try:
          services.update_plasma_state()
        except Exception:
          message = 'Error in call to update_plasma_state()'
          print message
          services.error(message)
          raise

# "Archive" output files in history directory

        try:
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception:
          message = 'Error in call to stage_output_files()'
          print message
          services.error(message)
          raise

        return

# ------------------------------------------------------------------------------
#
# RESTART function
# Gets restart files from restart directory
# Loads the global configuration parameters from the config file
#
# ------------------------------------------------------------------------------
        
    def restart(self, timeStamp):
        print 'model_RF_IC_3.restart() called'

        services = self.services
        workdir = services.get_working_dir()

        # Get restart files listed in config file.        
        restart_root = self.try_get_config_param(services,'RESTART_ROOT')
        restart_time = self.try_get_config_param(services,'RESTART_TIME')

        try:
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
        except Exception:
            message = 'Error in call to get_restart_files()'
            print message
            self.services.error(message)
            raise
        return 0

# ------------------------------------------------------------------------------
#
# STEP function
#
# Does nothing out of the ordinary for a component script 'step' function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'model_RF_IC_3.step() called'

        if (self.services == None) :
            services.error('Error in model_NB init (): No self.services')
            raise Exception('Error in model_NB init (): No self.services')
        services = self.services

    # Get global configuration parameters
        INPUT_STATE_FILE = self.try_get_config_param(services,'INPUT_STATE_FILE')
        cur_state_file = self.try_get_config_param(services,'CURRENT_STATE')
        cur_eqdsk_file = self.try_get_config_param(services,'CURRENT_EQDSK')
        cur_cql_file = self.try_get_config_param(services,'CURRENT_CQL')
        cur_dql_file = self.try_get_config_param(services,'CURRENT_DQL')

    # Get component-specific configuration parameters.
        BIN_PATH = self.try_get_config_param(services,'BIN_PATH')
        NPROC = self.try_get_config_param(services,'NPROC ')      

    # Copy plasma state files over to working directory
        try:
          services.stage_plasma_state()
        except Exception, e:
          print 'Error in call to stage_plasma_state()' , e
          services.error('Error in call to stage_plasma_state()')
          raise
      
    # Get input files  
        try:
          services.stage_input_files(self.INPUT_FILES)
        except Exception:
          print 'Error in call to stageInputFiles()'
          self.services.error('Error in call to stageInputFiles()')
          raise


# Call model_RF_IC
        RF_IC_bin = os.path.join(BIN_PATH, 'model_RF_IC_3')

        print 'Executing ', [RF_IC_bin, cur_state_file, 'STEP', timeStamp]
        cwd = services.get_working_dir()
        task_id  = services.launch_task(NPROC, cwd, RF_IC_bin, cur_state_file, cur_eqdsk_file,
        cur_cql_file, cur_dql_file, 'STEP', timeStamp)
        retcode = services.wait_task(task_id)
        partial_file = cwd + '/RF_IC_' + cur_state_file
        if (retcode != 0):
            message = 'Error executing ', RF_IC_bin
            self.services.error(message)
            raise Exception(message)
            return 1

# Update plasma state files in plasma_state work directory
        try:
            services.merge_current_plasma_state(partial_file)
        except Exception:
            message = 'Error in call merge_current_plasma_state()'
            print message
            services.error(message)
            raise

# "Archive" output files in history directory

        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception:
            message = 'Error in call to stage_output_files()'
            print message
            services.error(message)
            raise

        return

# ------------------------------------------------------------------------------
#
# checkpoint function
# Saves plasma state files to restart directory
#
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
            print 'model_RF_IC_2.checkpoint() called'
            services = self.services
            services.save_restart_files(timestamp, self.RESTART_FILES)
            return 0

# ------------------------------------------------------------------------------
#
# finalize function
#
# Does nothing
# ------------------------------------------------------------------------------



    def finalize(self, timestamp=0.0):
        print 'model_RF_IC_2 finalize() called'

# "Private"  methods

    # Try to get config parameter - wraps the exception handling for get_config_parameter()
    def try_get_config_param(self, services, param_name, optional=False):

        try:
            value = services.get_config_param(param_name)
            print param_name, ' = ', value
        except Exception :
            if optional: 
                print 'config parameter ', param_name, ' not found'
                value = None
            else:
                message = 'required config parameter ', param_name, ' not found'
                print message
                services.exception(message)
                raise
        
        return value
