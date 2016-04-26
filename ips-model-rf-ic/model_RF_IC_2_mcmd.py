#! /usr/bin/env python

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

        if (self.services == None) :
            services.error('Error in model_NB init (): No self.services')
            raise Exception('Error in model_NB init (): No self.services')
        services = self.services

    # Get global configuration parameters
        try:
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            cur_cql_file = self.services.get_config_param('CURRENT_CQL')
            cur_dql_file = self.services.get_config_param('CURRENT_DQL')
        except:
            print 'model_RF_IC_2: error in getting config parameters'
            services.error('model_RF_IC_2: error in getting config parameters')
            raise Exception, 'model_RF_IC_2: error in getting config parameters'

    # Get component-specific configuration parameters. Note: Not all of these are
    # used in 'init' but if any are missing we get an exception now instead of
    # later
        try:
            BIN_PATH = self.BIN_PATH
            RESTART_FILES = self.RESTART_FILES
            NPROC = self.NPROC
        except:
            print 'model_NB init: error getting component-specific config parameters'
            services.error('model_NB: error getting component-specific\
            config parameters')
            raise Exception, 'model_NB: error getting model_epa-specific\
            config parameters'

    # Copy plasma state files over to working directory
        try:
          services.stage_plasma_state()
        except Exception, e:
          print 'Error in call to stage_plasma_state()' , e
          services.error('Error in call to stage_plasma_state()')
          raise Exception, 'Error in call to stage_plasma_state()'
      
    # Get input files  
        try:
          services.stage_input_files(self.INPUT_FILES)
        except Exception, e:
          print 'Error in call to stageInputFiles()' , e
          self.services.error('Error in call to stageInputFiles()')
          raise Exception, 'Error in call to stageInputFiles()'

        RF_IC_bin = os.path.join(BIN_PATH, 'model_RF_IC_2_mcmd')

        print 'Executing ', [RF_IC_bin, cur_state_file, 'INIT', timeStamp]
        retcode = subprocess.call([RF_IC_bin, cur_state_file, cur_eqdsk_file,
        cur_cql_file, cur_dql_file, 'INIT', timeStamp])
        if (retcode != 0):
            print 'Error executing ', RF_IC_bin
            return 1

# Update plasma state files in plasma_state work directory
        try:
          services.update_plasma_state()
        except Exception, e:
          print 'Error in call to update_plasma_state()', e
          services.error('Error in call to update_plasma_state()')
          raise Exception, 'Error in call to update_plasma_state()'

# "Archive" output files in history directory
        try:
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
          print 'Error in call to stage_output_files()', e
          services.error('Error in call to stage_output_files()')
          raise Exception, 'Error in call to stage_output_files()'

        return 0

# ------------------------------------------------------------------------------
#
# RESTART function
# Gets restart files from restart directory
# Loads the global configuration parameters from the config file
#
# ------------------------------------------------------------------------------
        
    def restart(self, timeStamp):
      print 'model_RF_IC_2.restart() called'

      services = self.services
      workdir = services.get_working_dir()

    # Get restart files listed in config file.        
      try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
      except Exception, e:
            print 'Error in call to get_restart_files()' , e
            self.services.error('model_RF_IC_2.restart: error in call to get_restart_files()')
            raise Exception, 'model_RF_IC_2.restart: error in call to get_restart_files()'
      return 0

# ------------------------------------------------------------------------------
#
# STEP function
#
# Does nothing out of the ordinary for a component script 'step' function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'model_RF_IC_2_mcmd.step() called'

        if (self.services == None) :
            print 'Error in model_RF_IC:;step () : init() function not called before step().'
            return 1
        services = self.services

    # Get global configuration parameters
        try:
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            cur_cql_file = self.services.get_config_param('CURRENT_CQL')
            cur_dql_file = self.services.get_config_param('CURRENT_DQL')
        except:
            print 'model_RF_IC_2_step: error in getting config parameters'
            services.error('model_RF_IC_2: error in getting config parameters')
            raise Exception, 'model_RF_IC_2: error in getting config parameters'

    # Get component-specific configuration parameters.
        try:
            BIN_PATH = self.BIN_PATH
            NPROC = self.NPROC        
        except:
            print 'model_RF_IC_2_step: error getting component-specific config parameters'
            services.error('model_RF_IC_2_step: error getting component-specific\
            config parameters')
            raise Exception, 'model_RF_IC_2_step: error getting model_epa-specific\
            config parameters'

    # Copy plasma state files over to working directory
        try:
          services.stage_plasma_state()
        except Exception, e:
          print 'Error in call to stage_plasma_state()' , e
          services.error('Error in call to stage_plasma_state()')
          raise Exception, 'Error in call to stage_plasma_state()'
      
    # Get input files  
        try:
          services.stage_input_files(self.INPUT_FILES)
        except Exception, e:
          print 'Error in call to stageInputFiles()' , e
          self.services.error('Error in call to stageInputFiles()')
          raise Exception, 'Error in call to stageInputFiles()'

# Call model_RF_IC
        RF_IC_bin = os.path.join(BIN_PATH, 'model_RF_IC_2_mcmd')

        print 'Executing ', [RF_IC_bin, cur_state_file, 'STEP', timeStamp]
        cwd = services.get_working_dir()
        task_id  = services.launch_task(NPROC, cwd, RF_IC_bin, cur_state_file, cur_eqdsk_file,
        cur_cql_file, cur_dql_file, 'STEP', timeStamp)
        retcode = services.wait_task(task_id)
        partial_file = cwd + '/RF_IC_' + cur_state_file
        if (retcode != 0):
            print 'Error executing ', RF_IC_bin
            return 1

# Update plasma state
        services.merge_current_plasma_state(partial_file)
        print 'merged partial RF_IC update'
# "Archive" output files in history directory
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        return 0

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
