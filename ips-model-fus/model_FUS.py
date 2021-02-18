#! /usr/bin/env python

# version 3.0 7/29/2018 (Batchelor)  Merged model_FUS_2_mcmd to model_FUS
# version 2.0 3/7/11 (Batchelor)
# ------------------------------------------------------------------------------
# Added checkpoint/restart. Made exception handling more uniform across all of the
# model compoinents.  Added reference to config parameters (e.g. NPROC) that are used in
# the real components so the config files will be more uniform.  Changed the executable
# launch in STEP to services.launch_task so the model components can be used concurrently
# and will behave more like real components

# version 0.0 3/28/08 (Batchelor)
# ------------------------------------------------------------------------------
#
# FUS component script to drive model_FUS executable. The executable requires 4
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
from  component import Component
from get_IPS_config_parameters import get_global_param, get_component_param

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
        print 'model_FUS.init() called'

        if (self.services == None) :
            services.error('Error in model_FUS init (): No self.services')
            raise Exception('Error in model_FUS init (): No self.services')
        services = self.services

    # Get global configuration parameters
        try:
            cur_state_file = get_global_param(self, services,'CURRENT_STATE')
            cur_eqdsk_file = get_global_param(self, services,'CURRENT_EQDSK')
        except:
            print 'model_FUS: error in getting config parameters'
            services.error('model_FUS: error in getting config parameters')
            raise Exception, 'model_FUS: error in getting config parameters'

    # Get component-specific configuration parameters. Note: Not all of these are
    # used in 'init' but if any are missing we get an exception now instead of
    # later
        try:
            BIN_PATH = self.BIN_PATH
            RESTART_FILES = self.RESTART_FILES
            NPROC = self.NPROC
        except:
            print 'model_FUS init: error getting component-specific config parameters'
            services.error('model_FUS: error getting component-specific\
            config parameters')
            raise Exception, 'model_FUS: error getting model_epa-specific\
            config parameters'

    # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception, e:
          print 'Error in call to stage_state()' , e
          services.error('Error in call to stage_state()')
          raise Exception, 'Error in call to stage_state()'
      
    # Get input files  
        try:
          services.stage_input_files(self.INPUT_FILES)
        except Exception, e:
          print 'Error in call to stage_input_files()' , e
          self.services.error('Error in call to stage_input_files()')
          raise Exception, 'Error in call to stage_input_files()'

        FUS_bin = os.path.join(BIN_PATH, 'model_FUS')

        print 'Executing ', [FUS_bin, cur_state_file, 'INIT', timeStamp]
        
        try:
            retcode = subprocess.call([FUS_bin, cur_state_file, cur_eqdsk_file, 
            'INIT', timeStamp])
        except Exception:
            services.error(' error executing model_FUS_bin')
            raise Exception(' error executing model_FUS_bin')

# Update plasma state files in plasma_state work directory
        try:
          services.update_state()
        except Exception, e:
          print 'Error in call to update_state()', e
          services.error('Error in call to update_state()')
          raise Exception, 'Error in call to update_state()'

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
      print 'model_FUS.restart() called'

      services = self.services
      workdir = services.get_working_dir()

    # Get restart files listed in config file.        
      try:
            restart_root = get_global_param(self, services,'RESTART_ROOT')
            restart_time = get_global_param(self, services,'RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
      except Exception, e:
            print 'Error in call to get_restart_files()' , e
            self.services.error('model_FUS.restart: error in call to get_restart_files()')
            raise Exception, 'model_FUS.restart: error in call to get_restart_files()'
      return 0

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
            services.error('Error in model_FUS init (): No self.services')
            raise Exception('Error in model_FUS init (): No self.services')
        services = self.services

    # Get global configuration parameters
        try:
            cur_state_file = get_global_param(self, services,'CURRENT_STATE')
            cur_eqdsk_file = get_global_param(self, services,'CURRENT_EQDSK')
        except:
            print 'model_FUS: error in getting config parameters'
            services.error('model_FUS: error in getting config parameters')
            raise Exception, 'model_FUS: error in getting config parameters'

    # Get component-specific configuration parameters.
        try:
            BIN_PATH = self.BIN_PATH
            NPROC = self.NPROC
        except:
            print 'model_FUS init: error getting component-specific config parameters'
            services.error('model_FUS: error getting component-specific\
            config parameters')
            raise Exception, 'model_FUS: error getting model_epa-specific\
            config parameters'

    # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception, e:
          print 'Error in call to stage_state()' , e
          services.error('Error in call to stage_state()')
          raise Exception, 'Error in call to stage_state()'
      
    # Get input files  
        try:
          services.stage_input_files(self.INPUT_FILES)
        except Exception, e:
          print 'Error in call to stage_input_files()' , e
          self.services.error('Error in call to stage_input_files()')
          raise Exception, 'Error in call to stage_input_files()'

# Call model_FUS
        FUS_bin = os.path.join(BIN_PATH, 'model_FUS')

        print 'Executing ', [FUS_bin, cur_state_file, 'STEP', timeStamp]
        cwd = os.getcwd()
        task_id = services.launch_task(NPROC, cwd, FUS_bin, cur_state_file,
            cur_eqdsk_file, 'STEP', timeStamp)
        retcode = services.wait_task(task_id)
        partial_file = cwd + '/FUS_' + cur_state_file
        if(retcode != 0):
            services.error(' error executing model_FUS_bin')
            raise Exception(' error executing model_FUS_bin')
            
# Update plasma state
        services.merge_current_state(partial_file)
        print 'merged partial FUS update'
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
            print 'model_FUS.checkpoint() called'
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
        print 'model_FUS finalize() called'

