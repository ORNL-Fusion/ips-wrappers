#! /usr/bin/env python

# version 0.3 3/9/2011 (Batchelor)
# Eliminated the next_state_file. Now the source power levels as of the end of the time
# step are put into the full and partial plasmas states that get written 
# model_epa_ps_file_init.f90.  This is to be consistent with the new generic_driver.py that 
# no longer copies next_state to current_state at the beginning of each time step.

# version 0.2 1/11/2011 (Batchelor)
# This version also copies INPUT_EQDSK_FILE to the generic file name = 'input_eqdsk_file'

# version 0.1 10/18/2010 (Batchelor)
# This version will either do a merge of a partial plasma state in the 'step' function so that 
# The component can be used concurrently without overwriting other components plasma state
# inputs or it can do an update of the full plasma state in the 'step' function.  This is 
# controlled by the parameter STATE_WRITE_MODE = full/partial in the model_EPA section of the
# IPS config file.  The default (i.e. STATE_WRITE_MODE absent) is merge partial state
#
# In this version the python component copies the input eqdsk file to
# current_eqdsk_file and copies the input state file to the generic name = "input_state_file" 
#  
# The paths to the input plasma state file and input eqdsk files should be specified as input
# files in the [model_EPA] section of the simulation config file
# Also model_EPA config variables "INPUT_PLASMA_STATE_FILE" and "INPUT_EQDSK_FILE"  should be 
# defined  with the names of these input files.

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

class model_EPA(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

# ------------------------------------------------------------------------------
#
# init function
#
# model_EPA init function allocates plasma profiles and initializes rf sources
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp):
        print 'model_epa.init() called'

        if (self.services == None) :
            services.error('Error in model_EPA init (): No self.services')
            raise Exception('Error in model_EPA init (): No self.services')
        services = self.services

    # Get global configuration parameters
        try:
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        except:
            print 'model_epa: error in getting config parameters'
            services.error('model_epa: error in getting config parameters')
            raise Exception, 'model_epa: error in getting config parameters'

    # Get component-specific configuration parameters. Note: Not all of these are
    # used in 'init' but if any are missing we get an exception now instead of
    # later
        try:
            NPROC = self.NPROC
            BIN_PATH = self.BIN_PATH
            INPUT_FILES = self.INPUT_FILES
            OUTPUT_FILES = self.OUTPUT_FILES
            RESTART_FILES = self.RESTART_FILES
            BIN_PATH = self.BIN_PATH
            input_state_file = self.INPUT_STATE_FILE
            input_eqdsk_file = self.INPUT_EQDSK_FILE

        except:
            print 'model_epa init: error getting component-specific config parameters'
            services.error('model_epa: error getting component-specific\
            config parameters')
            raise Exception, 'model_epa: error getting model_epa-specific\
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
          services.stage_input_files(INPUT_FILES)
        except Exception, e:
          print 'Error in call to stageInputFiles()' , e
          services.error('Error in call to stageInputFiles()')
          raise Exception, 'Error in call to stageInputFiles()'
          
     # Copy input_state_file to generic name "input_state_file"
        try:
            subprocess.call(['cp', input_state_file, "input_state_file" ])
        except Exception, e:
          print 'Error in renaming input_state_file' , e
          services.error('Error in renaming input_state_file')
          raise Exception, 'Error in renaming input_state_file'
          raise Exception, 'Error in call to stageInputFiles()'
                    
     # Copy input_eqdsk_file to generic name "input_eqdsk_file"
        try:
            subprocess.call(['cp', input_eqdsk_file, "input_eqdsk_file" ])
        except Exception, e:
          print 'Error in renaming input_eqdsk_file' , e
          services.error('Error in renaming input_eqdsk_file')
          raise Exception, 'Error in renaming input_eqdsk_file'
        
     # Copy input_eqdsk_file to cur_eqdsk_file
        try:
            subprocess.call(['cp', input_eqdsk_file, cur_eqdsk_file ])
        except Exception, e:
          print 'Error in renaming input_eqdsk_file' , e
          services.error('Error in renaming input_eqdsk_file')
          raise Exception, 'Error in renaming input_eqdsk_file'
        
        model_epa_bin = os.path.join(self.BIN_PATH, 'model_epa_ps_file_init')

# Call model_epa
        epa_bin = os.path.join(BIN_PATH, 'model_epa_ps_file_init')
        print 'Executing ', ' '.join([epa_bin, cur_state_file, cur_eqdsk_file, 'INIT', timeStamp])
        cwd = os.getcwd()
        task_id = services.launch_task(NPROC, cwd, epa_bin, cur_state_file,
            cur_eqdsk_file, 'STEP', timeStamp)
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            print 'Error executing command: ', epa_bin
            sys.exit(1)

#         print 'Executing ', [model_epa_bin, cur_state_file, 'INIT', timeStamp]
#         
#         try:
#             retcode = subprocess.call([model_epa_bin, cur_state_file,
#             cur_eqdsk_file, 'INIT', timeStamp])
#         except Exception, e:
#             services.error(' error executing model_epa_bin')
#             raise Exception(' error executing model_epa_bin')

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
      print 'model_epa_ps_init.restart() called'

      services = self.services
      workdir = services.get_working_dir()

    # Get restart files listed in config file.        
      try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
      except Exception, e:
            print 'Error in call to get_restart_files()' , e
            self.services.error('model_epa_ps_init: error in call to get_restart_files()')
            raise Exception, 'model_epa_ps_init: error in call to get_restart_files()'

# ------------------------------------------------------------------------------
#
# STEP function
#
# Does nothing out of the ordinary for a component script 'step' function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'model_epa.step() called'

        if (self.services == None) :
            services.error('Error in model_epa step (): No self.services')
            raise Exception('Error in model_epa step (): No self.services')
        services = self.services

    # Get global configuration parameters
        try:
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        except:
            print 'model_epa: error in getting config parameters'
            services.error('model_epa: error in getting config parameters')
            raise Exception, 'model_epa: error in getting config parameters'

    # Get component-specific configuration parameters.
        try:
            NPROC = self.NPROC
            BIN_PATH = self.BIN_PATH
        except:
            print 'model_epa init: error getting component-specific config parameters'
            services.error('model_epa: error getting component-specific\
            config parameters')
            raise Exception, 'model_epa: error getting model_epa-specific\
            config parameters'

    # Get optional config parameter
        try:
            STATE_WRITE_MODE = self.STATE_WRITE_MODE
        except:
            STATE_WRITE_MODE = 'partial'

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

# Call model_epa
        model_epa_bin = os.path.join(BIN_PATH, 'model_epa_ps_file_init')
        print 'Executing ', [model_epa_bin, cur_state_file, 'STEP', timeStamp]
        cwd = services.get_working_dir()
        task_id  = services.launch_task(NPROC, cwd, model_epa_bin, cur_state_file,
            cur_eqdsk_file, 'STEP', timeStamp)
        retcode = services.wait_task(task_id)
        partial_file = cwd + '/EPA_' + cur_state_file
        if (retcode != 0):
            print 'Error executing ', model_epa_bin
            return 1

# Either update plasma state or merge plasma state depending on STATE_WRITE_MODE

        if STATE_WRITE_MODE == 'full':
        # Update plasma state files in plasma_state work directory
            try:
              services.update_plasma_state()
              print 'Updated model_epa plasma state data ', cur_state_file
            except Exception, e:
              print 'Error in call to update_plasma_state()', e
              services.error('Error in call to update_plasma_state()')
              raise Exception, 'Error in call to update_plasma_state()'
        else:
    # Merge partial plasma state containing updated EPA data
            try:
              cwd = services.get_working_dir()
              partial_file = cwd + '/EPA_' + cur_state_file
              services.merge_current_plasma_state(partial_file, logfile='log.update_state')
              print 'Merged model_epa plasma state data ', partial_file
            except Exception, e:
              print 'Error in call to merge_current_plasma_state(' , partial_file, ')'
              self.services.error('Error in call to merge_current_plasma_state')
              raise Exception, 'Error in call to merge_current_plasma_state'

    # Archive output files
        try:
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
          print 'Error in call to stage_output_files()', e
          self.services.error('Error in call to stage_output_files()')
          raise Exception, 'Error in call to stage_output_files()'
          
        return 0

# ------------------------------------------------------------------------------
#
# checkpoint function
# Saves plasma state files to restart directory
#
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
            print 'model_epa.checkpoint() called'
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
        print 'model_epa finalize() called'
