#! /usr/bin/env python

# Version 1.0 of new model_RF_IC 3/6/2018 (Batchelor)
# Eliminated a large number of previous versions of this component and reverted name from
# model_RF_IC_3 to model_RF_IC. 

# Version 4.0 3/6/2018 (Batchelor)
# Eliminated reference to cur_cql_file and cur_dql_file since they are not used here.
# Also generally cleaned up.  Now the fortran executable model_RF_IC requires 4
# commandline arguments:
#       1) path to the current plasma state file
#       2) path to the current plasma eqdsk file
#       3) action mode, i.e. one of: "INIT", "STEP", or "FINALIZE"
#       4) time stamp the time set by the driver component to which the simulation is 
#          supposed to advance.

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
# RF_IC component script to drive change_ICRF_profiles executable. The fortran executable
# model_RF_IC_3 requires 6
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
from ipsframework import Component

class model_RF_IC (Component):
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
            message = 'Error in model_RF_IC init (): No self.services'
            print message
            services.error(message)
            raise
        services = self.services

    # Get global configuration parameters
        cur_state_file = self.get_config_param(services,'CURRENT_STATE')
        cur_eqdsk_file = self.get_config_param(services,'CURRENT_EQDSK')

    # Get component-specific configuration parameters. Note: Not all of these are
    # used in 'init' but if any are missing we get an exception now instead of
    # later
        BIN_PATH = self.get_component_param(services, 'BIN_PATH')
        RESTART_FILES = self.get_component_param(services, 'RESTART_FILES')
        NPROC = self.get_component_param(services, 'NPROC')

    # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception, e:
          print 'Error in call to stage_state()' , e
          services.error('Error in call to stage_state()')
          raise
      
    # Get input files  
        try:
          services.stage_input_files(self.INPUT_FILES)
        except Exception, e:
          print 'Error in call to stage_input_files()' , e
          self.services.error('Error in call to stage_input_files()')
          raise

        RF_IC_bin = os.path.join(BIN_PATH, 'model_RF_IC')

    # Run model_RF_IC fortran
        cmd = [RF_IC_bin, cur_state_file, cur_eqdsk_file, 'INIT', timeStamp]
        print 'Executing = ', cmd
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
          event_comment =  cmd)
        retcode = subprocess.call(cmd)
        if (retcode != 0):
            logMsg = 'Error executing '.join(map(str, cmd))
            self.services.error(logMsg)
            raise Exception(logMsg)

# Update plasma state files in plasma_state work directory
        try:
          services.update_state()
        except Exception:
          message = 'Error in call to update_state()'
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
        print 'model_RF_IC.restart() called'

        if (self.services == None) :
            message = 'Error in model_RF_IC init (): No self.services'
            print message
            services.error(message)
            raise
        services = self.services
        workdir = services.get_working_dir()

        # Get restart files listed in config file.        
        restart_root = self.get_config_param(services,'RESTART_ROOT')
        restart_time = self.get_config_param(services,'RESTART_TIME')

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
        print 'model_RF_IC.step() called'

        if (self.services == None) :
            message = 'Error in model_RF_IC init (): No self.services'
            print message
            services.error(message)
            raise
        services = self.services

    # Get global configuration parameters
        cur_state_file = self.get_config_param(services,'CURRENT_STATE')
        cur_eqdsk_file = self.get_config_param(services,'CURRENT_EQDSK')

    # Get component-specific configuration parameters.
        BIN_PATH = self.get_component_param(services, 'BIN_PATH')
        NPROC = self.get_component_param(services, 'NPROC', optional = True)      

    # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception, e:
          print 'Error in call to stage_state()' , e
          services.error('Error in call to stage_state()')
          raise
      
    # Get input files  
        try:
          services.stage_input_files(self.INPUT_FILES)
        except Exception:
          print 'Error in call to stage_input_files()'
          self.services.error('Error in call to stage_input_files()')
          raise


# Call model_RF_IC
        RF_IC_bin = os.path.join(BIN_PATH, 'model_RF_IC')

        print 'Executing ', [RF_IC_bin, cur_state_file, 'STEP', timeStamp]
        cwd = services.get_working_dir()
        task_id  = services.launch_task(NPROC, cwd, RF_IC_bin, cur_state_file, cur_eqdsk_file,
        'STEP', timeStamp)
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            message = 'Error executing ', RF_IC_bin
            print message
            self.services.error(message)
            raise Exception(message)
            return 1
        partial_file = cwd + '/RF_IC_' + cur_state_file

# Update plasma state files in plasma_state work directory
        try:
            services.merge_current_state(partial_file)
        except Exception:
            message = 'Error in call merge_current_state()'
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
        print 'model_RF_IC.checkpoint() called'
        if (self.services == None) :
            message = 'Error in model_RF_IC init (): No self.services'
            print message
            services.error(message)
            raise
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
        print 'model_RF_IC finalize() called'

# ------------------------------------------------------------------------------
#
# "Private"  methods
#
# ------------------------------------------------------------------------------


    # Try to get config parameter - wraps the exception handling for get_config_parameter()
    def get_config_param(self, services, param_name, optional=False):

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

    # Try to get component specific config parameter - wraps the exception handling
    def get_component_param(self, services, param_name, optional=False):

        if hasattr(self, param_name):
            value = getattr(self, param_name)
            print param_name, ' = ', value
        elif optional:
            print 'optional config parameter ', param_name, ' not found'
            value = None
        else:
            message = 'required component config parameter ', param_name, ' not found'
            print message
            services.exception(message)
            raise
        
        return value

