#! /usr/bin/env python

"""
CQL3D_basic_component (Batchelor) 5/17/2018
A simple component script to run the CQL3D code from provided input files

"""
import shutil
import get_IPS_config_parameters as config
from ipsframework import Component

class CQL3D_basic (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

# ------------------------------------------------------------------------------
#
# init function
#
# Since there are no state files for this example and the IPS regards "init" as
# optional, there is no init function.
#
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#
# STEP function
#
# Runs the CQL3D code  
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print('CQL3D_basic.step() called')
        services = self.services

    # Get global configuration parameters (none for this example)
 
    # Get component-specific configuration parameters.
        BIN_PATH = config.get_component_param(self, services, 'BIN_PATH')
        NPROC = config.get_component_param(self, services, 'NPROC')
        EXECUTABLE = config.get_component_param(self, services, 'EXECUTABLE')

    # Copy state files over to working directory (none for this example)
      
    # Get input files  
        try:
          services.stage_input_files(self.INPUT_FILES)
        except Exception:
          print('Error in call to stage_input_files()')
          self.services.error('Error in call to stage_input_files()')
          raise
      
#     Launch CQL3D - N.B: Path to executable is in config parameter EXECUTABLE
        print('rf_CQL3D: launching CQL3D')
        cwd = services.get_working_dir()
        task_id = services.launch_task(self.NPROC, cwd, self.EXECUTABLE, logfile='log.CQL3D')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            print('Error executing command: ', CQL3D_bin)
            services.error('Error executing CQL3D')
            raise Exception('Error executing CQL3D')
        print('rf_CQL3D: finished CQL3D')

# "Archive" output files in history directory
        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception:
            message = 'Error in call to stage_output_files()'
            print(message)
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
        print('CQL3D_basic.checkpoint() called')
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
        print('CQL3D_basic finalize() called')
        