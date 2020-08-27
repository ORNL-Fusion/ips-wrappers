#! /usr/bin/env python

"""
GENRAY_basic_component (Batchelor) 5/17/2018
A simple component script to run the GENRAY code from provided input files

"""
import shutil
import get_IPS_config_parameters as config
from component import Component

class GENRAY_basic (Component):
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
# Runs the genray code  
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print('GENRAY_basic.step() called')
        services = self.services

    # Get global configuration parameters (none for this example)
 
    # Get component-specific configuration parameters.
        BIN_PATH = config.get_component_param(self, services, 'BIN_PATH')
        NPROC = config.get_component_param(self, services, 'NPROC')
        EXECUTABLE = config.get_component_param(self, services, 'EXECUTABLE')
        GENRAYNML = config.get_component_param(self, services, 'GENRAYNML')

    # Copy state files over to working directory (none for this example)
      
    # Get input files  
        try:
          services.stage_input_files(self.INPUT_FILES)
        except Exception:
          print('Error in call to stageInputFiles()')
          self.services.error('Error in call to stageInputFiles()')
          raise

    # Copy genray input namelist file to generic name "genray.dat"
        try:
            shutil.copyfile(GENRAYNML, 'genray.dat')
        except IOError as xxx_todo_changeme:
            (errno, strerror) = xxx_todo_changeme.args
            print('Error copying file %s to %s' % (GENRAYNML, 'genray.dat'), strerror)
            services.error('Error copying GENRAYNM -> genray.dat')
            raise Exception('Error copying GENRAYNM -> genray.dat')

      
#     Launch GENRAY - N.B: Path to executable is in config parameter EXECUTABLE
        print('rf_genray: launching GENRAY')
        cwd = services.get_working_dir()
        task_id = services.launch_task(self.NPROC, cwd, self.EXECUTABLE, logfile='log.genray')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            print('Error executing command: ', self.EXECUTABLE)
            services.error('Error executing genray')
            raise Exception('Error executing genray')
        print('rf_genray: finished GENRAY')

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
        print('GENRAY_basic.checkpoint() called')
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
        print('GENRAY_basic finalize() called')
        