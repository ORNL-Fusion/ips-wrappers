#! /usr/bin/env python

"""
RAYS_basic_component (Batchelor) 4/30/2024
A simple component script to run the RAYS code from provided input files

"""
import shutil
import get_IPS_config_parameters as config
from ipsframework import Component

class RAYS (Component):
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
# Runs the rays code
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print('RAYS_basic.step() called')
        services = self.services

    # Get global configuration parameters (none for this example)

    # Get component-specific configuration parameters.
        BIN_PATH = config.get_component_param(self, services, 'BIN_PATH')
        NPROC = config.get_component_param(self, services, 'NPROC')
        EXECUTABLE = config.get_component_param(self, services, 'EXECUTABLE')
        RAYS_NML = config.get_component_param(self, services, 'RAYS_NML')

    # Copy state files over to working directory (none for this example)

    # Get input files
        try:
          services.stage_input_files(self.INPUT_FILES)
        except Exception:
          print('Error in call to stage_input_files()')
          self.services.error('Error in call to stage_input_files()')
          raise

    # Copy rays input namelist file to generic name "rays.in" unless the NML name is rays.in
        if RAYS_NML != 'rays.in':
            try:
                shutil.copyfile(RAYS_NML, 'rays.in')
            except IOError as xxx_todo_changeme:
                (errno, strerror) = xxx_todo_changeme.args
                print('Error copying file %s to %s' % (RAYS_NML, 'rays.dat'), strerror)
                services.error('Error copying RAYSNM -> rays.dat')
                raise Exception('Error copying RAYSNM -> rays.dat')


#     Launch RAYS - N.B: Path to executable is in config parameter EXECUTABLE
        print('rf_rays: launching RAYS')
        cwd = services.get_working_dir()
        task_id = services.launch_task(self.NPROC, cwd, self.EXECUTABLE, logfile='log.rays')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            print('Error executing command: ', self.EXECUTABLE)
            services.error('Error executing rays')
            raise Exception('Error executing rays')
        print('rf_rays: finished RAYS')

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
        print('RAYS_basic.checkpoint() called')
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
        print('RAYS_basic finalize() called')
