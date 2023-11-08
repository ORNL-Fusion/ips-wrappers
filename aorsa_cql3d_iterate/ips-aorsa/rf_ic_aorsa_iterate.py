"""
Rough AORSA component. SF 10/2023
"""

from builtins import range
import sys
import os
import subprocess
import getopt
import shutil
import string
import h5py
import numpy as np
from netCDF4 import *
import scipy.io as spio
from  ipsframework import Component
from simple_file_editing_functions import put_lines

class aorsa (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp=0):
        print('\ntoric.init() called')

        services = self.services
        workdir = services.get_working_dir()

        # Get global configuration parameters
        cur_state_file = self.try_get_config_param(services,'CURRENT_STATE')
        cur_eqdsk_file = self.try_get_config_param(services,'CURRENT_EQDSK')
        cur_cql_file = self.try_get_config_param(services,'CURRENT_CQL')

        # Get component-specific configuration parameters. Note: Not all of these are
        # used in 'init' but if any are missing we get an exception now instead of
        # later
        try:
            NPROC = self.NPROC
            BIN_PATH = self.BIN_PATH
            AORSA_BIN=self.AORSA_BIN
            INPUT_DIR = self.INPUT_DIR
            INPUT_FILES = self.INPUT_FILES
            OUTPUT_FILES = self.OUTPUT_FILES
            RESTART_FILES = self.RESTART_FILES
        except:
            print('rf_ic_aorsa_iterate init: error getting aorsa-specific config parameters')
            services.error('rf_ic_aorsa_iterate: error getting aorsa-specific\
            config parameters')
            raise Exception('rf_ic_aorsa_iterate: error getting aorsa-specific\
            config parameters')
        
        # Copy plasma state files over to working directory
        try:
            services.stage_state()
        except Exception:
            logMsg = 'Error in call to stage_plasma_state()'
            self.services.exception(logMsg)
            raise 

        # Get input files
        try:
            services.stage_input_files(self.INPUT_FILES)
        except:
            logMsg = 'Error in call to stageInputFiles()'
            self.services.exception(logMsg)
            raise
        
        return 0
    
# ------------------------------------------------------------------------------
#
# RESTART function
# Gets restart files from restart directory
# Loads the global configuration parameters from the config file
#
# ------------------------------------------------------------------------------

    def restart(self, timeStamp):
        print('\naorsa.restart() called')

        services = self.services
        workdir = services.get_working_dir()

        # Get restart files listed in config file.
        try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
        except Exception:
            logMsg = 'Error in call to get_restart_files()'
            self.services.exception(logMsg)
            raise

        # Get global configuration parameters
        try:
            self.plasma_state_file = services.get_config_param('CURRENT_STATE')
            self.eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            self.toric_log = os.path.join(workdir, 'log.aorsa')
        except:
            logMsg = 'toric restart: error in getting config parameters'
            self.services.exception(logMsg)
            raise 
            
        return 0
    
# ------------------------------------------------------------------------------
#
# Prepare an input file and run AORSA
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp, **kwargs):
        """Take a step for the aorsa component.  Really a complete run."""
        print('\naorsa.step() called')

        # Ensure that everything has actually been initialized properly
        if (self.services == None):
            logMsg = 'Error in aorsa: step (): No self.services'
            self.services.error(logMsg)
            raise Exception(logMsg)
        services = self.services
        workdir = services.get_working_dir()

        # Copy plasma state files over to working directory
        try:
            services.stage_state()
        except:
            logMsg = 'Error in call to stage_plasma_state()'
            self.services.exception(logMsg)
            raise

        # Get input files
        try:
            services.stage_input_files(self.INPUT_FILES)
        except:
            logMsg = 'Error in call to stageInputFiles()'
            self.services.exception(logMsg)
            raise

        # Set preparation and process i/o script paths
        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_aorsa_input')
        #process_output_bin  = os.path.join(self.BIN_PATH, 'process_aorsa_output') SF not ready yet

        # Set physics exec bin locations
        aorsa_bin = self.AORSA_BIN
        # Get global configuration parameters
        cur_state_file = self.try_get_config_param(services,'CURRENT_STATE')
        cur_eqdsk_file = self.try_get_config_param(services,'CURRENT_EQDSK')
        cur_cql_file = self.try_get_config_param(services,'CURRENT_CQL')
        specs = self.try_get_config_param(services,'SPECS')
        arg_icrfpower = self.try_get_config_param(services,'RFPWR_IC')

        # Set working directory and output path
        cwd = os.getcwd()
        aorsa_log = os.path.join(workdir, 'log.aorsa')

        # Generate namelist for input preparation script
        nml_lines = ['&aorsa_prepare_nml \n']
        nml_lines.append(' cur_state_file = \"' + cur_state_file + '\",\n')
        nml_lines.append(' cur_geq_file = \"' + cur_eqdsk_file + '\",\n')
        nml_lines.append(' cur_cql_file = \"' + cur_cql_file + '\",\n')
        nml_lines.append(' arg_specs = \"' + specs + '\",\n')
        nml_lines.append(' arg_src_indx = 1,\n')
        nml_lines.append(' arg_nphi_indx = 0,\n')
        nml_lines.append(' /\n')
        put_lines('aorsa_prepare.nml',nml_lines)

        # Call aorsa prepare_input to generate aorsa2d.in
        if not os.path.isfile(prepare_input_bin):
            logMsg = 'Cannot find aorsa prepare_input binary: ' + prepare_input_bin
            self.services.error(logMsg)
            raise Exception(logMsg)

        aorsa_prep_log = open('log_prepare_aorsa_input','w')
            
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                                   event_comment =  prepare_input_bin)
        retcode = subprocess.call(prepare_input_bin, stdout = aorsa_prep_log,\
                                  stderr = subprocess.STDOUT)
        if (retcode != 0):
            logMsg = 'Error executing ' + prepare_input_bin
            self.services.error(logMsg)
            raise Exception(logMsg)
        
        cwd = services.get_working_dir()
        
        #run aorsa
        run_nproc = self.NPROC
        run_nppn = self.NPPN
        coresPerProc = int(2 * int(128/int(run_nppn)))
        task_id = services.launch_task(run_nproc, cwd, aorsa_bin, logfile=aorsa_log,\
                                           whole_nodes=True, task_ppn = run_nppn, task_cpp = coresPerProc)
        retcode = services.wait_task(task_id, timeout = 1600 ,delay=20.0)
        if (retcode != 0):
            services.error("AORSA run failed")
            raise Exception("AORSA run failed")

        #output process script would go here

        #copy diffusion coefficients to the correct names
        if specs.strip() == 'MIN':
            shutil.copyfile('out_cql3d.coef2', 'du0u0_input')
        elif specs.strip() == 'MIN+':
            shutil.copyfile('out_cql3d.coef2', 'du0u0_input_1')
            shutil.copyfile('out_cql3d.coef1', 'du0u0_input_2')
        
        #update the plasma state
        try:
            services.update_state()
        except Exception:
            logMsg = 'Error in call to update_plasma_state()'
            self.services.exception(logMsg)
            raise 

        #save the output files
        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception:
            logMsg = 'Error in call to stage_output_files()'
            self.services.exception(logMsg)
            raise 

        return 0
# ------------------------------------------------------------------------------
#
# checkpoint function
# Saves plasma state files to restart directory
#
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print('RF_IC_aorsa.checkpoint() called')
        services = self.services
        services.save_restart_files(timestamp, self.RESTART_FILES)


# ------------------------------------------------------------------------------
#
# FINALIZE function
# As of now it does nothing
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print('aorsa.finalize() called')


# ------------------------------------------------------------------------------
#
# "Private"  methods
#
# ------------------------------------------------------------------------------

    def try_get_config_param(self, services, param_name, optional=False, default=None):

        try:
            value = services.get_config_param(param_name)
            print(param_name, ' = ', value)
        except Exception:
            if optional:
                print('config parameter ', param_name, ' not found')
                value = default
            else:
                message = 'required config parameter ', param_name, ' not found'
                print(message)
                services.exception(message)
                raise

        return value

    # Try to get component specific config parameter - wraps the exception handling
    def try_get_component_param(self, services, param_name, optional=False):

        if hasattr(self, param_name):
            value = getattr(self, param_name)
            print(param_name, ' = ', value)
        elif optional:
            print('optional config parameter ', param_name, ' not found')
            value = None
        else:
            message = 'required component config parameter ', param_name, ' not found'
            print(message)
            services.exception(message)
            raise

        return value
