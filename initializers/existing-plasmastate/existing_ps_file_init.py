#! /usr/bin/env python

# version 0.0 9/20/2015 (Batchelor)
#--------------------------------------------------------------------------
#
# This init component generates an initial plasma state file from an existing plasma 
# state.  It takes the existing plasma state file as an input file in the
# simulation config file, copies it to current plasma state file and updates the plasma
# state.
#
# This version also supports RESTART as specified by the SIMULATION_MODE variable in
# the config file.  For a restart run the plasma state files are retrieved 
# by the framework from the path indicated by the INPUT_DIR config parameter in the 
# [existing_ps_file_init] section.  The new values of ps%t0 and ps%tfinal are written into
# CURRENT_STATE, and CURRENT_STATE is copied to PRIOR_STATE and NEXT_STATE if these are
# in the PLASMA_STATE_FILES list.  The state files are copied to the plasma state work
# directory by services.update_plasma_state().
#
# Nota Bene: For restart the plasma state files should be listed in the config file as  
# input files to the existing_ps_file_init component.
#
# N.B. The other plasma state files that in previous versions were produced by the
#      fortran code are now produced here. These files include:
#      prior_state file and next_state file as well as the dummy files: cur_cql_file
#      cur_eqdsk_file, cur_dql_file, and cur_jsdsk_file.
#
# N.B. Both ps%t0 and ps%t1 are set to the value time_stamp.  tinit and tfinal
#      are generated here from the TIME_LOOP variable in the
#      simulation config file.  Note that the initial t0 can be different from 
#      tinit, as might be needed in a restart.
#
# ------------------------------------------------------------------------------

import sys
import os
import subprocess
import getopt
import shutil
import string
from  component import Component
from netCDF4 import *

class existing_ps_file_init (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

# ------------------------------------------------------------------------------
#
# init function
#
# Does nothing.
#
# ------------------------------------------------------------------------------

    def init(self, timestamp=0.0):
        print (' ')
        print ('existing_ps_file_init.init() called')
        return

# ------------------------------------------------------------------------------
#
# step function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print (' ')
        print ('existing_ps_file_init.step() called')

        if (self.services == None) :
            logMsg = 'Error in existing_ps_file_init: step () : No services'
            self.services.error(logMsg)
            raise Exception(logMsg)

        services = self.services

# Check if this is a restart simulation

        try:
            mode = services.get_config_param('SIMULATION_MODE')
            if mode == 'RESTART':
                print 'existing_ps_file_init: RESTART'
            if mode not in ['RESTART', 'NORMAL']:
                logMsg = 'existing_ps_file_init: unrecoginzed SIMULATION_MODE: '+mode
                services.error(logMsg)
                raise Exception(logMsg)

        except:
            logMsg = 'existing_ps_file_init: No SIMULATION_MODE variable in config file, NORMAL assumed'
            self.services.exception(logMsg)
            raise
# ------------------------------------------------------------------------------
#
# RESTART simulation mode
#
# ------------------------------------------------------------------------------
            
        if mode == 'RESTART':
            # Get restart files listed in config file. Here just the plasma state files.
            try:
                restart_root = services.get_config_param('RESTART_ROOT')
                restart_time = services.get_config_param('RESTART_TIME')
                services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
            except:
                message = 'existing_ps_file_init: error getting restart config parameters '
                services.exception(message)
                raise
            
            cur_state_file = self.services.get_config_param('CURRENT_STATE')
    
            # Update ps%t0, ps%t1 and ps%tfinal. Note ps%tinit stays the same in the plasma
            # state file, but this value of tinit is the restart time from the config file
            ps = NetCDFFile(cur_state_file, 'r+')
            ps.variables['t0'].assignValue(float(tinit))
            ps.variables['t1'].assignValue(float(tinit))
            ps.variables['tfinal'].assignValue(float(tfinal))
            ps.close()
        
# ------------------------------------------------------------------------------
#
# NORMAL simulation mode
#
# ------------------------------------------------------------------------------
        
        else:

            ps_file_list = services.get_config_param('PLASMA_STATE_FILES').split(' ')
            print 'ps_file_list = ', ps_file_list

        # Get name of current plasma state file
            try:
                cur_state_file = self.services.get_config_param('CURRENT_STATE')
            except:
                message = 'existing_ps_file_init: error getting CURRENT_STATE_FILE name'
                services.exception(message)
                raise

         # Get name of current eqdsk file
            try:
                cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            except:
                message = 'existing_ps_file_init: error getting CURRENT_EQDSK_FILE name'
                services.exception(message)
                raise
       
        # Get component specific configuration parameters.
            try:
                INPUT_FILES = self.INPUT_FILES
                OUTPUT_FILES = self.OUTPUT_FILES
                RESTART_FILES = self.RESTART_FILES
                input_state_file = self.INPUT_STATE_FILE
                input_eqdsk_file = self.INPUT_EQDSK_FILE
            except:
                message = 'existing_ps_file_init: error getting config parameters'
                services.exception(message)
                raise
        
        # Get input files  
            try:
              INPUT_FILES = self.INPUT_FILES
              services.stage_input_files(INPUT_FILES)
            except:
                message = 'existing_ps_file_init: Error in call to stageInputFiles()'
                services.exception(message)
                raise
      
         # Copy input_state_file to current state file
            try:
                subprocess.call(['cp', input_state_file, cur_state_file ])
            except:
                message = 'existing_ps_file_init: Error in copying input_state_file to current state file'
                services.exception(message)
                raise

        #Open Plasma State file and put in t0, t1, tinit, and tfinal from config file
            plasma_state = Dataset(cur_state_file, 'r', format = 'NETCDF3_CLASSIC')
            timeloop = services.get_time_loop()
            t0 = timeloop[0]
            t1 = t0
            tfinal = timeloop[-1]
            plasma_state.variables['t0'] = t0
            plasma_state.variables['t1'] = t1
            plasma_state.variables['tinit'] = t0
            plasma_state.variables['tfinal'] = tfinal
            plasma_state.close()

         # Copy input_eqdsk_file to cur_eqdsk_file if there is one
            if input_eqdsk_file != '':
                try:
                    subprocess.call(['cp', input_eqdsk_file, cur_eqdsk_file ])
                except:
                    message =  'existing_ps_file_init: Error in renaming input_eqdsk_file' 
                    services.exception(message)
                    raise

          # Generate other state files as dummies so framework will have a complete set
#             for file in ps_file_list:
#                 print 'file = ', file
#                 subprocess.call(['touch', file])                

            try:
                cur_bc_file = services.get_config_param('CURRENT_BC')
                subprocess.call(['touch', cur_bc_file])
            except:
                logMsg = 'No CURRENT_BC file '
                self.services.info(logMsg)
                pass

            try:
                cur_cql_file = services.get_config_param('CURRENT_CQL')
                subprocess.call(['touch', cur_cql_file])
            except:
                logMsg = 'No CURRENT_CQL file '
                self.services.info(logMsg)
                pass
    
            try:
                cur_dql_file = services.get_config_param('CURRENT_DQL')
                subprocess.call(['touch', cur_dql_file])
            except:
                logMsg = 'No CURRENT_DQL file '
                self.services.info(logMsg)
                pass
    
            try:
                cur_jsdsk_file = services.get_config_param('CURRENT_JSDSK')
                subprocess.call(['touch', cur_jsdsk_file])
            except:
                logMsg = 'No CURRENT_JSDSK file '
                pass

        # Copy current plasma state to prior state and next state if they are in state
            try:
                prior_state_file = services.get_config_param('PRIOR_STATE')
                shutil.copyfile(cur_state_file, prior_state_file)
            except:
                logMsg = 'No PRIOR_STATE file '
                self.services.info(logMsg)
                pass

            try:
                next_state_file = services.get_config_param('NEXT_STATE')
                shutil.copyfile(cur_state_file, next_state_file)
            except:
                logMsg = 'No NEXT_STATE file '
                self.services.info(logMsg)
                pass

        # Update plasma state files in plasma_state work directory
            try:
              services.update_plasma_state()
            except:
                message = 'existing_ps_file_init: Error in call to update_plasma_state()'
                self.services.exception(message)
                raise

    # "Archive" output files in history directory
            try:
              services.stage_output_files(timeStamp, self.OUTPUT_FILES)
            except:
                message = 'existing_ps_file_init: Error in call to stage_output_files()', e
                self.services.exception(message)
                raise

            # For benefit of framework file handling generate dummy dakota.out file
            subprocess.call(['touch', 'dakota.out'])

# ------------------------------------------------------------------------------
#
# checkpoint function
#
# Saves plasma state files to restart directory
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print 'existing_ps_file_init.checkpoint() called'
        
        services = self.services
        services.stage_plasma_state()
        services.save_restart_files(timestamp, self.RESTART_FILES)
        

# ------------------------------------------------------------------------------
#
# finalize function
#
# Does nothing
# ------------------------------------------------------------------------------



    def finalize(self, timestamp=0.0):
        print 'existing_ps_file_init.finalize() called'
