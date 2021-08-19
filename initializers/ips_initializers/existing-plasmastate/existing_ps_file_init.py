#! /usr/bin/env python

# version 1.0 1/25/2016 (Batchelor)
#--------------------------------------------------------------------------
#
# This init component generates an initial plasma state file from an existing plasma 
# state file.  It takes the existing plasma state file as an input file in the
# simulation config file, copies it to current plasma state file and updates the plasma
# state.  The input state file must be specified in the [existing_ps_file_init] section
# of the config file as INPUT_STATE_FILE = <file name>.
#
# There should also be an INPUT_EQDSK_FILE = <file name> line in the config file
# which also is copied to the current eqdsk if the <file name> isn't blank.
#
# CAVEAT UTILITOR
# The plasma state files produced by this code will have have complete arrays allocated 
# for the components present when the input state file was generated.  Normally the 'init'
# functions of the physics component scripts will try to allocate the arrays for data 
# supplied by that component However arrays in Plasma State files can only be allocated 
# once. Therefore when using the pre-allocated state from the present code either:
#
# 1) The physics init should not re-allocate the arrays and should ensure that the array
#    dimensions of its initial data match those in the input state.
# or
# 2) The physics component init should instantiate a new plasma state structure, copy
#    the data from all components except itself to the new instance, should allocate its
#    own arrays in the new state with whatever dimensions it likes, and fill the initial
#    data.  Then the new state instance should be written out as the current plasma state.
#
# This version also supports RESTART as specified by the SIMULATION_MODE variable in
# the config file.  For a restart run the plasma state files are retrieved 
# by the framework from the path indicated by the INPUT_DIR config parameter in the 
# [existing_ps_file_init] section.  The new values of ps%t0 and ps%tfinal are written into
# CURRENT_STATE, and CURRENT_STATE is copied to PRIOR_STATE and NEXT_STATE if these are
# in the STATE_FILES list.  The state files are copied to the plasma state work
# directory by services.update_state().
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

# Working notes
#
# 1/25/2016 Changed old Scientific.IO.NetCDF format to netCDF4

import sys
import os
import subprocess
import getopt
import shutil
import string
from ipsframework import Component
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
            print 'Error in existing_ps_file_init: step () : No services'
            raise Exception('Error in existing_ps_file_init: step (): No services')
        services = self.services

# Get timeloop for simulation
        timeloop = services.get_time_loop()
        tlist_str = ['%.3f'%t for t in timeloop]
        t = tlist_str[0]
        tinit  = tlist_str[0]
        tfinal  = tlist_str[-1]
        
# Check if this is a restart simulation

        try:
            mode = services.get_config_param('SIMULATION_MODE')
            if mode == 'RESTART':
                print 'existing_ps_file_init: RESTART'
            if mode not in ['RESTART', 'NORMAL']:
                print 'existing_ps_file_init: unrecoginzed SIMULATION_MODE: ', mode
                print message
                services.exception(message)
                raise
        except Exception, e:
            print 'existing_ps_file_init: No SIMULATION_MODE variable in config file' \
                  ', NORMAL assumed', e
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
            except Exception, e:
                message = 'existing_ps_file_init: error getting restart config parameters '
                print message
                services.exception(message)
                raise e
            
            cur_state_file = self.services.get_config_param('CURRENT_STATE')
    
            # Update ps%t0, ps%t1 and ps%tfinal. Note ps%tinit stays the same in the plasma
            # state file, but this value of tinit is the restart time from the config file

            ps = Dataset(cur_state_file, 'r+', format = 'NETCDF3_CLASSIC')
            ps.variables['t0'].assignValue(float(tinit))
            ps.variables['t1'].assignValue(float(tinit))
            ps.variables['tfinal'].assignValue(float(tfinal))
            ps.close()
            
            print 'existing_ps_file_init RESTART: tinit  = ', tinit, '  tfinal = ', tfinal

        
# ------------------------------------------------------------------------------
#
# NORMAL simulation mode
#
# ------------------------------------------------------------------------------
        
        else:
            print 'existing_ps_file_init: simulation mode NORMAL'
            ps_file_list = services.get_config_param('STATE_FILES').split(' ')
            print 'ps_file_list = ', ps_file_list

        # Get name of current plasma state file
            try:
                cur_state_file = self.services.get_config_param('CURRENT_STATE')
            except Exception, e:
                message = 'existing_ps_file_init: error getting CURRENT_STATE_FILE name'
                print message
                services.exception(message)
                raise e

         # Get name of current eqdsk file
            try:
                cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            except:
                message = 'existing_ps_file_init: error getting CURRENT_EQDSK_FILE name'
                print message
                services.exception(message)
                raise
       
        # Get component specific configuration parameters.
            try:
                INPUT_FILES = self.INPUT_FILES
                OUTPUT_FILES = self.OUTPUT_FILES
                RESTART_FILES = self.RESTART_FILES
                INPUT_STATE_FILE = self.INPUT_STATE_FILE
                INPUT_EQDSK_FILE = self.INPUT_EQDSK_FILE
            except Exception, e:
                message = 'existing_ps_file_init: error getting config parameters'
                print message
                services.exception(message)
                raise e
        
        # Get input files  
            try:
              INPUT_FILES = self.INPUT_FILES
              services.stage_input_files(INPUT_FILES)
            except Exception, e:
                message = 'existing_ps_file_init: Error in call to stage_input_files()'
                print message
                services.exception(message)
                raise e
      
         # Copy INPUT_STATE_FILE to current state file
            try:
                subprocess.call(['cp', INPUT_STATE_FILE, cur_state_file ])
            except Exception, e:
                message = 'existing_ps_file_init: Error in copying INPUT_STATE_FILE to current state file'
                print message
                services.exception(message)
                raise e

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

         # Copy INPUT_EQDSK_FILE to cur_eqdsk_file if there is one
            if INPUT_EQDSK_FILE != '':
                try:
                    subprocess.call(['cp', INPUT_EQDSK_FILE, cur_eqdsk_file ])
                except Exception, e:
                    message =  'existing_ps_file_init: Error in renaming INPUT_EQDSK_FILE' 
                    print message
                    services.exception(message)
                    raise e

        # Generate other state files as dummies so framework will have a complete set
            for file in ps_file_list:
                print 'touching plasma state file = ', file
                try:
                    subprocess.call(['touch', file])
                except Exception, e:
                    print 'No file ', file

# Update plasma state files in plasma_state work directory
        try:
          services.update_state()
        except Exception, e:
            message = 'existing_ps_file_init: Error in call to update_state()'
            print message
            services.exception(message)
            raise e

# "Archive" output files in history directory
        try:
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
            message = 'existing_ps_file_init: Error in call to stage_output_files()', e
            print message
            services.exception(message)
            raise e

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
        services.stage_state()
        services.save_restart_files(timestamp, self.RESTART_FILES)
        

# ------------------------------------------------------------------------------
#
# finalize function
#
# Does nothing
# ------------------------------------------------------------------------------



    def finalize(self, timestamp=0.0):
        print 'existing_ps_file_init.finalize() called'

