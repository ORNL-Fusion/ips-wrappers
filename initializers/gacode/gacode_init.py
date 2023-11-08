#! /usr/bin/env python

# -------------------------------------------------------------------------------
#
# Initializes the plasma state from a GACODE file typically used for running
# TGYRO, CGYRO, etc.
#
# -------------------------------------------------------------------------------

import sys
import os
import subprocess
import getopt
import shutil
import string
import datetime
from ipsframework import Component
from netCDF4 import *
from simple_file_editing_functions import put_lines
from get_IPS_config_parameters import get_global_param, get_component_param

class gacode_init (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    # ---------------------------------------------------------------------------
    # init function
    #
    # does nothing
    # ---------------------------------------------------------------------------
    def init(self, timestamp=0.0):
        print (' ')
        print ('gacode_init.init() called')
        return

    # ---------------------------------------------------------------------------
    # step function
    #
    # Calls fortran executable init_gacode and reads in a GACODE file then
    # maps the GACODE file values into the plasma state file
    # ---------------------------------------------------------------------------
    def step(self, timestamp):
        print(' ')
        print('gacode_init.step() called')

        services = self.services
        
        # Get timeloop for simulation
        timeloop = services.get_time_loop()
        tlist_str = ['%.3f'%t for t in timeloop]
        t = tlist_str[0]
        tinit  = tlist_str[0]
        tfinal  = tlist_str[-1]

        # Check if this is a restart simulation
        try:
            simulation_mode = services.get_config_param('SIMULATION_MODE')
        except:
            logMsg = 'gacode_init: No SIMULATION_MODE variable in config \
                      file. Please set NORMAL or RESTART'
            self.services.exception(logMsg)
            raise
        
        if simulation_mode == 'RESTART':
            print('gacode_init: RESTART')
        if simulation_mode not in ['RESTART', 'NORMAL']:
            logMsg = 'gacode_init: unrecoginzed SIMULATION_MODE: ' + simulation_mode
            self.services.error(logMsg)
            raise ValueError(logMsg)

        # Restart simulation mode
        if simulation_mode == 'RESTART':
            # Get restart files listed in config file. Here just the plasma state files.
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            try:
                 services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
            except:
                logMsg = 'Error in call to get_restart_files()'
                self.services.exception(logMsg)
                raise

            cur_state_file = services.get_config_param('CURRENT_STATE')

            # Update ps%t0, ps%t1 and ps%tfinal.
            # Note ps%tinit stays the same in the plasma state file,
            # tinit from the config file timeloop is the restart time
            ps = Dataset(cur_state_file, 'r+', format = 'NETCDF3_CLASSIC')
            ps.variables['t0'].assignValue(float(tinit))
            ps.variables['t1'].assignValue(float(tinit))
            ps.variables['tfinal'].assignValue(float(tfinal))
            ps.close()

        else:
            print('gacode_init: simulation mode NORMAL')
            
            ps_file_list = services.get_config_param('PLASMA_STATE_FILES').split(' ')

            #Required params
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            cur_gacode_file = services.get_config_param('CURRENT_GACODE')

            #Optional params
            ic_rfpwr = self.try_get_config_param(services, 'RFPWR_IC',optional=True, default = '0.0')
            ic_addmin = self.try_get_config_param(services, 'ICRF_ADDMIN',optional=True, default = 'None')
            ic_minfrac = self.try_get_config_param(services, 'ICRF_MINFRAC',optional=True, default = '0.0')

            ic_min_is_thermal_arg = self.try_get_config_param(services, 'MIN_IS_THERMAL',optional=True, default = 'F')

            ic_antspec = self.try_get_config_param(services, 'ICRF_SPECFILE',optional=True)
            
            #Component params
            
            if ic_antspec != None:
                antspec = 'T'
                antspec_file = ic_antspec
            else:
                antspec = 'F'
                antspec_file = 'None'

            if ic_min_is_thermal_arg in ['T','True','TRUE','true']:
                ic_min_is_thermal = 'T'
            else:
                ic_min_is_thermal = 'F'

            if ic_addmin != 'None':
                if ic_addmin.strip() in ['He3','HE3','he3','he-3','He-3','HE-3']:
                    ic_addmin = 'He3'
                elif ic_addmin.strip() in ['H','h']:
                    ic_addmin = 'H'
                
            # Generate state files as dummies so framework has complete set
            for file in ps_file_list:
                print('touching plasma state file = ', file)
                try:
                    subprocess.call(['touch', file])
                except Exception:
                    print('No file ', file)

            try:
                services.stage_input_files(self.INPUT_FILES)
            except Exception:
                message = 'generic_ps_init: Error in staging input files'
                print(message)
                services.exception(message)
                raise

            # Copy eqdsk files to the plasma state eqdsk
            INPUT_EQDSK_FILE = get_component_param(self, services,  \
            'INPUT_EQDSK_FILE')
            try:
                subprocess.call(['cp', INPUT_EQDSK_FILE, cur_eqdsk_file ])
            except Exception:
                message = 'generic_ps_init: Error copying INPUT_EQDSK_FILE to CURRENT_EQDSK'
                print(message)
                services.exception(message)
                raise

            # Copy gacode files to the plasma state
            INPUT_GACODE_FILE = get_component_param(self, services, \
             'INPUT_GACODE_FILE')
            try:
                subprocess.call(['cp', INPUT_GACODE_FILE, cur_gacode_file ])
            except Exception:
                message = 'generic_ps_init: Error copying INPUT_GACODE_FILE to CURRENT_GACODE'
                print(message)
                services.exception(message)
                raise
            
            # Make namelist file for script
            nml_lines = ['&ps_init_nml\n']
            nml_lines.append(' cur_state_file = \"' + cur_state_file + '\",\n')
            nml_lines.append(' cur_eqdsk_file = \"' + cur_eqdsk_file + '\",\n')
            nml_lines.append(' cur_gacode_file = \"' + cur_gacode_file + '\",\n')
            nml_lines.append(' icrfpwr = ' + ic_rfpwr + ',\n')
            nml_lines.append(' addmin = "' + ic_addmin +'",\n')
            nml_lines.append(' minfrac = ' + ic_minfrac + ',\n')
            nml_lines.append(' ministhermal = ' + ic_min_is_thermal + ',\n')
            nml_lines.append(' antspec = ' + antspec + ',\n')
            nml_lines.append(' antspec_file = \"' + antspec_file + '\",\n')
            nml_lines.append('/\n')
            put_lines('gacode_init.nml', nml_lines)

            # Run initializer script to make plasma state file
            init_bin = os.path.join(self.BIN_PATH, 'gacode_init')
            init_log = open('log_gacode_init','w')
            print('Executing ', init_bin)
            retcode = subprocess.call(init_bin, stdout = init_log, \
                  stderr=subprocess.STDOUT)
            if (retcode != 0):
                print('Error executing ', init_bin)
                services.error('Error executing',init_bin)
                raise Exception('Error executing',init_bin)

            # Add run identifiers and time data to initialized file
            tokamak = services.get_config_param('TOKAMAK_ID')
            shot_number = services.get_config_param('SHOT_NUMBER')
            run_id = services.get_config_param('RUN_ID')

            timeloop = services.get_time_loop()
            t0 = timeloop[0]
            t1 = t0
            tfinal = timeloop[-1]

            # Put into current plasma state
            plasma_state = Dataset(cur_state_file, 'r+', format = 'NETCDF3_CLASSIC')
            plasma_state.variables['tokamak_id'][:] = stringtoarr(tokamak, 32)
            plasma_state.variables['RunID'][:] = stringtoarr(run_id, 32)
            plasma_state.variables['shot_number'][0] = shot_number
            plasma_state.variables['t0'][0] = t0
            plasma_state.variables['t1'][0] = t1
            plasma_state.variables['tinit'][0] = t0
            plasma_state.variables['tfinal'][0] = tfinal
            plasma_state.close()

        # Preserve initial plasma state file
        try:
            shutil.copyfile(cur_state_file, 'initial_PLASMA_STATE.nc')
        except Exception as e:
            print('Copy to initial_PLASMA_STATE file failed ', e)

        # For benefit of framework file handling generate dummy dakota.out file
        subprocess.call(['touch', 'dakota.out'])
        
            
        # Update plasma state
        try:
            services.update_state()
        except Exception as e:
            print('Error in call to update_state()', e)
            raise

        # Archive output files
        services.stage_output_files(timestamp, self.OUTPUT_FILES)


    def checkpoint(self, timestamp=0.0):
        print('gacode_init.checkpoint() called')

        services = self.services
        services.stage_state()
        services.save_restart_files(timestamp, self.RESTART_FILES)

    def finalize(self, timestamp=0.0):
        print('gacode_init.finalize() called')


    # "Private" methods (wrappers for IPS services)
    #---------------------------------------------------------------------------------
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
