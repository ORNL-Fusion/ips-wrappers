#! /usr/bin/env python

"""
generic_init.py  Batchelor (9-17-2018)

See working notes below

The Swiss army knife of Plasma State initializers.  It produces the intial CURRENT_STATE
and optionally the initial CURRENT_EQDSK.

This version combines several previous initializer routines and extends them.  There are
4 modes of initialization which must be specified by the config file variable INIT_MODE

INIT_MODE = touch_only
This mode only does a touch on all of the files listed as plasma state files so the
framework will have a complete set.  It does not actually put data in the plasma state file.

INIT_MODE = minimal
This is exactly the same as the previous minimal_state_init.py. It produces a CURRENT_STATE
that is empty except for some metadata:
time variables - ps%t0, ps%t1, ps%tinit, and ps%tfinal
simulation identifiers - ps%tokamak_id, ps%shot_number, ps%run_id.
ps%Global_label is set to run_id_tokamak_id_shot_number.
This data is set for all initialization modes, but for 'minimal' this is all the data
included.

INIT_MODE = existing_ps_file
This copies an existing input plasma state file and optionally an existing eqdsk file to
CURRENT_STATE and CURRENT_EQDSK.  If the config parameter GENERATE_EQDSK is set to 'True'
the CURRENT_EQDSK file is generated from equilibrium data in the INPUT_STATE_FILE.
The INPUT_STATE_FILE and INPUT_EQDSK_FILE must be specified in the config file.

INIT_MODE = mdescr
This initializes all machine description data from a plasma state machine description
file, e.g. <tokamak>.mdescr, as specified by config parameter MDESCR_FILE. In addition
if a shot configuration file config parameter, SCONFIG_FILE, is specified, the shot config
data is also loaded into CURRENT_STATE.  Machine description and shot configuration files
are namelist files that can be read and loaded using Plasma State subroutines ps_mdescr_read()
and ps_scongif_read().  Note:  machine description and shot configuration do not define
the MHD equilibrium, so the equilibrium must be specified during further component
initializations

INIT_MODE = mixed
This combines existing_ps_file and mdescr modes.  This copies an existing input plasma state
file and optionally an existing eqdsk file to CURRENT_STATE and CURRENT_EQDSK as in
existing_ps_file mode.  But initializations from MDESCR_FILE and SCONFIG_FILE are also added.
Caution is advised.  If the MDESCR_FILE or SCONFIG_FILE attempts to reallocate any of the
arrays already allocated in the CURRENT_STATE file a Plasma State error will occur.
Therefore one must provide in the simulation config file a list of which components are
to be initialized from mdescr or sconfig.  The config parameter is MDESCR_COMPONENTS.  For
example  "MDESCR_COMPONENTS = LH EC NBI"

IPS Components are

1  PLASMA   ! Thermal plasma parameters; fluid profile advance
2  EQ       ! MHD equilibrium
3  NBI      ! Neutral beam
4  IC       ! Ion cyclotron heating
5  LH       ! Lower Hybrid heating and current drive
6  EC       ! Electron cyclotron heating and current drive
7  RUNAWAY  ! Runaway electrons
8  FUS      ! Fusion product fast ions
9  RAD      ! Radiated Power & impurity transport
10  GAS     ! Neutral Gas sources & transport
11  LMHD    ! Linear MHD stability
12  RIPPLE  ! TF field ripple
13  ANOM    ! Anomalous transport

Except for possibly mode = existing_ps_file, all modes call on the fortran helper code
generic_ps_file_init.f90 to interact with the Plasma State. The fortran code is also used
in existing_ps_file mode to extract the CURRENT_EQDSK when GENERATE_EQDSK = true.

"""

# Version 6.2 (Batchelor 4/25/2020) 2to3 convert

# Version 6.2 (Batchelor 9/12/2019)
# Eliminated import get_lines, put_lines, edit_nml_file, get_global_param, and
# get_component_param.  Get these from /ips-wrappers/utilities.  Needs to be on
# PYTHON_PATH.

# Version 6.1 (Batchelor 8/26/2019)
# A new version of netCDF4 module has broken the syntax for setting variables.  The symptom
# is that the variables don't get set and the plasma state file is not updated on close().
# Fix is:
# plasma_state.variables['tokamak_id'] = tokamak ->
#       plasma_state.variables['tokamak_id'][:] = stringtoarr(tokamak, 32)
# plasma_state.variables['t0'] = t0 -> plasma_state.variables['t0'][0]
# N.B. Apparently ps.variables['t0'].assignValue(float(tinit)) also works

# Version 6.0 (Batchelor 9/17/2018)
# Implemented INIT_MODE = mixed

# Version 5.0 (Batchelor 7/29/2018)
# Eliminated all reference to NEXT_STATE

# Working notes for generic_ps_init.py 3/5/2018 (Batchelor)
# Added code to preserve the initial plasma state file generated during the init phase.
# This is done because the call to services.checkpoint_components(port_id_list, t) in the
# generic_driver.py overwrites the initial plasma state in the INIT work directory with
# the then current plasma state.
#
# Since nobody is using PRIOR_STATE and NEXT_STATE anymore I made it easier not to have
# them but maintained the capability to have them is somebody needs them.


#--------------------------------------------------------------------------
# Notes below are for the previous version, minimal_state_init.py
#--------------------------------------------------------------------------
# version 4.0 5/21/2010 (Batchelor)
#--------------------------------------------------------------------------
#
# This version supports both checkpoint and restart using framework functions
# A checkpoint() function is provided that is called by the framework
# services.checkpoint_components() fnction.  The restart process now uses framework
# function get_restart_files() in 'step'.  There is no 'restart' function for the
# INIT component.

# version 2.0 2/4/2010 (Batchelor)
#--------------------------------------------------------------------------
#
# The 'init' component produces a complete set of (almost) empty plasma state files
# and puts them in the plasma state work directory to be further populated by the
# 'init' functions of the other components.  This version only generates files that
# are called out as global config parameters in the config file.  As of now
# the files that it looks for include [CURRENT_STATE, PRIOR_STATE, NEXT_STATE,
# CURRENT_EQDSK, CURRENT_CQL, CURRENT_DQL, CURRENT_EQDSK].  We can always add more.
#
# This produces a CURRENT_STATE that is empty except for:
# time variables - ps%t0, ps%t1, ps%tinit, and ps%tfinal
# simulation identifiers - ps%tokamak_id, ps%shot_number, ps%run_id.
# ps%Global_label is set to run_id_tokamak_id_shot_number.
#
# This component drives the fortran executable generic_ps_init.f90 which uses
# Plasma State calls to generate CURRENT_STATE
#
# This version also supports RESTART as specified by the SIMULATION_MODE variable in
# the config file.  For a restart run the plasma state files are retrieved
# by the framework from the path indicated by the INPUT_DIR config parameter in the
# [generic_ps_init] section.  The new values of ps%t0 and ps%tfinal are written into
# CURRENT_STATE, and CURRENT_STATE is copied to PRIOR_STATE and NEXT_STATE if these are
# in the STATE_FILES list.  The state files are copied to the plasma state work
# directory by services.update_state().
#
# Nota Bene: For restart the plasma state files should be listed in the config file as
# input files to the generic_ps_init component.
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
import datetime
from ipsframework import Component
from netCDF4 import *
from simple_file_editing_functions import put_lines
from get_IPS_config_parameters import get_global_param, get_component_param

component_dict = {'PLASMA':10, 'EQ':2, 'NBI':9, 'IC':6, 'LH':7, 'EC':2,\
             'RUNAWAY':13, 'FUS':4, 'RAD':11, 'GAS':5, 'LMHD':8, 'RIPPLE':12, 'ANOM':1}


class generic_ps_init (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

# ------------------------------------------------------------------------------
#
# init function
#
# Does nothing.
#
# ------------------------------------------------------------------------------

    def init(self, timestamp=0.0):
        print (' ')
        print ('generic_ps_init.init() called')
        return

# ------------------------------------------------------------------------------
#
# step function
#
# Calls fortran executable init_empty_plasma_state and updates plasma state
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print (' ')
        print ('generic_ps_init.step() called')

        services = self.services

# Get timeloop for simulation
        timeloop = services.get_time_loop()
        tlist_str = ['%.3f'%t for t in timeloop]
        t = tlist_str[0]
        tinit  = tlist_str[0]
        tfinal  = tlist_str[-1]

# Check if this is a restart simulation
        simulation_mode = get_global_param(self, services, 'SIMULATION_MODE')

        if simulation_mode == 'RESTART':
            print('generic_ps_init: RESTART')
        if simulation_mode not in ['RESTART', 'NORMAL']:
            logMsg = 'generic_ps_init: unrecoginzed SIMULATION_MODE: ' + mode
            self.services.error(logMsg)
            raise ValueError(logMsg)

# ------------------------------------------------------------------------------
#
# RESTART simulation mode
#
# ------------------------------------------------------------------------------

        if simulation_mode == 'RESTART':
            # Get restart files listed in config file. Here just the plasma state files.
            restart_root = get_global_param(self, services, 'RESTART_ROOT')
            restart_time = get_global_param(self, services, 'RESTART_TIME')
            try:
                 services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
            except:
                logMsg = 'Error in call to get_restart_files()'
                self.services.exception(logMsg)
                raise

            cur_state_file = self.services.get_global_param('CURRENT_STATE')

            # Update ps%t0, ps%t1 and ps%tfinal.
            # Note ps%tinit stays the same in the plasma state file,
            # tinit from the config file timeloop is the restart time
            ps = Dataset(cur_state_file, 'r+', format = 'NETCDF3_CLASSIC')
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

            print('generic_ps_init: simulation mode NORMAL')
            nml_lines = ['&ps_init_nml\n']
            ps_file_list = get_global_param(self, services, 'STATE_FILES').split(' ')


            init_mode = get_component_param(self, services, 'INIT_MODE')
            nml_lines.append(' init_mode = \"' + init_mode + '\"\n')

        # Generate state files as dummies so framework will have a complete set
            for file in ps_file_list:
                print('touching plasma state file = ', file)
                try:
                    subprocess.call(['touch', file])
                except Exception:
                    print('No file ', file)
            if init_mode in ['touch_only', 'TOUCH_ONLY'] :
                # Update plasma state
                try:
                    services.update_state()
                except Exception as e:
                    print('Error in call to updatePlasmaState()', e)
                    raise
                return

            try:
                services.stage_input_files(self.INPUT_FILES)
            except Exception:
                message = 'generic_ps_init: Error in staging input files'
                print(message)
                services.exception(message)
                raise

            cur_state_file = get_global_param(self, services, 'CURRENT_STATE')
            cur_eqdsk_file = get_global_param(self, services, 'CURRENT_EQDSK')
            nml_lines.append(' cur_state_file = \"' + cur_state_file + '\"\n')
            nml_lines.append(' cur_eqdsk_file = \"' + cur_eqdsk_file + '\"\n')

            INPUT_EQDSK_FILE = get_component_param(self, services, 'INPUT_EQDSK_FILE', \
            optional = True)
            if (INPUT_EQDSK_FILE is None) or (len(INPUT_EQDSK_FILE) == 0):
                INPUT_EQDSK_FILE = ' '
            else:
                nml_lines.append(' input_eqdsk_file = \"' + INPUT_EQDSK_FILE + '\"\n')

                # If there is an INPUT_EQDSK_FILE copy it to CURRENT_EQDSK although
                # CURRENT_EQDSK will be overwritten with plasma state data if
                # GENERATE_EQDSK is True
                try:
                    subprocess.call(['cp', INPUT_EQDSK_FILE, cur_eqdsk_file ])
                except Exception:
                    message = 'generic_ps_init: Error copying INPUT_EQDSK_FILE to CURRENT_EQDSK'
                    print(message)
                    services.exception(message)
                    raise

# ------------------------------------------------------------------------------
            # init from existing plasma state file
            if init_mode in ['existing_ps_file', 'EXISTING_PS_FILE', 'mixed', 'MIXED'] :
                INPUT_STATE_FILE = get_component_param(self, services, 'INPUT_STATE_FILE')

                # Copy INPUT_STATE_FILE to current state file
                try:
                    subprocess.call(['cp', INPUT_STATE_FILE, cur_state_file ])
                except Exception:
                    message = 'generic_ps_init: Error in copying INPUT_STATE_FILE \
                        to current state file'
                    print(message)
                    services.exception(message)
                    raise

                # Generate cur_eqdsk_file from cur_state_file if GENERATE_EQDSK is True
                GENERATE_EQDSK = get_component_param(self, services, 'GENERATE_EQDSK', \
                optional = True)
                if GENERATE_EQDSK in ['true', 'TRUE', 'True']:
                    nml_lines.append(' generate_eqdsk = \"True\"')
                    nml_lines.append('/\n')
                    put_lines('generic_ps_init.nml', nml_lines)

                    init_bin = os.path.join(self.BIN_PATH, 'generic_ps_init')
                    print('Executing ', init_bin)
                    retcode = subprocess.call(init_bin)
                    if (retcode != 0):
                       print('Error executing ', init_bin)
                       raise

             # Copy INPUT_EQDSK_FILE, if there is one, to cur_eqdsk_file.
             # Nota Bene: If there is an INPUT_EQDSK_FILE specified in config this copy
             # will overwrite any eqdsk generated in the block above.
                if INPUT_EQDSK_FILE != ' ':
                    try:
                        subprocess.call(['cp', INPUT_EQDSK_FILE, cur_eqdsk_file ])
                    except Exception as e:
                        message =  'generic_ps_init: Error in copying input_eqdsk_file'
                        print(message)
                        services.exception(message)
                        raise e

# ------------------------------------------------------------------------------
            # init from machine description file and possibly sconfig file
            if init_mode in ['mdescr', 'MDESCR', 'mixed', 'MIXED'] :
                MDESCR_FILE = get_component_param(self, services, 'MDESCR_FILE')
                nml_lines.append(' mdescr_file = \"' + MDESCR_FILE + '\"\n')
                SCONFIG_FILE = get_component_param(self, services, 'SCONFIG_FILE', \
                optional = 'TRUE')

                if (SCONFIG_FILE is None) or (len(SCONFIG_FILE) == 0):
                   SCONFIG_FILE = ' '
                else:
                    nml_lines.append(' sconfig_file = \"' + SCONFIG_FILE + '\"\n')

                INPUT_EQDSK_FILE = get_component_param(self, services, 'INPUT_EQDSK_FILE', \
                optional = True)
                if (INPUT_EQDSK_FILE is None) or (len(INPUT_EQDSK_FILE) == 0):
                   INPUT_EQDSK_FILE = ' '
                else:
                   nml_lines.append(' input_eqdsk_file = \"' + INPUT_EQDSK_FILE + '\"\n')

# ------------------------------------------------------------------------------
			# For init_mode = mixed add input_state_file to namelist
            if init_mode in ['mixed', 'MIXED'] :
                nml_lines.append(' input_state_file = \"' + INPUT_STATE_FILE + '\"\n')

				# Retrieve list of IPS components which are to be initialized from
				# mdescr/sconfig and construct cclist for generic_ps_init.f90
                mdescr_components =  get_component_param(self, services, 'MDESCR_COMPONENTS')
                if isinstance(mdescr_components, type('str')):
                	mdescr_components = [mdescr_components]
                cclist = [0 for i in range(len(component_dict))]
                for comp in mdescr_components:
                    if comp not in list(component_dict.keys()):
                        message = 'generic_ps_init: Unknown IPS component ' + comp
                        print(message)
                        services.exception(message)
                        raise
                    cclist[component_dict[comp]-1] = 1
                print('cclist = ', cclist)
                cclist_string = ''
                for i in range(len(cclist)):
                    cclist_string = cclist_string + str(cclist[i]) + ', '
                nml_lines.append(' cclist = ' + cclist_string + '\n')

# ------------------------------------------------------------------------------
            # For 'minimal', 'mdescr' and 'mixed' modes generate namelist for the fortran
            # helper code generic_ps_init.f90 and execute it
            if init_mode in ['minimal', 'MINIMAL', 'mdescr', 'MDESCR', 'mixed', 'MIXED'] :
                nml_lines.append('/\n')
                print('nml_lines = ', nml_lines)
                put_lines('generic_ps_init.nml', nml_lines)
                print('Did put_lines')

                init_bin = os.path.join(self.BIN_PATH, 'generic_ps_init')
                print('Executing ', init_bin)
                retcode = subprocess.call(init_bin)
                if (retcode != 0):
                   print('Error executing ', init_bin)
                   raise

# ------------------------------------------------------------------------------
            # For all init init modes insert run identifiers and time data
            # (do it here in python instead of in minimal_state_init.f90 as before)
            # For minimal mode this is the only data in initial state
            tokamak = get_global_param(self, services, 'TOKAMAK_ID')
            shot_number = get_global_param(self, services, 'SHOT_NUMBER')
            run_id = get_global_param(self, services, 'RUN_ID')

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
            print('Error in call to updatePlasmaState()', e)
            raise

# "Archive" output files in history directory
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

# ------------------------------------------------------------------------------
#
# checkpoint function
#
# Saves plasma state files to restart directory
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print('generic_ps_init.checkpoint() called')

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
        print('generic_ps_init.finalize() called')
