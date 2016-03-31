#! /usr/bin/env python

# version 0.0 9/11/2015 (Batchelor)
#--------------------------------------------------------------------------
# This is a generalization of the minimal_state_init component in that it
# initializes the plasma state machine description data from an mdescr namelist
# file using the plasma state function ps_mdescr_read().  The actual initialization is
# done by the fortran helper code mdescr_state_init.f90
#--------------------------------------------------------------------------
#
# The 'init' component produces a complete set of (almost) empty plasma state files
# and puts them in the plasma state work directory to be further populated by the
# 'init' functions of the other components.  This version only generates files that
# are called out as global config parameters in the config file.  As of now
# the files that it looks for include [CURRENT_STATE, PRIOR_STATE, NEXT_STATE,
# CURRENT_EQDSK, CURRENT_CQL, CURRENT_DQL, CURRENT_EQDSK].  We can always add more.
#
# This produces a CURRENT_STATE that in addition to the machine description data has:
# time variables - ps%t0, ps%t1, ps%tinit, and ps%tfinal
# simulation identifiers - ps%tokamak_id, ps%shot_number, ps%run_id.
# ps%Global_label is set to run_id_tokamak_id_shot_number.
#
# This component drives the fortran executable mdescr_state_init.f90 which uses
# Plasma State calls to generate CURRENT_STATE
#
# This version also supports RESTART as specified by the SIMULATION_MODE variable in
# the config file.  For a restart run the plasma state files are retrieved
# by the framework from the path indicated by the INPUT_DIR config parameter in the
# [mdescr_state_init] section.  The new values of ps%t0 and ps%tfinal are written into
# CURRENT_STATE, and CURRENT_STATE is copied to PRIOR_STATE and NEXT_STATE if these are
# in the PLASMA_STATE_FILES list.  The state files are copied to the plasma state work
# directory by services.update_plasma_state().
#
# Nota Bene: For restart the plasma state files should be listed in the config file as
# input files to the mdescr_state_init component.
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
#from Scientific.IO.NetCDF import *
#import Numeric

class mdescr_state_init (Component):
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
        print ('mdescr_state_init.init() called')
        return

# ------------------------------------------------------------------------------
#
# step function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print (' ')
        print ('mdescr_state_init.step() called')

        if (self.services == None) :
            print 'Error in mdescr_state_init: step () : No services'
            raise Exception('Error in mdescr_state_init: step (): No services')
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
                print 'mdescr_state_init: RESTART'
            if mode not in ['RESTART', 'NORMAL']:
                print 'mdescr_state_init: unrecoginzed SIMULATION_MODE: ', mode
                raise
                return
        except Exception, e:
            print 'mdescr_state_init: No SIMULATION_MODE variable in config file' \
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
                print 'Error in call to get_restsrt_files()' , e
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

            # Get run identifiers
            shot_number = services.get_config_param('SHOT_NUMBER')
            run_id = services.get_config_param('RUN_ID')

            mdescr_file = self.MDESCR_FILE
            print 'mdescr_file =', mdescr_file, '  shot_number = ', shot_number,
            '   run_id = ', run_id

            # Generate initial plasma state files
            services.stage_input_files(self.INPUT_FILES)
            cur_state_file = services.get_config_param('CURRENT_STATE')

            init_bin = os.path.join(self.BIN_PATH, 'mdescr_state_init')
    
            print 'Executing ', [init_bin, cur_state_file]
            retcode = subprocess.call([init_bin, cur_state_file,
                mdescr_file, shot_number, run_id, tinit, tfinal, t])
            if (retcode != 0):
               print 'Error executing ', init_bin
               raise
               
            # Put in any shot configuration data the components need in order
            # to initialize themselves.

# The only component now that needs this is NUBEAM.  It expects TSC to have specified the
# beam species labels (H,D, T ...).  What is here is a quick and dirty way to get the
# 2 ITER beams set to 'D'.  This should be put into the NUBEAM section of the config file.
#             ps = NetCDFFile(cur_state_file, 'r+')
#             ps.variables['ps%nbion'][0] = 'D'
#             ps.variables['ps%nbion'][1] = 'D'
#             ps.close()

            # Generate other state files
            try:
                cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
                subprocess.call(['touch', cur_eqdsk_file])
            except Exception, e:
                print 'No CURRENT_EQDSK file ', e

            try:
                cur_cql_file = services.get_config_param('CURRENT_CQL')
                subprocess.call(['touch', cur_cql_file])
            except Exception, e:
                print 'No CURRENT_CQL file ', e

            try:
                cur_dql_file = services.get_config_param('CURRENT_DQL')
                subprocess.call(['touch', cur_dql_file])
            except Exception, e:
                print 'No CURRENT_DQL file ', e

            try:
                cur_jsdsk_file = services.get_config_param('CURRENT_JSDSK')
                subprocess.call(['touch', cur_jsdsk_file])
            except Exception, e:
                print 'No CURRENT_JSDSK file ', e

        # For benefit of framework file handling generate dummy dakota.out file
        subprocess.call(['touch', 'dakota.out'])

        # Copy current plasma state to prior state and next state
        try:
            prior_state_file = services.get_config_param('PRIOR_STATE')
            shutil.copyfile(cur_state_file, prior_state_file)
        except Exception, e:
            print 'No PRIOR_STATE file ', e

        try:
            next_state_file = services.get_config_param('NEXT_STATE')
            shutil.copyfile(cur_state_file, next_state_file)
        except Exception, e:
            print 'No NEXT_STATE file ', e


# Update plasma state
        try:
            services.update_plasma_state()
        except Exception, e:
            print 'Error in call to updatePlasmaState()', e
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
        print 'mdescr_state_init.checkpoint() called'

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
        print 'mdescr_state_init.finalize() called'
