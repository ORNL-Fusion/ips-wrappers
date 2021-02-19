#! /usr/bin/env python

"""
MCMD version of TORIC component.  Slightly adapted from rf_ic_toric.py

"""
# Working notes:  DBB 2/9/2011
# Eliminated 'torica.inp' as an input file.  The toric code requires a 'torica.inp'. However
# in the IPS this file is written either by the 'do_toric_init_abr.f90' during 'init' or by
# the 'prepare_toric_input_abr.f90' code during 'step'  Therefore I'm eliminating 'torica.inp'as
# an input file in the config files and eliminating reference to it in this component.
#
# Also put in the capability to define 'INPUT_SUFFIX' in the config file. If present, this
# component copies the input file machine.inp_<INPUT_SUFFIX> to generic file
# 'name machine.inp'. Most of the components now use this for flexibility in maintaining
# input files.
#
# Working notes:  DBB 5/24/2010
# Added checkpoint function.  Saves RESTART_FILES as listed in config file using framework
# services.save_restart_files()
#
# Working notes:  DBB 4/28/2010
# Added restart function
#
# Going back and adding exception handling
#
#      Note: the framework doesn't require you to define global parameters <CURRENT_STATE>,
#      <CURRENT_EQDSK>, etc in the config file, it's just a convenience.  But this
#      component uses these parmeters and will break if you don't define them there.
#
#      Note: This component gets the AORSA executable file name and path from AORSA config
#      parameter <AORSA_BIN>, so this must be defined in the AORSA section of the config file.
#
#      About zero ICRF power:  We don't want to run and ICRF code when the RF power is zero.
#      In the past we set the icrf source profiles to zero in this script and did
#      ps_store_plasma_state.  However with MCMD the STEP function needs to produce a partial
#      only plasma state containing the ICRF data (i.e. ps_write_update_file not
#      ps_store_plasma_state).  This really has to done with plasma state fortran code not not
#      with the python netcdf interface.So I wrote a simple code called zero_RF_IC_power to set
#      all ICRF source profiles in plasma state to zero and then write a partial plasma state.
#      For now this code lives in: /ips/trunk//components/rf/model_RF_IC and it gets built and
#      installed by the Makefile there.

import sys
import os
import subprocess
import getopt
import shutil
import string
from ipsframework import Component
from netCDF4 import *
from get_IPS_config_parameters import get_global_param, get_component_param

class toric (Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print(('Created %s' % (self.__class__)))

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------


    def init(self, timeStamp=0):
        print('toric.init() called')

        services = self.services
        workdir = services.get_working_dir()

    # Get global configuration parameters
        cur_state_file = get_global_param(self, services, 'CURRENT_STATE')
        cur_eqdsk_file = get_global_param(self, services, 'CURRENT_EQDSK')

        try:
             self.toric_log = os.path.join(workdir, 'log.toric')
        except:
            logMsg = 'rf_ic_toric_mcmd: error in getting log.toric path'
            self.services.exception(logMsg)
            raise 
        toric_log = self.toric_log


      # Copy plasma state files over to working directory
        try:
            services.stage_state()
        except Exception:
            logMsg = 'Error in call to stage_state()'
            self.services.exception(logMsg)
            raise 

      # Get input files
        try:
            services.stage_input_files(self.INPUT_FILES)
        except:
            logMsg = 'Error in call to stage_input_files()'
            self.services.exception(logMsg)
            raise

      # Copy machine.inp_<suffix> to generic file name -> machine.inp if there is
      # a suffix
        try:
            suffix = self.INPUT_SUFFIX
            have_suffix = True
        # If suffix is not empty put an underscore in front of it.
            if len(suffix) > 0:
                print(('INPUT_SUFFIX = ', suffix))
                suffix = '_' + suffix
        # If suffix is empty you don't really have one
            else:
                have_suffix = False
        except:
            have_suffix = False
            pass

        # If there is a non-empty suffix, copy to generic filename
        if have_suffix:
            try:
                shutil.copyfile('machine.inp' + suffix, 'machine.inp')
            except IOError as xxx_todo_changeme:
                (errno, strerror) = xxx_todo_changeme.args
                print(('Error copying file %s to %s' % ('machine.inp' + suffix, 'machine.inp', strerror)))
                logMsg = 'Error copying machine.inp_<suffix> -> machine.inp'
                services.exception(logMsg)
                raise

      # run TORIC init
        do_input = os.path.join(self.BIN_PATH, 'do_toric_init_abr')
        retcode = subprocess.call([do_input,cur_state_file])
        if (retcode != 0):
            logMsg = 'Error in call to toric_init'
            self.services.error(logMsg)
            raise Exception(logMsg)

      # Update plasma state files in plasma_state work directory
        try:
            services.update_state()
        except Exception:
            logMsg = 'Error in call to update_state()'
            self.services.exception(logMsg)
            raise 

      # Archive output files
      # N.B.  do_toric_init does not produce a complete set of TORIC output
      #       files.  This causes an error in stage_output_files().  To
      #       solve this we generate a dummy set of output files here with
      #       system call 'touch'
        for file in self.OUTPUT_FILES.split():
            print(('touching ', file))
            subprocess.call(['touch', file])
      # Now stage them
        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception:
            logMsg = 'Error in call to stage_output_files()'
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
        print('toric.restart() called')

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
        cur_state_file = get_global_param(self, services, 'CURRENT_STATE')
        cur_eqdsk_file = get_global_param(self, services, 'CURRENT_EQDSK')

        try:
             self.toric_log = os.path.join(workdir, 'log.toric')
        except:
            logMsg = 'rf_ic_toric_mcmd: error in getting log.toric path'
            self.services.exception(logMsg)
            raise 

        return 0

# ------------------------------------------------------------------------------
#
# STEP function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        """Take a step for the toric component.  Really a complete run."""
        print('toric.step() called')

        if (self.services == None) :
            logMsg = 'Error in toric: step (): No self.services'
            self.services.error(logMsg)
            raise Exception(logMsg)
        services = self.services

    # Get global configuration parameters
        cur_state_file = get_global_param(self, services, 'CURRENT_STATE')
        cur_eqdsk_file = get_global_param(self, services, 'CURRENT_EQDSK')

      # Copy plasma state files over to working directory
        try:
            services.stage_state()
        except:
            logMsg = 'Error in call to stage_state()'
            self.services.exception(logMsg)
            raise 

      # Get input files
        try:
            services.stage_input_files(self.INPUT_FILES)
        except:
            logMsg = 'Error in call to stage_input_files()'
            self.services.exception(logMsg)
            raise 

      # Copy machine.inp_<suffix> to generic file name -> machine.inp if there is
      # a suffix
        try:
            suffix = self.INPUT_SUFFIX
            have_suffix = True
        # If suffix is not empty put an underscore in front of it.
            if len(suffix) > 0:
                print(('INPUT_SUFFIX = ', suffix))
                suffix = '_' + suffix
        # If suffix is empty you don't really have one
            else:
                have_suffix = False
        except:
            have_suffix = False
            pass

        # If there is a non-empty suffix, copy to generic filename
        if have_suffix:
            try:
                shutil.copyfile('machine.inp' + suffix, 'machine.inp')
            except IOError as xxx_todo_changeme1:
                (errno, strerror) = xxx_todo_changeme1.args
                print(('Error copying file %s to %s' % ('machine.inp' + suffix,
                'machine.inp', strerror)))
                logMsg = 'Error copying machine.inp_<suffix> -> machine.inp'
                services.exception(logMsg)
                raise 

        prepare_input = os.path.join(self.BIN_PATH, 'prepare_toric_input_abr')
        process_output  = os.path.join(self.BIN_PATH, 'process_toric_output_mcmd')

        zero_RF_IC_power = self.ZERO_IC_POWER_BIN
        toric_bin = self.TORIC_BIN
        prepare_eqdsk  = self.GEQXPL_BIN

#        cur_state_file = self.plasma_state_file
#        cur_eqdsk_file = self.eqdsk_file
        toric_log = self.toric_log
        cwd = os.getcwd()

# Check if ICRF power is zero (or effectively zero).  If true don't run toric just
# run zero_RF_IC_power fortran code
        print('cur_state_file = ', cur_state_file)
        ps = Dataset(cur_state_file, 'r', format='NETCDF3_CLASSIC')
        rf_power = ps.variables['power_ic'][:]
        ps.close()
        total_rf_power = sum(rf_power)
        print(('Total IC power = ', total_rf_power))
        zero_rf_power_bin_string = 'ZERO_IC_POWER_BIN'
        power_ic = rf_power
        print('power = ', power_ic)
        if(-0.02 < power_ic < 0.02):
            print(zero_RF_IC_power)
            retcode = subprocess.call([zero_RF_IC_power, cur_state_file])
            if (retcode != 0):
                logMsg = 'Error executing ' + prepare_input
                self.services.error(logMsg)
                raise Exception(logMsg)

            # N.B. zero_RF_IC_power does not produce a complete set of TORIC output
            #      files.  This causes an error in stage_output_files().  To
            #      solve this we generate a dummy set of output files here with
            #      system call 'touch'
            for file in self.OUTPUT_FILES.split():
                subprocess.call(['touch', file])

# Check if ICRF power is negative.  If true don't run toric just
# retain power from previous time step i.e. leave sources untouched in the state.
# However power_ic needs to be reset back to positive

        elif( power_ic < -0.02):
            print('continuing power from previous time step')
            ps.variables['power_ic'].assignValue(-power_ic)
            ps.close()

    # Or actually run TORIC

        else:

            if not os.path.isfile(prepare_input):
                logMsg = 'Cannot fine TORIC prepare_input binary: ' + prepare_input
                self.services.error(logMsg)
                raise Exception(logMsg)

            # Call TORIC prepare_input to generate torica.inpp
            retcode = subprocess.call([prepare_input, cur_state_file]) #, cur_eqdsk_file])
            if (retcode != 0):
                logMsg = 'Error executing ' + prepare_input
                self.services.error(logMsg)
                raise Exception(logMsg)

            # Call xeqdsk_setup to generate eqdsk.out file
            print(('prepare_eqdsk', prepare_eqdsk, cur_eqdsk_file))

            retcode = subprocess.call([prepare_eqdsk, \
                                       '@equigs_gen', '/g_filename='+cur_eqdsk_file,\
                                       '/equigs_filename=equigs.data'])
            if (retcode != 0):
                logMsg = 'Error in call to prepare_eqdsk'
                self.services.error(logMsg)
                raise Exception(logMsg)

            # Launch TORIC executable
            print(('toric processors = ', self.NPROC))
            cwd = services.get_working_dir()
            task_id = services.launch_task(self.NPROC, cwd, toric_bin, logfile=toric_log)
            retcode = services.wait_task(task_id)
            if (retcode != 0):
                logMsg = 'Error executing command: ' + toric_bin
                self.services.error(logMsg)
                raise Exception(logMsg)

            # Call process_output
            # First rename default fort.* to expected names by component method as of toric5 r918 from ipp
            os.rename('fort.9','toric_cfg.nc')
            os.rename('fort.21','toric.nc')
            retcode = subprocess.call([process_output, cur_state_file])
            if (retcode != 0):
                logMsg = 'Error executing' + process_output
                self.services.error(logMsg)
                raise Exception(logMsg)


# Merge partial plasma state containing updated IC data
        try:
            partial_file = cwd + '/RF_IC_' + cur_state_file
            services.merge_current_state(partial_file, logfile='log.update_state')
            print(('merged TORIC plasma state data ', partial_file))
        except:
            logMsg = 'Error in call to merge_current_state(' + partial_file + ')'
            self.services.exception(logMsg)
            raise 

      # Archive output files
        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except:
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
        print('rf_ic_toric.checkpoint() called')
        services = self.services
        services.save_restart_files(timestamp, self.RESTART_FILES)


# ------------------------------------------------------------------------------
#
# FINALIZE function
# As of now it does nothing
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print('toric.finalize() called')
