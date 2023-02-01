#! /usr/bin/env python

"""
TORIC component.  Adapted from rf_lh_torlh.py. (7-24-2015)

"""
# Working notes: DBB 9-29-2021
# Converting from rf_lh_torlh.py to rf_lh_toric.py.  Essentially this is just global search
# and replace.  This reverses change mentioned in Working notes: DBB 5-14-2016 below.

# Working notes:  DBB 10-5-2017 
# Because of random crashes on EDISON, we had previously introduced config parameter 
# TORLH_TIME_LIMIT so that if TORLH does crash it won't just sit there and burn up the 
# whole allocation.  Now making that an optional config parameter.
# Also now adding capability to do multiple tries of TORLH if it crashes or times out.
# To use this set optional config parameter NUM_TORLH_TRIES > 1.

# Working notes:  DBB 9-22-2017
# Added config parameter NPROC_QLDCE so that torlh can run in qldce mode with different
# number of processors than in toric mode.  Added config parameter TORLH_TIME_LIMIT to set
# a time limit for an individual run of TORLH using Wael's new optional arguments to 
# services.wait_task()
#
# Working notes: DBB 7-28-2017
# TORLH runs in either of two modes.  
# toric_Mode = 'toric' -> normal torlh run that solves wave equation
# toric_Mode = 'tqldce' -> does not solve wave equation, generates QL operator file toric_qlde.cdf
# This is controlled by optional keyword argument to STEP function 'toric_mode'
#
# There are 2 additional keyword args to STEP that set the parameters 'inumin' and 'isol'.
# These 3 parameters are communicated as command line args to the prepare_torlh_input_abr.f90
# code, which in turn writes them into the 'torica.inp' file, which is the actual input file
# to TORLH.
#  
# The optional command line arguments to prepare_torlh_input_abr.f90 are:
# arg_toric_Mode = 'toric' or 'qldce'
# arg_inumin_Mode = Maxwell (sets inumin = (0,0,0,0)) or nonMaxwell (sets inumin = (3,0,0,0))
# arg_isol_Mode = 0 or 1 (Normally = 1)
#
# The defaults are arg_toric_Mode = 'toric', arg_inumin_Mode = 'Maxwell', arg_isol_Mode = 1
# So if no keyword args are provided to the call to STEP the component functions as a normal
# TORLH RF component.
#
# Working notes: DBB 5-11-2017
# Modified to pick up enorm as global config parameter and pass to prepare_torlh_input
# through a command line argument
#
# Nota Bene: This component uses services.update_plasma_state() which overwrites all state
# files. To use this in a concurrent simulation should use merge_plasma_state instead.

# Working notes: DBB 8-29-2016
# Adapting to replicate more of the functionality of the TORLH/CQL3D iteration script
# Adding code to do IDL plots and to run cql3d_mapin
# Adding optional config parameters DO_IDL_PLOTS and RUN_CQL3D_MAPIN to config file

# Working notes: DBB 8-28-2016
# Have not yet developed a process_torlh_output code.  For the torlh/CQL3D coupling
# no lower hybrid data needs to go back into plasma state.  So for now have just commented
# out the calls to process_torlh_output and the merge_current_plasma_state.

# Working notes: DBB 5-14-2016
# Changed all IC references to LH.  Assuming torlh works just the same as TORIC.  
# Can modify later.

# Old notes from RF_LH_toric_abr_mcmd.py (Maybe delete them when torlh is working)
# Working notes:  DBB 2/9/2011
# Eliminated 'torlha.inp' as an input file.  The torlh code requires a 'torlha.inp'. However
# in the IPS this file is written either by the 'do_torlh_init_abr.f90' during 'init' or by
# the 'prepare_torlh_input_abr.f90' code during 'step'  Therefore I'm eliminating 'torlha.inp'as
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
#      About zero LH power:  We don't want to run and LH code when the RF power is zero.
#      In the past we set the LH source profiles to zero in this script and did
#      ps_store_plasma_state.  However with MCMD the STEP function needs to produce a partial
#      only plasma state containing the LH data (i.e. ps_write_update_file not
#      ps_store_plasma_state).  This really has to done with plasma state fortran code not not
#      with the python netcdf interface.So I wrote a simple code called zero_RF_LH_power to set
#      all LH source profiles in plasma state to zero and then write a partial plasma state.
#      For now this code lives in: /ips/trunk//components/rf/model_RF_LH and it gets built and
#      installed by the Makefile there.

from builtins import range
import sys
import os
import subprocess
import getopt
import shutil
import string
from netCDF4 import *
from  component import Component

class toric (Component):

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
        BIN_PATH = self.try_get_component_param(services,'BIN_PATH')
        RESTART_FILES = self.try_get_component_param(services,'RESTART_FILES')
        NPROC = self.try_get_component_param(services,'NPROC')
        NPROC_QLDCI = self.try_get_component_param(services,'NPROC_QLDCI', \
                                optional = True)

        self.TORIC_TIME_LIMIT = self.try_get_component_param(services,'TORIC_TIME_LIMIT', \
                                optional = True)
        if self.TORIC_TIME_LIMIT == None: 
            self.TORIC_TIME_LIMIT = -1
            
        self.NUM_TORIC_TRIES = self.try_get_component_param(services,'NUM_TORIC_TRIES', \
                                optional = True)
        if self.NUM_TORIC_TRIES == None: 
            self.NUM_TORIC_TRIES = 1
        
        toric_log = os.path.join(workdir, 'log.toric')

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

      # Copy machine.inp_<suffix> to generic file name -> machine.inp if there is
      # a suffix
        try:
            suffix = self.INPUT_SUFFIX
            have_suffix = True
        # If suffix is not empty put an underscore in front of it.
            if len(suffix) > 0:
                print('INPUT_SUFFIX = ', suffix)
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
                print('Error copying file %s to %s' % ('machine.inp' + suffix, 'machine.inp', strerror))
                logMsg = 'Error copying machine.inp_<suffix> -> machine.inp'
                services.exception(logMsg)
                raise

        do_init = os.path.join(self.BIN_PATH, 'do_toric_init_abr')
        retcode = subprocess.call([do_init, cur_state_file])
        if (retcode != 0):
            logMsg = 'Error in call to toric_init'
            self.services.error(logMsg)
            raise Exception(logMsg)
            
      # Update plasma state files in plasma_state work directory
        try:
            services.update_state()
        except Exception:
            logMsg = 'Error in call to update_plasma_state()'
            self.services.exception(logMsg)
            raise 

      # Archive output files
      # N.B.  do_toric_init does not produce a complete set of toric output
      #       files.  This causes an error in stage_output_files().  To
      #       solve this we generate a dummy set of output files here with
      #       system call 'touch'
        for file in self.OUTPUT_FILES.split():
            print('touching ', file)
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
        print('\ntoric.restart() called')

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
            self.toric_log = os.path.join(workdir, 'log.toric')
        except:
            logMsg = 'toric restart: error in getting config parameters'
            self.services.exception(logMsg)
            raise 
            
        return 0

# ------------------------------------------------------------------------------
#
# Run TORIC with toricmode, inumin, and isol as defined in argument list
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp, **kwargs):
        """Take a step for the toric component.  Really a complete run."""
        print('\ntoric.step() called')

        if (self.services == None):
            logMsg = 'Error in toric: step (): No self.services'
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

      # Copy machine.inp_<suffix> to generic file name -> machine.inp if there is
      # a suffix
        try:
            suffix = self.INPUT_SUFFIX
            have_suffix = True
        # If suffix is not empty put an underscore in front of it.
            if len(suffix) > 0:
                print('INPUT_SUFFIX = ', suffix)
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
                print('Error copying file %s to %s' % ('machine.inp' + suffix,
                'machine.inp', strerror))
                logMsg = 'Error copying machine.inp_<suffix> -> machine.inp'
                services.exception(logMsg)
                raise 

        prepare_input = os.path.join(self.BIN_PATH, 'prepare_toric_input_abr')
        process_output  = os.path.join(self.BIN_PATH, 'process_toric_output_mcmd')

        zero_RF_IC_power = self.ZERO_IC_POWER_BIN
        toric_bin = self.TORIC_BIN
        prepare_eqdsk  = self.GEQXPL_BIN

    # Get global configuration parameters
        cur_state_file = self.try_get_config_param(services,'CURRENT_STATE')
        cur_eqdsk_file = self.try_get_config_param(services,'CURRENT_EQDSK')
        cur_cql_file = self.try_get_config_param(services,'CURRENT_CQL')
        # enorm which is used here and in cql3d
        arg_enorm = 'None'
        arg_enorm = self.try_get_config_param(services,'ENORM', optional = True)

        toric_log = os.path.join(workdir, 'log.toric')
        cwd = os.getcwd()

# Check if IC power is zero (or effectively zero).  If true don't run toric just
# run zero_RF_IC_power fortran code
        print('cur_state_file = ', cur_state_file)
#         ps = NetCDFFile(cur_state_file, 'r')
#         power_ic = ps.variables['power_ic'].getValue()[0]
#         ps.close()
        ps = Dataset(cur_state_file, 'r+', format = 'NETCDF3_CLASSIC')
        power_ic = ps.variables['power_ic'][0]
        ps.close()
        
        print('power = ', power_ic)
        if(-0.0001 < power_ic < 0.0001):
            print(zero_RF_IC_power)
            services.send_portal_event(event_type = 'COMPONENT_EVENT',\
              event_comment =  'running ' + zero_RF_IC_power)
            retcode = subprocess.call([zero_RF_IC_power, cur_state_file])
            if (retcode != 0):
                logMsg = 'Error executing ' + prepare_input
                self.services.error(logMsg)
                raise Exception(logMsg)

            # N.B. zero_RF_IC_power does not produce a complete set of toric output
            #      files.  This causes an error in stage_output_files().  To
            #      solve this we generate a dummy set of output files here with
            #      system call 'touch'
            for file in self.OUTPUT_FILES.split():
                subprocess.call(['touch', file])

# Check if IC power is negative.  If true don't run toric just
# retain power from previous time step i.e. leave sources untouched in the state.
# However power_ic needs to be reset back to positive

        elif( power_ic < -0.02):
            print('continuing power from previous time step')
            ps.variables['power_ic'].assignValue(-power_ic)
            ps.close()
# ------------------------------------------------------------------------------                

    # Or actually run toric

        else:
            if not os.path.isfile(prepare_input):
                logMsg = 'Cannot find toric prepare_input binary: ' + prepare_input
                self.services.error(logMsg)
                raise Exception(logMsg)

            # Call toric prepare_input to generate toric.inp

            arg_toric_Mode = kwargs.get('toric_Mode', 'toric')
            arg_isol_Mode = kwargs.get('isol_Mode', '1')           
            arg_inumin_Mode = kwargs.get('inumin_Mode', 'Maxwell')
            if arg_toric_Mode == 'qldci':
                toric_log = os.path.join(workdir, 'log.toric_qldci')
            
            cmd_prepare_input = [prepare_input, cur_state_file, arg_toric_Mode,\
                      arg_inumin_Mode,arg_isol_Mode, arg_enorm]
                      
            print('running = ', cmd_prepare_input)
            services.send_portal_event(event_type = 'COMPONENT_EVENT',\
              event_comment =  cmd_prepare_input)
            retcode = subprocess.call(cmd_prepare_input)
            if (retcode != 0):
                logMsg = 'Error executing ' + prepare_input
                self.services.error(logMsg)
                raise Exception(logMsg)

            # Call xeqdsk_setup to generate eqdsk.out file
            cmd_eqdsk = [prepare_eqdsk, '@equigs_gen', '/g_filename='+cur_eqdsk_file,\
                                       '/equigs_filename=equigs.data']
            print('running', cmd_eqdsk)
            services.send_portal_event(event_type = 'COMPONENT_EVENT',\
              event_comment =  cmd_eqdsk)
            retcode = subprocess.call(cmd_eqdsk)
            if (retcode != 0):
                logMsg = 'Error in call to prepare_eqdsk'
                self.services.error(logMsg)
                raise Exception(logMsg)

            cwd = services.get_working_dir()

            # For toric mode and nonMaxwellian run ABJ
            # SF I originally wanted the dielectric representation to be its own 
            if arg_toric_Mode == 'toric' and arg_inumin_Mode == 'nonMaxwell':
                print('\nRunning ABJ')

                #Copy the cql3d output file to the expected file name for abj script
                try:
                    subprocess.call(['cp', cur_cql_file, 'cql3d.cdf' ])
                except Exception:
                    message = 'rf_ic_toric_iterate: Error copying CURRENT_CQL_FILE to cql3d.cdf'
                    print(message)
                    services.exception(message)
                    raise              

                #Run preparation script for ABJ input
                #log_file = open('log_prepare_abj_input', 'w')
                #prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_abj_input')

                #command = prepare_input_bin + ' ' + cur_state_file + ' ' + arg_enorm
                #services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                #                           event_comment = command)
                #subprocess.call(command.split(), stdout = log_file,\
                #                  stderr = subprocess.STDOUT)
                #if (retcode != 0):
                #    print('Error executing abj prep ', prepare_input_bin)
                #    services.error('Error executing abj prep')
                #    raise Exception('Error executing abj prep')

                #shutil.copyfile('ABJ.inp_new','ABJ.inp')
                #shutil.copyfile('ABJ_driver.inp_new','ABJ_driver.inp')
        
                # Run ABJ
                abj_bin = self.ABJ_BIN
                cmd_abj=self.ABJ_BIN
                nproc_abj = self.NPROC_ABJ
                abj_log = os.path.join(workdir, 'log.abj')
                task_id = services.launch_task(nproc_abj,cwd,abj_bin, logfile=abj_log)
                retcode = services.wait_task(task_id, timeout = 1800.0, delay = 60.)
                if (retcode != 0):
                    services.error("ABJ Failed")
                    raise Exception("ABJ Failed")
                print('Finished ABJ')

            # Launch toric executable
            # Set number of processors depending on toric_mode
            run_nproc = self.NPROC
            if arg_toric_Mode == 'qldci':
                run_nproc = self.NPROC_QLDCI

            print('arg_toric_Mode = ', arg_toric_Mode, '   toric processors = ', run_nproc)
            time_limit = float(self.TORIC_TIME_LIMIT)
            # Try to launch TORIC multiple times if TORIC_TRIES > 1 in config file
            for i in range(int(self.NUM_TORIC_TRIES)):
                print(' TORIC try number ', i + 1)
                task_id = services.launch_task(run_nproc, cwd, toric_bin, logfile=toric_log)
                retcode = services.wait_task(task_id, timeout = time_limit, delay = 60.)
                if (retcode == 0):
                    break
            else:
                services.error("TORIC failed after %d trials" % int(self.NUM_TORIC_TRIES))
                raise Exception("TORIC failed after %d trials" % int(self.NUM_TORIC_TRIES))


#             print 'arg_toric_Mode = ', arg_toric_Mode, '   toric processors = ', run_nproc
#             task_id = services.launch_task(run_nproc, cwd, toric_bin, logfile=toric_log)
#             time_limit = float(self.TORIC_TIME_LIMIT)
#             retcode = services.wait_task(task_id, timeout = time_limit, delay = 60.)
#             if (retcode != 0):
#                 logMsg = 'Error executing command: ' + toric_bin
#                 self.services.error(logMsg)
#                 raise Exception(logMsg)

            # SF removed for now this data is kept in log.toric right now
            # Preserve torica.out from run to distinguish toric mode = 'toric' from 'qldce'
            #new_file_name = 'torica_' + arg_toric_Mode + '.out'
            #try:
            #    shutil.copyfile('torica.out', new_file_name)
            #except IOError as xxx_todo_changeme2:
            #    (errno, strerror) = xxx_todo_changeme2.args
            #    logMsg =  'Error copying file %s to %s' % ('torica.out', new_file_name\
            #            , strerror)
            #    print(logMsg)
            #    services.exception(logMsg)
            #    raise

            #SF removed mapin direct mapping now used
            # For qldce mode need to also run mapin
            #if arg_toric_Mode == 'qldci':
            #    mapin_bin = self.try_get_component_param(services,'MAPIN_BIN')
            #    print('\nRunning ' + mapin_bin)
            #    services.send_portal_event(event_type = 'COMPONENT_EVENT', \
            #         event_comment = 'running ' + mapin_bin)
            #    retcode = subprocess.call([mapin_bin])
            #    if (retcode != 0):
            #        logMsg = 'Error executing ' + mapin_bin
            #        self.services.error(logMsg)
            #        raise Exception(logMsg)

# For toric mode merge partial plasma state containing updated IC data
        if arg_toric_Mode == 'toric':
            try:
                partial_file = cwd + '/RF_IC_' + cur_state_file
                # No process_output code yet
                #services.merge_current_plasma_state(partial_file, logfile='log.update_state')
                #print 'merged toric plasma state data ', partial_file
                print('No process_output code yet, so no plasma state merge')
            except:
                logMsg = 'Error in call to merge_current_plasma_state(' + partial_file + ')'
                self.services.exception(logMsg)
                raise 

      # Update plasma state files in plasma_state work directory
        try:
            services.update_state()
        except Exception:
            logMsg = 'Error in call to update_plasma_state()'
            self.services.exception(logMsg)
            raise 

      # Archive output files
        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except:
            logMsg = 'Error in call to stage_output_files()'
            self.services.exception(logMsg)
            raise 

        if 'last_pwr' in kwargs:
            ps = Dataset(cur_state_file, 'w', format = 'NETCDF3_CLASSIC')
            ps.variables['power_ic'][0] = saved_pwr
            ps.close()

        return 0

# ------------------------------------------------------------------------------
#
# checkpoint function
# Saves plasma state files to restart directory
#
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print('RF_IC_toric.checkpoint() called')
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


# ------------------------------------------------------------------------------
#
# "Private"  methods
#
# ------------------------------------------------------------------------------

    def try_get_config_param(self, services, param_name, optional=False):

        try:
            value = services.get_config_param(param_name)
            print(param_name, ' = ', value)
        except Exception:
            if optional:
                print('config parameter ', param_name, ' not found')
                value = None
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
