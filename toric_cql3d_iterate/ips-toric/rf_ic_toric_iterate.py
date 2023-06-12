#! /usr/bin/env python

"""
TORIC component. SF cleanup 05/2023
"""

from builtins import range
import sys
import os
import subprocess
import getopt
import shutil
import string
import numpy as np
from netCDF4 import *
from  ipsframework import Component
from simple_file_editing_functions import put_lines
import toric_tools as tt

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
        try:
            NPROC = self.NPROC
            NPROC_QLDCI = self.NPROC_QLDCI
            NPROC_ABJ = self.NPROC_ABJ
            NPPN_QLDCI = self.NPPN_QLDCI
            NPPN_ABJ = self.NPPN_ABJ
            BIN_PATH = self.BIN_PATH
            TORIC_BIN=self.TORIC_BIN
            GEQXPL_BIN = self.GEQXPL_BIN
            ABJ_BIN = self.ABJ_BIN
            INPUT_DIR = self.INPUT_DIR
            INPUT_FILES = self.INPUT_FILES
            OUTPUT_FILES = self.OUTPUT_FILES
            RESTART_FILES = self.RESTART_FILES
        except:
            print('rf_ic_toric_iterate init: error getting cql3d-specific config parameters')
            services.error('rf_ic_toric_iterate: error getting toric-specific\
            config parameters')
            raise Exception('rf_ic_toric_iterate: error getting toric-specific\
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

        # SF reworked initializers do not require do_toric_init to run anymore
            
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

        # Ensure that everything has actually been initialized properly
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

        # Set preparation and process i/o script paths
        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_toric_input')
        process_output_bin  = os.path.join(self.BIN_PATH, 'process_toric_output')

        # Set physics exec bin locations
        toric_bin = self.TORIC_BIN
        prepare_eqdsk  = self.GEQXPL_BIN
        
        # Component specific config params
        self.TORIC_TIME_LIMIT = self.try_get_component_param(services,'TORIC_TIME_LIMIT', \
                                optional = True)
        if self.TORIC_TIME_LIMIT == None: 
            self.TORIC_TIME_LIMIT = "1e6" #some arbitrary big number
        self.NUM_TORIC_TRIES = self.try_get_component_param(services,'NUM_TORIC_TRIES', \
                                optional = True)
        if self.NUM_TORIC_TRIES == None: 
            self.NUM_TORIC_TRIES = "1"
        self.FORCE_DEFAULTS = self.try_get_component_param(services,'FORCE_DEFAULTS', \
                                optional = True)
        if self.FORCE_DEFAULTS == None: 
            self.FORCE_DEFAULTS = "0"

        time_limit = int(self.TORIC_TIME_LIMIT)
        num_tries  = int(self.NUM_TORIC_TRIES)
        arg_do_py_plots = self.DO_PY_PLOTS
        arg_force_defaults = self.FORCE_DEFAULTS
        
        # Get global configuration parameters
        cur_state_file = self.try_get_config_param(services,'CURRENT_STATE')
        cur_eqdsk_file = self.try_get_config_param(services,'CURRENT_EQDSK')
        cur_cql_file = self.try_get_config_param(services,'CURRENT_CQL')
        specs = self.try_get_config_param(services,'SPECS')
        arg_icrfpower = self.try_get_config_param(services,'RFPWR_IC')
        #optional global config params
        arg_enorm = '0.0'
        arg_nsurf_FP = '60'
        arg_rhoFPlo  = '0.01'
        arg_rhoFPhi  = '0.95'
        arg_enorm = self.try_get_config_param(services,'ENORM', optional = True, default = '0.0')
        arg_nsurf_FP = self.try_get_config_param(services,'NSURF_FP', optional = True, default = '60')
        arg_rhoFPlo = self.try_get_config_param(services,'RHO_FPLO', optional = True, default = '0.01')
        arg_rhoFPhi = self.try_get_config_param(services,'RHO_FPHI', optional = True, default = '0.95')

        # Get kwargs
        arg_toric_Mode = kwargs.get('toric_Mode', 'toric')
        arg_isol_Mode = kwargs.get('isol_Mode', '1')           
        arg_inumin_Mode = kwargs.get('inumin_Mode', 'Maxwell')
        arg_nphi_indx = kwargs.get('nphi_indx','0')

        # Set working directory and output path
        cwd = os.getcwd()
        toric_log = os.path.join(workdir, 'log.toric')
        if arg_toric_Mode == 'qldci':
            toric_log = os.path.join(workdir, 'log.toric_qldci')
        if arg_toric_Mode == 'qldci1':
            toric_log = os.path.join(workdir, 'log.toric_qldci1')
        if arg_toric_Mode == 'qldci2':
            toric_log = os.path.join(workdir, 'log.toric_qldci2')


        # Generate namelist for input preparation script
        nml_lines = ['&toric_prepare_nml \n']
        nml_lines.append(' cur_state_file = \"' + cur_state_file + '\",\n')
        nml_lines.append(' cur_geq_file = \"' + cur_eqdsk_file + '\",\n')    
        nml_lines.append(' arg_toric_mode = \"' + arg_toric_Mode + '\",\n')
        nml_lines.append(' arg_inumin_mode = \"' + arg_inumin_Mode + '\",\n')
        nml_lines.append(' arg_enorm = ' + arg_enorm + ',\n')
        nml_lines.append(' arg_specs = \"' + specs + '\",\n')
        nml_lines.append(' arg_nphi_indx = ' + arg_nphi_indx + ',\n')
        nml_lines.append(' arg_nfpsurf = ' + arg_nsurf_FP + ',\n')
        nml_lines.append(' arg_fplo = ' + arg_rhoFPlo + ',\n')
        nml_lines.append(' arg_fphi = ' + arg_rhoFPhi + ',\n')
        nml_lines.append(' force_defaults = ' + arg_force_defaults + ',\n')
        nml_lines.append(' /\n')
        put_lines('toric_prepare.nml',nml_lines)

        # Call toric prepare_input to generate toric.inp
        if not os.path.isfile(prepare_input_bin):
            logMsg = 'Cannot find toric prepare_input binary: ' + prepare_input_bin
            self.services.error(logMsg)
            raise Exception(logMsg)

        toric_prep_log = open('log_prepare_toric_input','w')
            
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                                   event_comment =  prepare_input_bin)
        retcode = subprocess.call(prepare_input_bin, stdout = toric_prep_log,\
                                  stderr = subprocess.STDOUT)
        if (retcode != 0):
            logMsg = 'Error executing ' + prepare_input_bin
            self.services.error(logMsg)
            raise Exception(logMsg)

        # Call xeqdsk_setup to generate eqdsk.out file
        cmd_eqdsk = [prepare_eqdsk, '@equigs_gen', '/g_filename='+cur_eqdsk_file,\
                                       '/equigs_filename=equigs.data']
        print('running', cmd_eqdsk)
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                                   event_comment =  cmd_eqdsk)
        eq_prep_log = open('log_geqxpl','w')
        retcode = subprocess.call(cmd_eqdsk, stdout = eq_prep_log,\
                                       stderr = subprocess.STDOUT)
        if (retcode != 0):
            logMsg = 'Error in call to prepare_eqdsk'
            self.services.error(logMsg)
            raise Exception(logMsg)

        cwd = services.get_working_dir()
        
        # For toric mode and nonMaxwellian run ABJ
        if arg_toric_Mode == 'toric' and arg_inumin_Mode == 'nonMaxwell':
            print('\nRunning ABJ')

            #Run preparation script for ABJ input
            log_file = open('log_prepare_abj_input', 'w')
            prepare_abj_input_bin = os.path.join(self.BIN_PATH, 'prepare_abj_input')

            command = prepare_abj_input_bin + ' ' + cur_state_file + ' ' + arg_enorm + ' ' + specs
            services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                                           event_comment = command)
            subprocess.call(command.split(), stdout = log_file,\
                            stderr = subprocess.STDOUT)
            if (retcode != 0):
                print('Error executing abj prep ', prepare_input_bin)
                services.error('Error executing abj prep')
                raise Exception('Error executing abj prep')

            #Copy prepared files to input file names 
            shutil.copyfile('ABJ.inp_new','ABJ.inp')
            shutil.copyfile('ABJ_driver.inp_new','ABJ_driver.inp')
                
            # Run ABJ
            abj_bin = self.ABJ_BIN
            cmd_abj=self.ABJ_BIN
            nproc_abj = int(self.NPROC_ABJ)
            nppn_abj = int(self.NPPN_ABJ)
            cpp_abj = int(2 * int(128/int(nppn_abj)))
            abj_log = os.path.join(workdir, 'log.abj')
            task_id = services.launch_task(nproc_abj,cwd,abj_bin, \
                                           logfile=abj_log,whole_nodes=True, \
                                           task_ppn = nppn_abj, task_cpp= cpp_abj)
            retcode = services.wait_task(task_id, timeout = 1800.0, delay = 10.0)
            if (retcode != 0):
                services.error("ABJ Failed")
                raise Exception("ABJ Failed")
            print('Finished ABJ')

        # Launch toric executable

        # Set number of processors depending on toric_mode
        run_nproc = self.NPROC
        run_nppn = self.NPPN
        if arg_toric_Mode == 'qldci' or arg_toric_Mode == 'qldci1' or arg_toric_Mode == 'qldci2':
            run_nproc = self.NPROC_QLDCI
            run_nppn = self.NPPN_QLDCI

        print('arg_toric_Mode = ', arg_toric_Mode, '   toric processors = ', run_nproc)
        coresPerProc = int(2 * int(128/int(run_nppn)))
            
        # Try to launch TORIC multiple times if TORIC_TRIES > 1 in config file
        for i in range(int(self.NUM_TORIC_TRIES)):
            print(' TORIC try number ', i + 1)
                
            task_id = services.launch_task(run_nproc, cwd, toric_bin, logfile=toric_log,\
                                           whole_nodes=True, task_ppn = run_nppn, task_cpp = coresPerProc)
            retcode = services.wait_task(task_id, timeout = time_limit,delay=15.0)
            if (retcode == 0):
                break
            else:
                services.error("TORIC failed after %d trials" % int(self.NUM_TORIC_TRIES))
                raise Exception("TORIC failed after %d trials" % int(self.NUM_TORIC_TRIES))

        if (arg_toric_Mode == 'toric') and (arg_do_py_plots=='1'):
            print('Plotting toric outputs')

            #find the min spec
            ps = Dataset(cur_state_file, 'r+', format = 'NETCDF3_CLASSIC')
            rfmin_indx = np.copy(ps.variables['rfmin_to_alla'][0])
            ps.close()
            rfmin_indx = rfmin_indx-1

            #run the toric plotting routine
            ICRun = tt.toric_analysis("fort.21")
            ICRun.mode='ICRF'
            ICRun.threeplots(min_indx = rfmin_indx)

        #SF modify this later, need to modify this filename when dealing with multiple toroidal modes/freq sources
        cur_toric_file = 'fort.9'
        nml_lines = ['&toric_process_nml \n']
        nml_lines.append(' cur_state_file = \"' + cur_state_file + '\",\n')
        nml_lines.append(' toric_output_file = \"' + cur_toric_file + '\",\n')    
        nml_lines.append(' /\n')
        put_lines('toric_process.nml',nml_lines)

        # Call toric prepare_input to generate toric.inp
        if not os.path.isfile(process_output_bin):
            logMsg = 'Cannot find toric process_output binary: ' + process_output_bin
            self.services.error(logMsg)
            raise Exception(logMsg)

        toric_proc_log = open('log_process_toric_output','w')
            
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                                   event_comment =  process_output_bin)

        retcode = subprocess.call(process_output_bin, stdout = toric_proc_log,\
                                  stderr = subprocess.STDOUT)

        if (retcode != 0):
            logMsg = 'Error executing ' + process_output_bin
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
        arg_save_output = kwargs.get('save_output', 'True')
        if arg_save_output == 'True':
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
