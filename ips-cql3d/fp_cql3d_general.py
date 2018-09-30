#! /usr/bin/env python

# fp_cql3d_general.py, Version 0.2 Batchelor 5-16-2017
# Added code to pick up enorm from global parameter in config file and add to argument list
# of prepare_cql3d_input.f90

# fp_cql3d_general.py, Version 0.1 Batchelor 3-2-2017
# Added cur_cql_file only to update_plasma_state for iteration (not full update).
# Plasma state variables are already updated with merge_current_plasma_state.  So this 
# component can be used concurrently.

# Version 0.0, fp_cql3d_general.py, BH 08/15/2012
# General (LH,EC,IC,RW,NBI, etc) cql3d component .py file.
# Output is controlled by the input namelist file, cqlinput.

# Derived from rf_genray_LH.py.
# Changes:
#  *replace genray.in ==> cqlinput
#  *rf_genray_LH ==> fp_cql3d_general
#  *replace RF_GENRAY ==> FP_CQL3D,  rf_genray ==> fp_cql3d
#  *cql3dnml (equiv. of genraynml) is not used.  Adjust code,
#   along with other changes for prepare_cql3d_input/process_cql3d_output
#   command line arguments..
#  *replace GENRAY ==> CQL3D, genray ==> cql3d
#  *check remaining uses of _lh, lh, replace with _general
#
#
# LH version 0.0 11/19/2010 notes:
# This version distinguishes between the different RF components that GENRAY can 
# implement.  In this case the LH component.  So far the only place this affects
# is in the writing the partial plasma state that the framework needs to merge.
#
# Note: To merge plasma states the IPS framework expects the component to
# produce a partial state file with the component-specific name:
# <component>_<curr_state_filename>.  GENRAY works for EC, LH, and IC, i.e.
# it can implement several different components.  Therefore process_genray_output_mcmd.f90
# writes a partial state file with generic name 'RF_GENRAY_PARTIAL_STATE' and delegetes
# to this python component script the copying of this to the proper 
# component-specific update file name.  In this case RF_LH_<curr_state_filename>.

# version 1.0 9/29/2010 (Batchelor):
# This version adds checkpoint and restart functions and makes the exception 
# handling more uniform.

# version 0.0 3/1/2010 (Batchelor)



# Adjusted coding (BH)
# ------------------------------------------------------------------------------
#
# FP_CQL3D component script to drive GENRAY ray tracing code.
#
# Workflow:
# 1) 'init' function:
#       Stage input files: cqlinput_<suffix>
#       Copy to generic input file name -> cqlinput
#       Copy current plasma state file to generic name -> cur_state.cdf
#       Launch prepare_cql3d_input with command line arg 'init'.  This does
#       plasma state intitialization. Redirect output to log file
#       Copy generic state file cur_state.cdf -> current plasma state file 
#       Update plasma state and stage output files

# 2) 'step' function:
#       Copy current plasma state file to generic name -> cur_state.cdf
#       Launch prepare_cql3d_input with command line arg 'step'.  This gets
#       data needed from plasma state and writes a new cqlinput file.
#       Copy cqlinput to cqlinput_<new name> for archive.
#       Delete log file.
#       Launch cql3d executable, redirect output to new log file
#       Launch process_cql3d_output    
#       Copy generic state file cur_state.cdf -> current plasma state file 
#       Update plasma state and stage output files
#
# N.B. There is a good explanation of the command line arguments and program
#      operation in the prepare_cql3d_input.f90 source.
# ------------------------------------------------------------------------------

import sys
import os
import subprocess
import getopt
import shutil
import string
from  component import Component
from Numeric import *                        #Use numpy instead?? BH
from Scientific.IO.NetCDF import *           #Use scipy.io.netcdf implentation??  BH
from get_IPS_config_parameters import get_global_param, get_component_param

class cql3d(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        self.firstTime  = True
        self.curTime = -1.0
        self.prevTime = -1.0
        print 'Created %s' % (self.__class__)

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp=0):
        print 'cql3d.init() called'

        services = self.services

    # Get global configuration parameters
        cur_state_file = get_global_param(self, services, 'CURRENT_STATE')
        cur_eqdsk_file = get_global_param(self, services, 'CURRENT_EQDSK')
        cur_dql_file = get_global_param(self, services, 'CURRENT_DQL')
        cur_cql_file = get_global_param(self, services,'CURRENT_CQL')
        cur_ImChizz_inp_file = get_global_param(self, services,'CURRENT_ImChizz_inp', optional = True)

    # Get component-specific configuration parameters. Note: Not all of these are
    # used in 'init' but if any are missing we get an exception now instead of
    # later

        NPROC = get_component_param(self, services, 'NPROC')
        BIN_PATH = get_component_param(self, services, 'BIN_PATH')
        INPUT_FILES = get_component_param(self, services, 'INPUT_FILES')
        OUTPUT_FILES = get_component_param(self, services, 'OUTPUT_FILES')
        RESTART_FILES = get_component_param(self, services, 'RESTART_FILES')
        BIN_PATH = get_component_param(self, services, 'BIN_PATH')
        CQL3D_BIN = get_component_param(self, services, 'CQL3D_BIN')
        cql3d_mode = get_component_param(self, services, 'CQL3D_MODE')
        cql3d_output = get_component_param(self, services, 'CQL3D_OUTPUT')
        cql3d_nml = get_component_param(self, services, 'CQL3D_NML')
        nsteps_str = get_component_param(self, services, 'NSTEPS_STR')
        deltat_str = get_component_param(self, services, 'DELTAT_STR')
        ps_add_nml = get_component_param(self, services, 'PS_ADD_NML')

        # enorm which is used here and in cql3d
        arg_enorm = get_component_param(self, services, 'ENORM', optional = True)

    # Copy plasma state files over to working directory
        try:
          services.stage_plasma_state()
        except Exception, e:
          print 'Error in call to stage_plasma_state()' , e
          services.error('Error in call to stage_plasma_state()')
          raise Exception, 'Error in call to stage_plasma_state()'
        
    # Get input files  
        try:
          services.stage_input_files(INPUT_FILES)
        except Exception, e:
          print 'Error in call to stageInputFiles()' , e
          services.error('Error in call to stageInputFiles()')
          raise Exception, 'Error in call to stageInputFiles()'
            
    # Copy cqlinput_<suffix> to generic file name -> cqlinput if there is
    # a suffix
        try:
          suffix = self.INPUT_SUFFIX
          have_suffix = True
        # If suffix is not empty put an underscore in front of it.
          if len(suffix) > 0:
              print 'cql3d INPUT_SUFFIX = ', suffix      
              suffix = '_' + suffix
        # If suffix is empty you don't really have one
          else:
              have_suffix = False
        except:
          have_suffix = False
      
      # If there is a non-empty suffix, copy to generic filename 'cqlinput'
        if have_suffix:      
          try:
              shutil.copyfile('cqlinput' + suffix, 'cqlinput')
          except IOError, (errno, strerror):
              print 'Error copying file %s to %s' % ('cqlinput' + suffix, 
              'cqlinput', strerror)
              services.error('Error copying cqlinput_<suffix> -> cqlinput')
              raise Exception, 'Error copying cqlinput_<suffix> -> cqlinput'

    # Copy current plasma state file to generic name -> cur_state.cdf
        try:
            shutil.copyfile(cur_state_file, 'cur_state.cdf')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % (cur_state_file, 'cur_state.cdf', strerror)
            services.error('Error copying cur_state_file -> cur_state.cdf')
            raise Exception, 'Error copying cur_state_file -> cur_state.cdf'

    # Copy current eqdsk file to generic name -> eqdsk
        try:
            shutil.copyfile(cur_eqdsk_file, 'eqdsk')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % (cur_eqdsk_file, 'eqdsk', strerror)
            services.error('Error copying cur_eqdsk_file -> eqdsk')
            raise Exception, 'Error copying cur_eqdsk_file -> eqdsk'

    # Copy current Dql file to generic name -> genray.nc
        if cur_dql_file != 'genray.nc':
			try:
				shutil.copyfile(cur_dql_file, 'genray.nc')
			except IOError, (errno, strerror):
				print 'Error copying file %s to %s' % (cur_dql_file, 'genray.nc', strerror)
				services.error('Error copying cur_dql_file -> genray.nc')
				raise Exception, 'Error copying cur_dql_file -> genray.nc'

        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_cql3d_input')

    # Call prepare_input - init
        print 'fp_cql3d: calling prepare_input init'        
        log_file = open('log_prepare_cql3d_input_init', 'w')
        ips_mode = 'init'

    # ptb: Set restart depending on whether the distribution function file from cql3d
    # already exist
        restart = 'disabled'
    #ptb&SS    if os.path.isfile('./cql3d.nc'): restart = 'enabled'
    #ptb&SS    if os.path.isfile('./cql3d.nc'): shutil.copyfile('cql3d.nc','distrfunc.nc')

        if os.path.isfile('./cql3d.nc'): 
            if os.stat('./cql3d.nc')[6] != 0:
                 restart = 'enabled'
                 shutil.copyfile('cql3d.nc','distrfunc.nc')

    # ptb: End of ptb hack

    # ptb:    command = prepare_input_bin + ' ' + ips_mode + ' ' + cql3d_mode  + ' ' +\
    # ptb:    cql3d_output + ' ' + cql3d_nml + ' ' + nsteps_str + ' ' + ps_add_nml
        command = prepare_input_bin + ' ' + ips_mode + ' ' + cql3d_mode  + ' ' +\
        cql3d_output + ' ' + cql3d_nml + ' ' + restart + ' ' + nsteps_str + ' ' +\
        ' ' + deltat_str + ' ' + ps_add_nml
        if arg_enorm != None:
        	command = command  + ' ' + arg_enorm
        
        print 'running = ', command
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
          event_comment =  command)
        
        retcode = subprocess.call(command.split(), stdout = log_file,\
                                  stderr = subprocess.STDOUT)
        if (retcode != 0):
            print 'Error executing cql3d init ', prepare_input_bin
            services.error('Error executing cql3d init')
            raise Exception, 'Error executing cql3d init'

    # Copy generic cur_state.cdf -> current plasma state file
        try:
            shutil.copyfile('cur_state.cdf', cur_state_file)
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % ('cur_state.cdf', cur_state_file, strerror)
            services.error('Error copying cur_state.cdf -> cur_state_file')
            raise Exception, 'Error copying cur_state.cdf -> cur_state_file'

    # Update plasma state files in plasma_state work directory
        try:
          services.update_plasma_state()
        except Exception, e:
          print 'Error in call to update_plasma_state()', e
          services.error('Error in call to update_plasma_state()')
          raise Exception, 'Error in call to update_plasma_state()'
     
    # Archive output files
    # N.B.  prepare_cql3d_input in init mode does not produce a complete set 
    #       of ourput files.  This causes an error in stage_output_files().
    #       To solve this we generate a dummy set of output files here with 
    #       system call 'touch'
        for file in self.OUTPUT_FILES.split():
          print 'touching ', file
          subprocess.call(['touch', file])
    # Now stage them      
        try:
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
          print 'Error in call to stage_output_files()', e
          services.error('Error in call to stage_output_files()')
          raise Exception, 'Error in call to stage_output_files()'

        return 0
# ------------------------------------------------------------------------------
#
# RESTART function
# Gets restart files from restart directory
# Loads the global configuration parameters from the config file
#
# ------------------------------------------------------------------------------
        
    def restart(self, timeStamp):
      print 'cql3d.restart() called'

      services = self.services
      if (services == None) :
         print 'Error in cql3d: step (): No services'
         services.error('Error in cql3d: step (): No services')
         raise Exception, 'Error in cql3d: step (): No services'

#ptb&SS      services = self.services

    # Get restart files listed in config file.        
      try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
      except Exception, e:
            print 'Error in call to get_restart_files()' , e
            services.error('cql3d: error in call to get_restart_files()')
            raise Exception, 'cql3d: error in call to get_restart_files()'

    # Get global configuration parameters
      try:
            self.plasma_state_file = services.get_config_param('CURRENT_STATE')
            self.eqdsk_file = services.get_config_param('CURRENT_EQDSK')
      except:
            print 'cql3d restart: error in getting config parameters'
            services.error('cql3d restart: error in getting config parameters')
            raise Exception, 'cql3d restart: error in getting config parameters'

      return 0
        
# ------------------------------------------------------------------------------
#
# STEP function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'cql3d.step() called'

        if (self.services == None) :
           print 'Error in cql3d: step (): No services'
           raise Exception, 'Error in cql3d: step (): No services'
        services = self.services

        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_cql3d_input')
        process_output_bin  = os.path.join(self.BIN_PATH, 'process_cql3d_output')
        cql3d_bin = os.path.join(self.BIN_PATH, 'cql3d')

    # Copy plasma state files over to working directory
        try:
          services.stage_plasma_state()
        except Exception, e:
          print 'Error in call to stage_plasma_state()' , e
          services.error('Error in call to stage_plasma_state()')
          raise Exception, 'Error in call to stage_plasma_state()'

    # Get global configuration parameters
        cur_state_file = get_global_param(self, services, 'CURRENT_STATE')
        cur_eqdsk_file = get_global_param(self, services, 'CURRENT_EQDSK')
        cur_dql_file = get_global_param(self, services, 'CURRENT_DQL')
        cur_cql_file = get_global_param(self, services,'CURRENT_CQL')
        cur_ImChizz_inp_file = get_global_param(self, services,'CURRENT_ImChizz_inp', optional = True)
        
        print 'CURRENT_CQL = ', cur_cql_file
    # Copy current plasma state file to generic name -> cur_state.cdf
        try:
            shutil.copyfile(cur_state_file, 'cur_state.cdf')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % (cur_state_file, 'cur_state.cdf', strerror)
            raise

    # Copy current eqdsk file to generic name -> eqdsk
        try:
            shutil.copyfile(cur_eqdsk_file, 'eqdsk')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % (cur_eqdsk_file, 'eqdsk', strerror)
            services.error('Error copying cur_eqdsk_file -> eqdsk')
            raise Exception, 'Error copying cur_eqdsk_file -> eqdsk'

# Check if LHRF power is zero (or effectively zero).  If true don't run Genray

        ps = NetCDFFile(cur_state_file, 'r')
        power_lh = ps.variables['power_lh'].getValue()[0]
        ps.close()
        print 'power = ', power_lh
        if(power_lh > 1.0E-04):

    # Copy current Dql file to generic name -> genray.nc
          if cur_dql_file != 'genray.nc':
			try:
				shutil.copyfile(cur_dql_file, 'genray.nc')
			except IOError, (errno, strerror):
				print 'Error copying file %s to %s' % (cur_dql_file, 'genray.nc', strerror)
				services.error('Error copying cur_dql_file -> genray.nc')
				raise Exception, 'Error copying cur_dql_file -> genray.nc'

          cql3d_mode = get_component_param(self, services, 'CQL3D_MODE')
          cql3d_output = get_component_param(self, services, 'CQL3D_OUTPUT')
          cql3d_nml = get_component_param(self, services, 'CQL3D_NML')
          nsteps_str = get_component_param(self, services, 'NSTEPS_STR')
          deltat_str = get_component_param(self, services, 'DELTAT_STR')
          ps_add_nml = get_component_param(self, services, 'PS_ADD_NML')
          cur_ImChizz_inp_file = get_global_param(self, services,'CURRENT_ImChizz_inp', optional = True)

        # enorm which is used here and in cql3d
          arg_enorm = get_component_param(self, services, 'ENORM', optional = True)
          
    # Call prepare_input - step
          print 'fp_cql3d step: calling prepare_input'
          
          log_file = open('log_prepare_cql3d_input_step', 'w')
          ips_mode = 'step'

    # ptb: Set restart depending on whether the distribution function file from cql3d
    # already exist
          restart = 'disabled'
    #ptb&SS      if os.path.isfile('./cql3d.nc'): restart = 'enabled'
    #ptb&SS      if os.path.isfile('./cql3d.nc'): shutil.copyfile('cql3d.nc','distrfunc.nc')

          if os.path.isfile('./cql3d.nc'): 
              if os.stat('./cql3d.nc')[6] != 0:
                   restart = 'enabled'
                   shutil.copyfile('cql3d.nc','distrfunc.nc')


    # ptb: End of ptb hack

    # ptb:    command = prepare_input_bin + ' ' + ips_mode + ' ' + cql3d_mode  + ' ' +\
    # ptb:    cql3d_output + ' ' + cql3d_nml + ' ' + nsteps_str + ' ' + ps_add_nml
          command = prepare_input_bin + ' ' + ips_mode + ' ' + cql3d_mode  + ' ' +\
          cql3d_output + ' ' + cql3d_nml + ' ' + restart + ' ' + nsteps_str + ' ' +\
          ' ' + deltat_str + ' ' + ps_add_nml
          if arg_enorm != None:
        	command = command  + ' ' + arg_enorm

          print 'running', command
          services.send_portal_event(event_type = 'COMPONENT_EVENT',\
              event_comment =  command)

          retcode = subprocess.call(command.split(), stdout = log_file,\
                                    stderr = subprocess.STDOUT)
          if (retcode != 0):
              print 'Error executing cql3d: ', prepare_input_bin
              services.error('Error executing cql3d prepare_input')
              raise Exception, 'Error executing cql3d prepare_input'
        
#     Launch cql3d - N.B: Path to executable is in config parameter CQL3D_BIN
          print 'fp_cql3d: launching cql3d'
# ptb: Need to first copy the cqlinput_new file to cqlinput
          shutil.copyfile('cqlinput_new', 'cqlinput')
          cwd = services.get_working_dir()
          task_id = services.launch_task(self.NPROC, cwd, self.CQL3D_BIN, logfile='log.cql3d')
          retcode = services.wait_task(task_id)
          if (retcode != 0):
              print 'Error executing command: ', cql3d_bin
              services.error('Error executing cql3d')
              raise Exception, 'Error executing cql3d'

    # Call process_output - step
          print 'fp_cql3d step: calling process_output'
          
          log_file = open('log_process_cql3d_output', 'w')
          mode = 'step'
# ptb;          command = process_output_bin + ' ' +  cql3d_mode
          command = process_output_bin + ' ' +  cql3d_output    

          print 'running', command
          services.send_portal_event(event_type = 'COMPONENT_EVENT',\
              event_comment =  command)

          retcode = subprocess.call(command.split(), stdout = log_file,\
                                    stderr = subprocess.STDOUT)                                  
          if (retcode != 0):
              print 'Error executing cql3d init ', process_output_bin
              services.error('Error executing cql3d process_output')
              raise Exception, 'Error executing cql3d process_output'
        
    # Copy generic cql3d partial plasma state file -> FP_CQL3D_cur_state_file  [correct??, BH]
          try:
              partial_file = cwd + '/FP_CQL3D_' + cur_state_file
              shutil.copyfile('FP_CQL3D_PARTIAL_STATE', partial_file )
          except IOError, (errno, strerror):
              print 'Error copying file %s to %s' % ('cur_state.cdf', cur_state_file, strerror)
              raise


# Merge partial plasma state containing updated CQL3D data  [BH]
          try:
             services.merge_current_plasma_state(partial_file, logfile='log.update_state')
             print 'merged CQL3D plasma state data ', partial_file
          except Exception, e:
             print 'Error in call to merge_current_plasma_state(' , partial_file, ')'
             self.services.error('Error in call to merge_current_plasma_state')
             raise Exception, 'Error in call to merge_current_plasma_state'

      # Update plasma state files in plasma_state work directory, but only cur_cql_file
      # This way it can be used concurrently without overwriting other plasma state files.
          print 'CURRENT_CQL = ', cur_cql_file
          try:
#            services.update_plasma_state(plasma_state_files = cur_cql_file)
            services.update_plasma_state([cur_cql_file, cur_ImChizz_inp_file])
          except Exception:
            logMsg = 'Error in call to update_plasma_state()'
            self.services.exception(logMsg)
            raise 
            
    # Archive output files
          try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
          except Exception, e:
            print 'Error in call to stage_output_files()', e
            services.error('Error in call to stage_output_files()')
            raise Exception, 'Error in call to stage_output_files()'
        
#     rename the log file so that it is not appended next step
#        os.rename('log.cql3d', this_logfile)

          return 0

        return 0  # return on "zero" LHRF power condition
# ------------------------------------------------------------------------------
#
# checkpoint function
# Saves plasma state files to restart directory
#
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
      print 'cql3d.checkpoint() called'
      if (self.services == None) :
         print 'Error in cql3d: checkpoint(): No services'
         services.error('Error in cql3d: checkpoint(): No services')
         raise Exception, 'Error in cql3d: checkpoint(): No services'
      services = self.services
      services.save_restart_files(timestamp, self.RESTART_FILES)

      return 0
        
# ------------------------------------------------------------------------------
#
# finalize function
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print 'cql3d.finalize() called'

# ------------------------------------------------------------------------------
#
# Internal utility functions
#
# ------------------------------------------------------------------------------

    # This doesn't seem to actually be used - ask Lee
    def fix_namelist(self, in_file, out_file):
        fd = open(in_file, 'r')
        line_list = fd.readlines()

        buf = ''.join(line_list)
        out_buf = ''

        in_string = False
        for c in buf:
            string = ''
            if in_string:
                if (c == '\''):
                    out_buf = out_buf + c
                    in_string = False
                elif (c == '\n'):
                    pass
                else :
                    out_buf = out_buf + c
            else:
                if (c == '\''):
                    in_string = True
                out_buf = out_buf + c
        ofd = open(out_file,'w')
        ofd.write(out_buf)
        ofd.close()
        fd.close()
        return