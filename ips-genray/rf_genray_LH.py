#! /usr/bin/env python

# LH version 0.0 11/19/2010
# This version distinguishes between the different RF components that GENRAY can 
# implement.  In this case the LH component.  So far the only place this affects
# is in the writing the partial plasma state that the framework needs to merge.
#
# Note: To merge plasma states the IPS framework expects the component to
# produce a partial state file with the component-specific name:
# <component>_<cur_state_filename>.  GENRAY works for EC, LH, and IC, i.e.
# it can implement several different components.  Therefore process_genray_output_mcmd.f90
# writes a partial state file with generic name 'RF_GENRAY_PARTIAL_STATE' and delegetes
# to this python component script the copying of this to the proper 
# component-specific update file name.  In this case RF_LH_<curr_state_filename>.

# version 1.0 9/29/2010 (Batchelor)
# This version adds checkpoint and restart functions and makes the exception 
# handling more uniform.

# version 0.0 3/1/2010 (Batchelor)

# ------------------------------------------------------------------------------
#
# RF_GENRAY component script to drive GENRAY ray tracing code.
#
# Workflow:
# 1) 'init' function:
#       Stage input files: genray.in_<suffix>
#       Copy to generic input file name -> genray.in
#       Copy current plasma state file to generic name -> cur_state.cdf
#       Launch prepare_genray_input with command line arg 'init'.  This does
#       plasma state intitialization. Redirect output to log file
#       Copy generic state file cur_state.cdf -> current plasma state file 
#       Update plasma state and stage output files

# 2) 'step' function:
#       Copy current plasma state file to generic name -> cur_state.cdf
#       Launch prepare_genray_input with command line arg 'step'.  This gets
#       data needed from plasma state and writes a new genray.in file.
#       Copy genray.in to genray.in_<new name> for archive.
#       Delete log file.
#       Launch genray executable, redirect output to new log file
#       Launch process_genray_output    
#       Copy generic state file cur_state.cdf -> current plasma state file 
#       Update plasma state and stage output files
#
# N.B. There is a good explanation of the command line arguments and program
#      operation in the prepare_genray_input.f90 source.
# ------------------------------------------------------------------------------

import sys
import os
import subprocess
import getopt
import shutil
import string
from ipsframework import Component
from Numeric import *
from Scientific.IO.NetCDF import *

class genray_LH(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp=0):
        print 'genray.init() called'

        services = self.services

    # Get global configuration parameters
        try:
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            #BH: this refers to rf data
            cur_dql_file = services.get_config_param('CURRENT_DQL')
            #BH: Could be used in future?
            #cql_file = services.get_config_param('CURRENT_CQL')
        except:
            print 'rf_genray_LH: error getting config parameters CURRENT_STATE CURRENT_EQDSK'
            services.error('rf_genray_LH: error getting config parameters CURRENT_STATE CURRENT_EQDSK')
            raise Exception, 'rf_genray_LH: error getting config parameters CURRENT_STATE CURRENT_EQDSK'

    # Get component-specific configuration parameters. Note: Not all of these are
    # used in 'init' but if any are missing we get an exception now instead of
    # later
        try:
            NPROC = self.NPROC
            BIN_PATH = self.BIN_PATH
            INPUT_FILES = self.INPUT_FILES
            OUTPUT_FILES = self.OUTPUT_FILES
            RESTART_FILES = self.RESTART_FILES
            BIN_PATH = self.BIN_PATH
            GENRAY_BIN = self.GENRAY_BIN
            RFMODE = self.RFMODE
            ISOURCE_STRING = self.ISOURCE_STRING
            GENRAYNML = self.GENRAYNML
            ADJ_READ = self.ADJ_READ
            PS_ADD_NML = self.PS_ADD_NML
        except:
            print 'rf_genray_LH init: error getting genray-specific config parameters'
            services.error('rf_genray_LH: error getting genray-specific\
            config parameters')
            raise Exception, 'rf_genray_LH: error getting genray-specific\
            config parameters'

    # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception, e:
          print 'Error in call to stage_state()' , e
          services.error('Error in call to stage_state()')
          raise Exception, 'Error in call to stage_state()'
        
    # Get input files  
        try:
          services.stage_input_files(INPUT_FILES)
        except Exception, e:
          print 'Error in call to stage_input_files()' , e
          services.error('Error in call to stage_input_files()')
          raise Exception, 'Error in call to stage_input_files()'
            
    # Copy genray.in_<suffix> to generic file name -> genray.in if there is
    # a suffix
        try:
          suffix = self.INPUT_SUFFIX
          have_suffix = True
        # If suffix is not empty put an underscore in front of it.
          if len(suffix) > 0:
              print 'genray INPUT_SUFFIX = ', suffix      
              suffix = '_' + suffix
        # If suffix is empty you don't really have one
          else:
              have_suffix = False
        except:
          have_suffix = False
      
      # If there is a non-empty suffix, copy to generic filename 'genray.in'
        if have_suffix:      
          try:
              shutil.copyfile('genray.in' + suffix, 'genray.in')
          except IOError, (errno, strerror):
              print 'Error copying file %s to %s' % ('genray.in' + suffix, 
              'genray.in'), strerror
              services.error('Error copying genray.in_<suffix> -> genray.in')
              raise Exception, 'Error copying genray.in_<suffix> -> genray.in'

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

    # Copy current dql file to generic name -> genray.nc
        try:
            shutil.copyfile(cur_dql_file, 'genray.nc')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % (cur_eqdsk_file, 'genray.nc'), strerror
            services.error('Error copying cur_dql_file -> genray.nc')
            raise Exception, 'Error copying cur_dql_file -> genray.nc'

        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_genray_input')

        rfmode = self.RFMODE
        isource_string = self.ISOURCE_STRING
        genraynml = self.GENRAYNML
        adj_read = self.ADJ_READ
        ps_add_nml = self.PS_ADD_NML

    # Call prepare_input - init
        print 'rf_genray: calling prepare_input init'        
        log_file = open('log_prepare_genray_input_init', 'w')
        mode = 'init'
        command = prepare_input_bin + ' ' + mode + ' ' +  rfmode + ' ' +\
        isource_string + ' ' + genraynml + ' ' + adj_read + ' ' + ps_add_nml
        
        retcode = subprocess.call(command.split(), stdout = log_file,\
                                  stderr = subprocess.STDOUT)
        if (retcode != 0):
            print 'Error executing genray init ', prepare_input_bin
            services.error('Error executing genray init')
            raise Exception, 'Error executing genray init'

    # Copy generic cur_state.cdf -> current plasma state file
        try:
            shutil.copyfile('cur_state.cdf', cur_state_file)
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % ('cur_state.cdf', cur_state_file, strerror)
            services.error('Error copying cur_state.cdf -> cur_state_file')
            raise Exception, 'Error copying cur_state.cdf -> cur_state_file'

    # Update plasma state files in plasma_state work directory
        try:
          services.update_state()
        except Exception, e:
          print 'Error in call to update_state()', e
          services.error('Error in call to update_state()')
          raise Exception, 'Error in call to update_state()'
     
    # Archive output files
    # N.B.  prepare_genray_input in init mode does not produce a complete set 
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
      print 'genray.restart() called'

      if (self.services == None) :
         print 'Error in genray: step (): No services'
#Removed by PB:    services.error('Error in genray: step (): No services')
         raise Exception, 'Error in genray: step (): No services'
      services = self.services

    # Get restart files listed in config file.        
      try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
      except Exception, e:
            print 'Error in call to get_restart_files()' , e
            services.error('genray: error in call to get_restart_files()')
            raise Exception, 'genray: error in call to get_restart_files()'

    # Get global configuration parameters
      try:
            self.plasma_state_file = services.get_config_param('CURRENT_STATE')
            self.eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            self.dql_file = services.get_config_param('CURRENT_EQDSK')
      except:
            print 'genray restart: error in getting config parameters'
            services.error('genray restart: error in getting config parameters')
            raise Exception, 'genray restart: error in getting config parameters'

      return 0
        
# ------------------------------------------------------------------------------
#
# STEP function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'genray.step() called'

        if (self.services == None) :
           print 'Error in genray: step (): No services'
           raise Exception, 'Error in genray: step (): No services'
        services = self.services

        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_genray_input')
        process_output_bin  = os.path.join(self.BIN_PATH, 'process_genray_output_mcmd')
        genray_bin = os.path.join(self.BIN_PATH, 'genray')

    # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception, e:
          print 'Error in call to stage_state()' , e
          services.error('Error in call to stage_state()')
          raise Exception, 'Error in call to stage_state()'

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_dql_file = services.get_config_param('CURRENT_DQL') 
       
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

        rfmode = self.RFMODE
        isource_string = self.ISOURCE_STRING
        genraynml = self.GENRAYNML
        adj_read = self.ADJ_READ
        ps_add_nml = self.PS_ADD_NML

# Check if LHRF power is zero (or effectively zero).  If true don't run Genray

        ps = NetCDFFile(cur_state_file, 'r')
        power_lh = ps.variables['power_lh'].getValue()[0]
        ps.close()
        print 'power = ', power_lh
        if(power_lh > 1.0E+03):

    # Call prepare_input - step
          print 'rf_genray step: calling prepare_input'
          
          log_file = open('log_prepare_genray_input_step', 'w')
          mode = 'step'
          command = prepare_input_bin + ' ' + mode + ' ' +  rfmode + ' ' +\
          isource_string + ' ' + genraynml + ' ' + adj_read + ' ' + ps_add_nml
          
          retcode = subprocess.call(command.split(), stdout = log_file,\
                                  stderr = subprocess.STDOUT)
          if (retcode != 0):
              print 'Error executing genray: ', prepare_input_bin
              services.error('Error executing genray prepare_input')
              raise Exception, 'Error executing genray prepare_input'
        
#     Launch genray - N.B: Path to executable is in config parameter GENRAY_BIN
          print 'rf_genray: launching genray'
          cwd = services.get_working_dir()
          task_id = services.launch_task(self.NPROC, cwd, self.GENRAY_BIN, logfile='log.genray')
          retcode = services.wait_task(task_id)
          if (retcode != 0):
              print 'Error executing command: ', genray_bin
              services.error('Error executing genray')
              raise Exception, 'Error executing genray'

    # Copy current genray.nc file to Plasma State -> CURRENT_DQL
          try:
              shutil.copyfile('genray.nc', cur_dql_file)
          except IOError, (errno, strerror):
              print 'Error copying file %s to %s' % ('genray.nc', cur_dql_file), strerror
              services.error('Error copying genray.nc -> cur_dql_file')
              raise Exception, 'Error copying genray.nc -> cur_dql_file'

    # Call process_output - step
          print 'rf_genray step: calling process_output'
        
          log_file = open('log_process_genray_output', 'w')
          mode = 'step'
          command = process_output_bin + ' ' +  rfmode + ' ' + isource_string
        
          retcode = subprocess.call(command.split(), stdout = log_file,\
                                  stderr = subprocess.STDOUT)                                  
          if (retcode != 0):
              print 'Error executing genray init ', process_output_bin
              services.error('Error executing genray process_output')
              raise Exception, 'Error executing genray process_output'
        
    # Copy generic genray partial plasma state file -> LH_cur_state_file
          try:
              partial_file = cwd + '/RF_LH_' + cur_state_file
              shutil.copyfile('RF_GENRAY_PARTIAL_STATE', partial_file )
          except IOError, (errno, strerror):
              print 'Error copying file %s to %s' % ('cur_state.cdf', cur_state_file, strerror)
              raise


# Merge partial plasma state containing updated LH data
          try:
              services.merge_current_state(partial_file, logfile='log.update_state')
              print 'merged GENRAY plasma state data ', partial_file
          except Exception, e:
              print 'Error in call to merge_current_state(' , partial_file, ')'
              self.services.error('Error in call to merge_current_state')
              raise Exception, 'Error in call to merge_current_state'
            
    # Archive output files
          try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
          except Exception, e:
            print 'Error in call to stage_output_files()', e
            services.error('Error in call to stage_output_files()')
            raise Exception, 'Error in call to stage_output_files()'
# ptb tried this
    # Update plasma state files in plasma_state work directory
          try:
            services.update_state()
          except Exception, e:
            print 'Error in call to update_state()', e
            services.error('Error in call to update_state()')
            raise Exception, 'Error in call to update_state()'
# ptb tried this
        
#     rename the log file so that it is not appended next step
#        os.rename('log.genray', this_logfile)

          return 0

        return 0  # return on "zero" LHRF power condition
# ------------------------------------------------------------------------------
#
# checkpoint function
# Saves plasma state files to restart directory
#
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
      print 'genray.checkpoint() called'
      if (self.services == None) :
         print 'Error in genray: checkpoint(): No services'
         services.error('Error in genray: checkpoint(): No services')
         raise Exception, 'Error in genray: checkpoint(): No services'
      services = self.services
      services.save_restart_files(timestamp, self.RESTART_FILES)

      return 0
        
# ------------------------------------------------------------------------------
#
# finalize function
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print 'genray.finalize() called'

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

