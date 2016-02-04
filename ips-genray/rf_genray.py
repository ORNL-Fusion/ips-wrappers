#! /usr/bin/env python

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
from  component import Component
from Numeric import *
from Scientific.IO.NetCDF import *

class genray(Component):
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

        
    # Copy input files to working directory
        try:
            services.stage_input_files(self.INPUT_FILES)
        except Exception, e:
            print 'Error in call to stageInputFiles()' , e
            
    # Copy genray.in_<suffix> to generic file name -> genray.in        
        suffix = self.INPUT_SUFFIX
        # check if the suffix is empty.  If not put an underscore in front of it.
        if len(suffix) > 0: suffix = '_' + suffix
        print 'INPUT_SUFFIX = ', suffix
        try:
            shutil.copyfile('genray.in' + suffix, 'genray.in')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % ('genray.in' + suffix, 
            'genray.in', strerror)
            raise


        workdir = services.get_working_dir()

    # Copy needed state files to working directory
        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'Error in call to services.stage_plasma_state()', e

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        #cql_file = services.get_config_param('CURRENT_CQL')

    # Copy current plasma state file to generic name -> cur_state.cdf
        try:
            shutil.copyfile(cur_state_file, 'cur_state.cdf')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % (cur_state_file, 'cur_state.cdf', strerror)
            raise

        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_genray_input')
        files = self.INPUT_FILES.split()

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
                                  
        # retcode = subprocess.call([prepare_input_bin, 'init', rfmode, isource_string, \
#                                    genraynml, adj_read, ps_add_nml])
        if (retcode != 0):
            print 'Error executing genray init ', prepare_input_bin
            raise

    # Copy generic cur_state.cdf -> current plasma state file
        try:
            shutil.copyfile('cur_state.cdf', cur_state_file)
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % ('cur_state.cdf', cur_state_file, strerror)
            raise

#     Update plasma state
        try:
            services.update_plasma_state()
        except Exception, e:
            print 'Error in call to updatePlasmaState()', e
            raise
            
    # Stage output files
    
        # Touch output files so framework won't complain about missing files
        for file in self.OUTPUT_FILES.split():
            print file
            subprocess.call(['touch',file])
            
        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
            print 'Error in call to stage_output_files()', e
            raise

        return

        
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
         raise Exception, 'Error in genray: step (): No services'
      services = self.services
      workdir = services.get_working_dir()

    # Get restart files listed in config file.        
      try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
      except Exception, e:
            print 'Error in call to get_restart_files()' , e
            raise

    # Get global configuration parameters
      try:
            self.plasma_state_file = services.get_config_param('CURRENT_STATE')
            self.eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            self.toric_log = os.path.join(workdir, 'log.genray')
      except:
            print 'genray restart: error in getting config parameters'
            self.services.error('error in getting config parameters')
            raise Exception, 'error in getting config parameters'

      return 0
        
# ------------------------------------------------------------------------------
#
# STEP function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'genray.step() called'

        if (self.services == None) :
            print 'Error in genray: step () no self.services'
            raise Exception, 'Error in genray: step () no self.services'
            
        services = self.services

        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_genray_input')
        process_output_bin  = os.path.join(self.BIN_PATH, 'process_genray_output')
        genray_bin = os.path.join(self.BIN_PATH, 'genray')

    # Copy needed state files to working directory
        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'Error in call to services.stage_plasma_state()', e

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        
    # Copy current plasma state file to generic name -> cur_state.cdf
        try:
            shutil.copyfile(cur_state_file, 'cur_state.cdf')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % (cur_state_file, 'cur_state.cdf', strerror)
            raise Exception, 'Error copying current plasma state file to generic name'

        files = self.INPUT_FILES.split()

        rfmode = self.RFMODE
        isource_string = self.ISOURCE_STRING
        genraynml = self.GENRAYNML
        adj_read = self.ADJ_READ
        ps_add_nml = self.PS_ADD_NML

    # Call prepare_input - step
        print 'rf_genray step: calling prepare_input'
        
        log_file = open('log_prepare_genray_input_step', 'w')
        mode = 'step'
        command = prepare_input_bin + ' ' + mode + ' ' +  rfmode + ' ' +\
        isource_string + ' ' + genraynml + ' ' + adj_read + ' ' + ps_add_nml
        
        retcode = subprocess.call(command.split(), stdout = log_file,\
                                  stderr = subprocess.STDOUT)
                                  
        # retcode = subprocess.call([prepare_input_bin, 'init', rfmode, isource_string, \
#                                    genraynml, adj_read, ps_add_nml])
        if (retcode != 0):
            print 'Error executing genray init ', prepare_input_bin
            raise Exception, 'Error executing genray init '
        
#     setup file list so as to get antenna file name
        files = self.INPUT_FILES.split()
        
#     Launch genray - N.B: Path to executable is in config parameter GENRAY_BIN
        print 'rf_genray: launching genray'
        cwd = services.get_working_dir()
        task_id = services.launch_task(self.NPROC, cwd, self.GENRAY_BIN, logfile='log.genray')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            msg = 'Error executing command: ' + genray_bin
            print msg
            raise Exception, msg

    # Call process_output - step
        print 'rf_genray step: calling process_output'
        
        log_file = open('log_prepare_genray_output', 'w')
        mode = 'step'
        command = process_output_bin + ' ' +  rfmode + ' ' + isource_string
        
        retcode = subprocess.call(command.split(), stdout = log_file,\
                                  stderr = subprocess.STDOUT)                                  
        if (retcode != 0):
            msg = 'Error executing genray init ' + process_output_bin
            print msg
            raise Exception, msg

#     Call process_output and copy files over
#         print 'rf_genray: calling process_output step'
#         retcode = subprocess.call([process_output_bin, rfmode, isource_string])
#         if (retcode != 0):
#             print 'Error executing',  process_output_bin
#             sys.exit(1)
#         this_logfile = 'log.genray_'+'%s' %timeStamp
        
# Merge partial plasma state containing updated IC data
        try:
          partial_file = cwd + '/RF_IC_' + cur_state_file
          services.merge_current_plasma_state(partial_file, logfile='log.update_state')
          print 'merged TORIC plasma state data ', partial_file
        except Exception, e:
          print 'Error in call to merge_current_plasma_state(' , partial_file, ')'
          self.services.error('Error in call to merge_current_plasma_state')
          raise Exception, 'Error in call to merge_current_plasma_state'

        # Archive output files
        try:
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
          print 'Error in call to stage_output_files()', e
          self.services.error('Error in call to stage_output_files()')
          raise Exception, 'Error in call to stage_output_files()'
                    
#     rename the log file so that it is not appended next step
        os.rename('log.genray', this_logfile)

        return 0

# ------------------------------------------------------------------------------
#
# checkpoint function
# Saves plasma state files to restart directory
#
# ------------------------------------------------------------------------------

   def checkpoint(self, timestamp=0.0):
        print 'rf_genray.checkpoint() called'
        services = self.services
        services.save_restart_files(timestamp, self.RESTART_FILES)
        
        
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
