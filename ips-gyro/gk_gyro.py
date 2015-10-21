#! /usr/bin/env python

"""
gk_gyro.py version 0.0 Batchelor 2-24-2015

GK_GYRO component script to drive GYRO binary code. 

Present version 0.0 invokes the bash shell script "gyro" rather than launching the
gyro binary directly.

"""
# ------------------------------------------------------------------------------
#
#
# Workflow:
# 1) 'init' function:
#       Stage input files: batch.src_<suffix.>, input.gyro_<suffix>
#       Copy to generic input file names -> batch.src, input.gyro

# 2) 'step' function:
#       If parameters are programmed from config file modify input.gyro
#       Launch gyro executable, redirect output to new log file
#
# ------------------------------------------------------------------------------

import sys
import os
import subprocess
import getopt
import shutil
import string
from  component import Component
from netCDF4 import *
from numpy import *

#---------------------------------------------------------------------------------------
# Some functions for editing input file
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# Open an input file and return the lines
def get_lines(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    return lines

#---------------------------------------------------------------------------------------
# Open an output file and write lines into it
def put_lines(filename, lines):
    file = open(filename, 'w')
    file.writelines(lines)
    file.close()

#---------------------------------------------------------------------------------------
# edit gyro input file
def edit_gyro_input_file_0(lines, var, val):

    # Find the line in the input containing 'var='
    var_line_number = -1
    for i in range(len(lines)):
        line = lines[i]
        if '=' in line:
            split_line = line.split('=')
            #print 'split_line = ', split_line
            if (split_line[0].strip()).lower() == var.lower():
                var_line_number = i

    if var_line_number == -1:
        message = 'Could not find variable ', var, ' in input file lines'
        print message
        raise Exception(message)
        
    # Get rid of newline if present.
    
    if lines[var_line_number][-1] == '\n':
        lines[var_line_number] = lines[var_line_number][:-1]
            
    # Insert line with new values at lines[i]
    lines[var_line_number] = lines[var_line_number].split('=')[0] + '=' + str(val) + '\n'
    print 'New ', lines[var_line_number]
        
    return lines

#----------------------------------------------------------------------------------------------        
# gk_gyro component
#----------------------------------------------------------------------------------------------

class gk_gyro(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

#----------------------------------------------------------------------------------------------        
# Utility functions
#----------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# edit gyro input file
    def edit_gyro_input_file(self, lines):
    
        # Go thorough lines in input file
        for i in range(len(lines)):
            line = lines[i]
            if '=' in line:

                # Get rid of newline if present
                if line[-1] == '\n':
                    line = line[:-1] 
                           
                split_line = line.split('=')            
                var = split_line[0].strip()
            
                # If this variable is programmed from config file update it
                if hasattr(self, var):
                    print 'edit_gyro_input_file: Changing'
                    print line
                    print 'to'
                    val = getattr(self,var)
                    lines[i] = var + '=' + str(val) + '\n'
                    print lines[i]
        
        return lines


# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp=0):
        print 'gk_gyro.init() called'

        if (self.services == None) :
           print 'Error in gk_gyro: init(): No services'
           raise Exception, 'Error in gk_gyro: init(): No services'

        services = self.services
        INPUT_FILES = self.INPUT_FILES
        
    # If batch.src is to be programmed from config file put code here

    # Get input files  
        try:
          services.stage_input_files(INPUT_FILES)
        except Exception, e:
          print 'Error in call to stageInputFiles()' , e
          services.error('Error in call to stageInputFiles()')
          raise Exception, 'Error in call to stageInputFiles()'

    # Copy input.gyro_<suffix> to generic file name -> input.gyro if there is
    # a suffix
        try:
          suffix = self.INPUT_SUFFIX
          have_suffix = True
        # If suffix is not empty put an underscore in front of it.
          if len(suffix) > 0:
              print 'gyro INPUT_SUFFIX = ', suffix      
              suffix = '_' + suffix
        # If suffix is empty you don't really have one
          else:
              have_suffix = False
        except:
          have_suffix = False

      # If there is a non-empty suffix, copy batch script to generic filename 'batch.src'
#         if have_suffix:      
#           try:
#               shutil.copyfile('batch.src' + suffix, 'batch.src')
#           except IOError, (errno, strerror):
#               print 'Error copying file %s to %s' % ('batch.src' + suffix, 
#               'batch.src'), strerror
#               services.error('Error copying batch.src_<suffix> -> batch.src')
#               raise Exception, 'Error copying batch.src_<suffix> -> batch.src'
#           log_file = open('log_sub_proc_call_init', 'w')
#           err_file = open('err_sub_proc_call_init', 'w')
#           command = 'chmod -v a+x batch.src'
#           retcode = subprocess.call(command.split(), shell = True, stdout = log_file,\
#                               stderr = err_file) 
#                               #stderr = subprocess.STDOUT) 
      
      # If there is a non-empty suffix, copy to generic filename 'input.gyro'
        if have_suffix:      
          try:
              shutil.copyfile('input.gyro' + suffix, 'input.gyro')
          except IOError, (errno, strerror):
              print 'Error copying file %s to %s' % ('input.gyro' + suffix, 
              'input.gyro'), strerror
              services.error('Error copying input.gyro_<suffix> -> input.gyro')
              raise Exception, 'Error copying input.gyro_<suffix> -> input.gyro'

    # If parameters are programmed from config file modify parameters in input.gyro
        programming = False
        try:
            programming = self.PROGRAMMING
            print 'programming = ', programming
        except:
            pass
        
        if programming == 'True':
            in_file_name = 'input.gyro'
            print 'in_file_name = ', in_file_name
            input_lines = get_lines(in_file_name)
            
#             try:
#                 DLNNDR_ELECTRON = self.DLNNDR_ELECTRON
#                 print "DLNNDR_ELECTRON = ", DLNNDR_ELECTRON, " programmed from config file"
#             except:
#                 message = 'gk_gyro init: error getting DLNNDR_ELECTRON from config file'
#                 print message
#                 services.error(message)
#                 raise Exception, message
#                         
#             var = 'DLNNDR_ELECTRON'
#             value = DLNNDR_ELECTRON
#             lines = edit_gyro_input_file(input_lines, var, value)

            lines = self.edit_gyro_input_file(input_lines)
            
            put_lines(in_file_name, lines)

    # Archive output files
    # Nothing to do, init mode for this component does not produce output files.

        return 0
# ------------------------------------------------------------------------------
#
# RESTART function
# Gets restart files from restart directory
# Loads the global configuration parameters from the config file
#
# ------------------------------------------------------------------------------
        
    def restart(self, timeStamp):
      print 'gyro.restart() called'

      if (self.services == None) :
         print 'Error in gyro: restart(): No services'
         raise Exception, 'Error in gyro: restart(): No services'
      services = self.services

    # Get restart files listed in config file.        
      try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
      except Exception, e:
            print 'Error in call to get_restart_files()' , e
            services.error('gyro: error in call to get_restart_files()')
            raise Exception, 'gyro: error in call to get_restart_files()'

      return 0
        
# ------------------------------------------------------------------------------
#
# STEP function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'gyro.step() called'

        if (self.services == None) :
           print 'Error in gyro: step (): No services'
           raise Exception, 'Error in gyro: step (): No services'
        services = self.services


    # If parameters are to be programmed from config file for multiple steps, put 
    # coding here
        
    # Launch gyro
        log_file = open('log_sub_proc_call', 'w')
#        command = './batch.src_gyro_test'
        command = './batch.src'
        retcode = subprocess.call(command.split(), shell = True, stdout = log_file,\
                              stderr = subprocess.STDOUT)
        if (retcode != 0):
          print 'Error executing gyro: ', command
          services.error('Error executing gyro '+ command)
          raise Exception, 'Error executing gyro '+ command

    # Process output

    # Archive output files
#         try:
#           services.stage_output_files(timeStamp, self.OUTPUT_FILES)
#         except Exception, e:
#           print 'Error in call to stage_output_files()', e
#           services.error('Error in call to stage_output_files()')
#           raise Exception, 'Error in call to stage_output_files()'
        
        return 0

# ------------------------------------------------------------------------------
#
# checkpoint function
# Saves plasma state files to restart directory
#
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
      print 'gyro.checkpoint() called'
      if (self.services == None) :
         print 'Error in gyro: checkpoint(): No services'
         services.error('Error in gyro: checkpoint(): No services')
         raise Exception, 'Error in gyro: checkpoint(): No services'
      services = self.services
      services.save_restart_files(timestamp, self.RESTART_FILES)

      return 0
        
# ------------------------------------------------------------------------------
#
# finalize function
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print 'gyro.finalize() called'
        return 0

