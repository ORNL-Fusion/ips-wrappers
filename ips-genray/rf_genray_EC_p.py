#! /usr/bin/env python

"""
Programmable version: Batchelor (6/2/2013)

This version allows the aiming angles and power for multiple launchers to be set and
changed over time from within the simulation config file.  If this option is used, any 
programming of the launcher power from the EPA component is over-ridden.  This behavior 
is triggered by the presence a variable N_LAUNCHERS_PROGRAMMED in the [rf_genray] section 
of the simulation config file.  If N_LAUNCHERS_PROGRAMMED = 0, or is absent, the launcher
powers are taken from the current plasma state file.  
If > 0 the component looks for space delimited lists of time points for parameter changes 
and the associated parameters.  The config lists should be named: 
LAUNCHER1_TIMES,LAUNCHER1_alfast, LAUNCHER1_betast, LAUNCHER1_powtot_MW, 
LAUNCHER2_TIMES,LAUNCHER2_alfast, LAUNCHER2_betast, LAUNCHER2_powtot_MW,  etc

The parameter changes take effect the first time the simulation time at
the beginning of a time step (ps%t0) is >= LAUNCHERx_TIME then LAUNCHERx parameters are
changed.  LAUNCHERx_TIMES don't have to match simulation time step beginnings, but the 
changes won't take place until the beginning of the next time step.

Note: For ease of typing, enter powers in config file in MW, they get converted to Watts
in this component.  

The mechanism for changing the power is to modify the genray.in just before launching
the genray_EC step.  N_LAUNCHERS_PROGRAMMED must match the number of launchers found in
the input genray.in file.
"""
# 11/25/2014 DBB
# In many of the components we have coding to allow the user to add a suffix to an input 
# file name and then copy it to a generic file name, for example genray.in_myFile ->
# genray.in.  However Bob Harvey has this functionality in the prepare_genray_input.f90
# code.  He has "genraynml" as command line arg to prepare_genray_input, which then
# writes to the generic name genray.in.  So to avoid confusion I have removed the coding
# providing the suffix processing as used in other components.

# Also found a bug in the exception handling coding for file copies, which I fixed 
# (I think).

# EC version 0.0 1/10/2011
# This version distinguishes between the different RF components that GENRAY can 
# implement.  In this case the EC component.  So far the only place this affects
# is in the writing the partial plasma state that the framework needs to merge.
#
# Note: To merge plasma states the IPS framework expects the component to
# produce a partial state file with the component-specific name:
# <component>_<curr_state_filename>.  GENRAY works for EC, EC, and IC, i.e.
# it can implement several different components.  Therefore process_genray_output_mcmd.f90
# writes a partial state file with generic name 'RF_GENRAY_PARTIAL_STATE' and delegetes
# to this python component script the copying of this to the proper 
# component-specific update file name.  In this case RF_EC_<curr_state_filename>.

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
from  component import Component
from Numeric import *
from Scientific.IO.NetCDF import *

#---------------------------------------------------------------------------------------
# Some functions for editing a fortran namelist file
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
# Edit fortran namelist file
#---------------------------------------------------------------------------------------

# This function allows changing the variable values in a fortran namelist file. 
# Arguments:
#
# lines = List of character strings representing the of lines from the input namelist file.
#    For example obtained from the get_lines() function above
#
# var = String containing the name of the variable to be changed.
#
# value = String or list of strings containing data to go into var.  The editing that can
#     be done is pretty rudimentary.  If var is a vector then value must contain the whole
#     new vector.  For example it doesn't support changing an array element or array 
#     slice.
#
# separator = Optional argument specifying the separator between array elements.  The
#     the default is comma
#

def edit_nml_file(lines, var, values, separator = ','):

    # Find the line in the namelist containing 'var = ' and get rid of newline if present.
    var_line_number = -1
    for i in range(len(lines)):
        line = lines[i]
        if '=' in line:
            split_line = line.split('=')
            #print 'split_line = ', split_line
            if (split_line[0].strip()).lower() == var.lower():
                var_line_number = i

    if var_line_number == -1:
        message = 'Could not find variable ', var, ' in namelist lines'
        print message
        raise Exception(message)
        
    #print 'var_line_number = ', var_line_number
    #print 'lines[var_line_number] = ', lines[var_line_number]
    
    if lines[var_line_number][-1] == '\n':
        lines[var_line_number] = lines[var_line_number][:-1]

    # Try to find out how many lines of text go with this variable.  So find the next
    # line with either an '=' sign or a '/' signifying the end of the namelist group.
    var_lines = 1
    test = False
    while test == False:
        next_iine_no = var_line_number + var_lines
        next_line = lines[next_iine_no]
        if '=' in next_line:   # Could get fooled by = in a quoted string
            test = True
            eq_index = next_line.find('=') # so check if quote before =
            single_quote_index = next_line.find("'")
            if single_quote_index > -1 and single_quote_index < eq_index:
                test = False
            double_quote_index = next_line.find('"')
            if double_quote_index > -1 and double_quote_index < eq_index:
                test = False
        elif next_line[-1] == '/':  # At end of line means end of group
            test = True
        else:
            var_lines += 1
            
    # Insert line with new values at lines[i]
    lines[var_line_number] = lines[var_line_number].split('=')[0] + ' = '
    for val in values:
        lines[var_line_number] = lines[var_line_number] + str(val) + ', '
    lines[var_line_number] = lines[var_line_number] + '\n'
    print 'New ', lines[var_line_number]
        
    return lines[:var_line_number + 1] + lines[var_line_number + var_lines:]



class genray_EC(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

#----------------------------------------------------------------------------------------------        
# Utility functions
#----------------------------------------------------------------------------------------------


    def  piecewise_constant(self, x, x_points, y_points):
    
    # Function that is piecewise constant within zones defined by x_points.
    # x = function argument
    # x_points = monotonically increasing list of x values
    # y_points = list of function values at the x_points, must be same length as x_points
    # 
    # If x < x_points[0] returns zero
    # If x >= x_points[i] - eps and x < x_points[i+1] returns y_points[i]
    # If x > x_points[-1] returns y_points[-1]
    #
    # Since floating equality is always problematic, x is considered to be in the ith
    # zone if it is within -eps of x[i] where eps is small.

        eps = 10.0**(-10)
        
        if len(x_points) != len(y_points):
            raise Exception('piecewise_constant: len(x_points != len(y_points)')
              
        if x < x_points[0]:
            return 0.0
        elif x > x_points[-1] - eps:
            return y_points[-1]
        else:
            j_up = 1
            while x > x_points[j_up] - eps:  # Find the index of the x_point just above x.
                j_up = j_up + 1
            return y_points[j_up-1]      
        
#----------------------------------------------------------------------------------------------

    def time_points_and_parameters(self, n_launchers):
    
    # Does initializations for set_genray_EC_parameters() defined below.
    # Gets the ECH launcher programming parameters from the [EC_launcher] section of the 
    # config file:
    # LAUNCHERx_TIME, LAUNCHERx_alfast, LAUNCHERx_betast, LAUNCHERx_powtot_MW
    # Does a little checking for consistency and returns a 4 element list of lists
    # [times_list, alfast_list, betast_list, powtot_list]
    # Where times_list and  parameers lists are length = n_launchers and contain
    # lists with the time points and associated parameters for each of the launchers.  
    # Thus the shape of the returned list is (4,n_launchers, len(launcherx_TIME)) for 
    # x = 1 to n_launchers. It also converts powers to Watts so as to be consistent with 
    # plasma state.
    
        times_list = []
        alfast_list = []
        betast_list = []
        powtot_list = []
        for i in range(n_launchers):
            launcher_name = 'LAUNCHER' + str(i+1)

            # get time points for this launcher
            name = launcher_name + '_TIMES'
            time_expr = 'self.' + name + '.split()'
            try:
                str_time_points = eval(time_expr)
            except:
                message = 'error in getting EC_launcher config parameters ' + name
                print message
                self.services.error(message)
                raise Exception, message
            print name, ' = ', str_time_points

            # get alfasts for this launcher
            name = launcher_name + '_alfast'
            expr = 'self.' + name + '.split()'
            try:
                str_alfast = eval(expr)
            except:
                message = 'error in getting EC_launcher config parameters ' + name
                print message
                self.services.error(message)
                raise Exception, message
            print name, ' = ', str_alfast
 
            # get betasts for this launcher
            name = launcher_name + '_betast'
            expr = 'self.' + name + '.split()'
            try:
                str_betast = eval(expr)
            except:
                message = 'error in getting EC_launcher config parameters ' + name
                print message
                self.services.error(message)
                raise Exception, message
            print name, ' = ', str_betast

            # get powtots for this launcher
            name = launcher_name + '_powtot_MW'
            expr = 'self.' + name + '.split()'
            try:
                str_powtot_MW = eval(expr)
            except:
                message = 'error in getting EC_launcher config parameters ' + name
                print message
                self.services.error(message)
                raise Exception, message
            print name, ' = ', str_powtot_MW
                
            # check that times and parameters have same length
            if len(str_time_points) != len(str_powtot_MW):
                message = 'error: ', launcher_name, \
                ' time points and parameters have different lengths' 
                print message
                self.services.error(message)
                raise Exception, message
                
            # Convert to float, convert powers to Watts, and append to full lists
            times_list.append([float(x) for x in str_time_points])
            alfast_list.append([float(x) for x in str_alfast])
            betast_list.append([float(x) for x in str_betast])
            powtot_list.append([10.0**6*float(x) for x in str_powtot_MW])
            
        return [times_list, alfast_list, betast_list, powtot_list]

    #-----------------------------------------------------------------------------
    def set_genray_EC_parameters(self, input_file,  t, times_parameters_list):
    

        times_list = times_parameters_list[0]
        alfast_list = times_parameters_list[1]
        betast_list = times_parameters_list[2]
        powtot_list = times_parameters_list[3]

        # Get input_file
        file = open(input_file, 'r')
        nml= file.readlines()
        file.close()
        
        # get parameters for each launcher at this time
        n_launchers = len(times_list)
        alfast_t = []
        betast_t = []
        powtot_t = []
        for i in range(n_launchers):
            alfast_launcher_t = self.piecewise_constant(t, times_list[i], alfast_list[i])
            alfast_t.append(alfast_launcher_t)
            betast_launcher_t = self.piecewise_constant(t, times_list[i], betast_list[i])
            betast_t.append(betast_launcher_t)
            powtot_launcher_t = self.piecewise_constant(t, times_list[i], powtot_list[i])
            powtot_t.append(powtot_launcher_t)

        # change parameters in namelist file
        lines = get_lines(input_file)
        lines = edit_nml_file(lines, 'alfast', alfast_t, separator = ',')
        lines = edit_nml_file(lines, 'betast', betast_t, separator = ',')
        lines = edit_nml_file(lines, 'powtot', powtot_t, separator = ',')
        put_lines(input_file, lines)

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp=0):
        print 'genray.init() called'

        services = self.services
        global programming, times_parameters_list
        
    # Get global configuration parameters
        try:
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            #cql_file = services.get_config_param('CURRENT_CQL')
        except:
            print 'rf_genray_EC: error getting config parameters CURRENT_STATE CURRENT_EQDSK'
            services.error('rf_genray_EC: error getting config parameters CURRENT_STATE CURRENT_EQDSK')
            raise Exception, 'rf_EC_gneray: error getting config parameters CURRENT_STATE CURRENT_EQDSK'

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
            print 'rf_genray_EC init: error getting genray-specific config parameters'
            services.error('rf_genray_EC: error getting genray-specific\
            config parameters')
            raise Exception, 'rf_genray_EC: error getting genray-specific\
            config parameters'

        # Get [rf_genray_EC] programming configuration parameters, if present
        n_launchers = 0
        programming = False
        try:
            n_launchers = int(self.N_LAUNCHERS_PROGRAMMED)
        except:
            print '\nCould not get launcher programming from config file'

        if n_launchers > 0: 
            programming = True
            times_parameters_list = self.time_points_and_parameters(n_launchers)
        else:
            print '\nUsing ECH launcher parmeters in plasma state'


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
            
    # Copy current plasma state file to generic name -> cur_state.cdf
        try:
            shutil.copyfile(cur_state_file, 'cur_state.cdf')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % (cur_state_file, 'cur_state.cdf'), strerror
            services.error('Error copying cur_state_file -> cur_state.cdf')
            raise Exception, 'Error copying cur_state_file -> cur_state.cdf'

    # Copy current eqdsk file to generic name -> eqdsk
        try:
            shutil.copyfile(cur_eqdsk_file, 'eqdsk')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % (cur_eqdsk_file, 'eqdsk'), strerror
            services.error('Error copying cur_eqdsk_file -> eqdsk')
            raise Exception, 'Error copying cur_eqdsk_file -> eqdsk'

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

    # If parameters are programmed from config file set parameters in genray.in
        if programming == True:
            # Get t0 from plasma state
            ps = NetCDFFile(cur_state_file, 'r')
            t0 = ps.variables['t0'].getValue()
            ps.close()
            # set parameters those for time = t0
            self.set_genray_EC_parameters("genray.in", t0, times_parameters_list)

    # Copy generic cur_state.cdf -> current plasma state file
        try:
            shutil.copyfile('cur_state.cdf', cur_state_file)
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % ('cur_state.cdf', cur_state_file), strerror
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
      print 'rf_genray_EC_p'

      if (self.services == None) :
         print 'Error in genray: step (): No services'
         raise Exception, 'Error in genray: step (): No services'
      services = self.services
      global programming, times_parameters_list

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
      except:
            print 'genray restart: error in getting config parameters'
            services.error('genray restart: error in getting config parameters')
            raise Exception, 'genray restart: error in getting config parameters'

        # Get [rf_genray_EC] programming configuration parameters, if present
      n_launchers = 0
      programming = False
      try:
          n_launchers = int(self.N_LAUNCHERS_PROGRAMMED)
      except:
          print '\nCould not get launcher programming from config file'

      if n_launchers > 0: 
          programming = True
          times_parameters_list = self.time_points_and_parameters(n_launchers)
      else:
          print '\nUsing ECH launcher parmeters in plasma state'


      return 0
        
# ------------------------------------------------------------------------------
#
# STEP function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'genray.step() called'
        print 'rf_genray_EC_p'

        if (self.services == None) :
           print 'Error in genray: step (): No services'
           raise Exception, 'Error in genray: step (): No services'
        services = self.services
        global programming, times_parameters_list

        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_genray_input')
        process_output_bin  = os.path.join(self.BIN_PATH, 'process_genray_output_mcmd')
        genray_bin = os.path.join(self.BIN_PATH, 'genray')

    # Copy plasma state files over to working directory
        try:
          services.stage_plasma_state()
        except Exception, e:
          print 'Error in call to stage_plasma_state()' , e
          services.error('Error in call to stage_plasma_state()')
          raise Exception, 'Error in call to stage_plasma_state()'

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        
    # Copy current plasma state file to generic name -> cur_state.cdf
        try:
            shutil.copyfile(cur_state_file, 'cur_state.cdf')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % (cur_state_file, 'cur_state.cdf'), strerror
            raise

    # Copy current eqdsk file to generic name -> eqdsk
        try:
            shutil.copyfile(cur_eqdsk_file, 'eqdsk')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % (cur_eqdsk_file, 'eqdsk'), strerror
            services.error('Error copying cur_eqdsk_file -> eqdsk')
            raise Exception, 'Error copying cur_eqdsk_file -> eqdsk'

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
        if (retcode != 0):
            print 'Error executing genray: ', prepare_input_bin
            services.error('Error executing genray prepare_input')
            raise Exception, 'Error executing genray prepare_input'

    # If parameters are programmed from config file set parameters in genray.in
        if programming == True:
            # Get t0 from plasma state
            ps = NetCDFFile(cur_state_file, 'r')
            t0 = ps.variables['t0'].getValue()
            ps.close()
            # set parameters those for time = t0
            self.set_genray_EC_parameters("genray.in", t0, times_parameters_list)
        
#     Launch genray - N.B: Path to executable is in config parameter GENRAY_BIN
        print 'rf_genray: launching genray'
        cwd = services.get_working_dir()
        task_id = services.launch_task(self.NPROC, cwd, self.GENRAY_BIN, logfile='log.genray')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            print 'Error executing command: ', genray_bin
            services.error('Error executing genray')
            raise Exception, 'Error executing genray'

    # Call process_output - step
        print 'rf_genray step: calling process_output'
        
        log_file = open('log_process_genray_output', 'w')
        mode = 'step'
        command = process_output_bin + ' ' +  rfmode + ' ' + isource_string

        print 'running', command
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
            event_comment =  command)
        
        retcode = subprocess.call(command.split(), stdout = log_file,\
                                  stderr = subprocess.STDOUT)                                  
        if (retcode != 0):
            print 'Error executing genray init ', process_output_bin
            services.error('Error executing genray process_output')
            raise Exception, 'Error executing genray process_output'
        
    # Copy generic genray partial plasma state file -> EC_cur_state_file
        try:
            partial_file = cwd + '/RF_EC_' + cur_state_file
            shutil.copyfile('RF_GENRAY_PARTIAL_STATE', partial_file )
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % ('cur_state.cdf', cur_state_file), strerror
            raise


# Merge partial plasma state containing updated EC data
        try:
            services.merge_current_plasma_state(partial_file, logfile='log.update_state')
            print 'merged GENRAY plasma state data ', partial_file
        except Exception, e:
            print 'Error in call to merge_current_plasma_state(' , partial_file, ')'
            self.services.error('Error in call to merge_current_plasma_state')
            raise Exception, 'Error in call to merge_current_plasma_state'
            
    # Archive output files
        try:
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
          print 'Error in call to stage_output_files()', e
          services.error('Error in call to stage_output_files()')
          raise Exception, 'Error in call to stage_output_files()'
        
#     rename the log file so that it is not appended next step
#        os.rename('log.genray', this_logfile)

        return 0

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

