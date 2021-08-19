#! /usr/bin/env python

"""
Multi-frequency, config file programmable version: Batchelor (9/3/2018)
Supports both RFMODE = EC and LH.  Launcher power and launcher aiming programming is only 
implemented for ECH, so far.

"""

# Working notes:

# 9/3/2018 DBB
# This version is adapted from the previous version genray_EC_p.py.  It retains the
# programmability from the config file that was in genray_EC_p.py and the detection of zero
# RF power that was in the 6/16/2014 version but that seemed to get dropped from the
# genray_EC_p.py version.  The programming was implemented for ECH aiming angles, but no
# such k-spectrum programming is implemented for Lower Hybrid, yet.
# Switched NetCDF4.py from obsolete Scientific.IO.NetCDF.  Removed namelist editing routines
# to /wrappers/utilities/simple_file_editing_functions.py.

# 11/25/2014 DBB
# In many of the IPS components we have coding to allow the user to add a suffix to an input 
# file name and then copy it to a generic file name, for example genray.in_myFile ->
# genray.in.  However Bob Harvey has this functionality in the prepare_genray_input.f90
# code.  He has "genraynml" as command line arg to prepare_genray_input, which then
# writes to the generic name genray.in.  So to avoid confusion I have removed the coding
# providing the suffix processing as used in other components.

# Also found a bug in the exception handling coding for file copies, which I fixed 
# (I think).

# 6/16/2014 DBB
# About zero RF power:  We don't want to run an RF code when the RF power is zero.
# The STEP function needs to produce a partial plasma state only containing the RF data
# (i.e. ps_write_update_file not ps_store_plasma_state).  This really is done with
# plasma state fortran code not with the python netcdf interface.
# So I wrote a simple code called zero_RF_power to set
# all RF source profiles in plasma state to zero and then write a partial plasma state.
# For now this code lives in the GENRAY component directory and it gets built and
# installed by the Makefile there.  The zero_RF_EC_power.f90 code is essentially the same
# as the zero_RF_IC_power.f90, and the modifications to this python component are the
# same as were needed for the rf_ic_toric_abr_mcmd.py component.


# Programmable version: Batchelor (6/2/2013)
# This version allows the aiming angles and power for multiple launchers to be set and
# changed over time from within the simulation config file.  If this option is used, any 
# programming of the launcher power from the EPA component is over-ridden.  This behavior 
# is triggered by the presence a variable N_LAUNCHERS_PROGRAMMED in the [rf_genray] section 
# of the simulation config file.  If N_LAUNCHERS_PROGRAMMED = 0, or is absent, the launcher
# powers are taken from the current plasma state file.  
# If > 0 the component looks for space delimited lists of time points for parameter changes 
# and the associated parameters.  The config lists should be named: 
# LAUNCHER1_TIMES,LAUNCHER1_alfast, LAUNCHER1_betast, LAUNCHER1_powtot_MW, 
# LAUNCHER2_TIMES,LAUNCHER2_alfast, LAUNCHER2_betast, LAUNCHER2_powtot_MW,  etc
# 
# The parameter changes take effect the first time the simulation time at
# the beginning of a time step (ps%t0) is >= LAUNCHERx_TIME then LAUNCHERx parameters are
# changed.  LAUNCHERx_TIMES don't have to match simulation time step beginnings, but the 
# changes won't take place until the beginning of the next time step.
# 
# Note: For ease of typing, enter powers in config file in MW, they get converted to Watts
# in this component.  
# 
# The mechanism for changing the power is to modify the genray.in just before launching
# the genray_EC step.  N_LAUNCHERS_PROGRAMMED must match the number of launchers found in
# the input genray.in file.

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
from ipsframework import Component
#from Numeric import *
from netCDF4 import *
from simple_file_editing_functions import get_lines, put_lines, edit_nml_file

class genray(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

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
                print(message)
                self.services.error(message)
                raise Exception(message)
            print(name, ' = ', str_time_points)

            # get alfasts for this launcher
            name = launcher_name + '_alfast'
            expr = 'self.' + name + '.split()'
            try:
                str_alfast = eval(expr)
            except:
                message = 'error in getting EC_launcher config parameters ' + name
                print(message)
                self.services.error(message)
                raise Exception(message)
            print(name, ' = ', str_alfast)
 
            # get betasts for this launcher
            name = launcher_name + '_betast'
            expr = 'self.' + name + '.split()'
            try:
                str_betast = eval(expr)
            except:
                message = 'error in getting EC_launcher config parameters ' + name
                print(message)
                self.services.error(message)
                raise Exception(message)
            print(name, ' = ', str_betast)

            # get powtots for this launcher
            name = launcher_name + '_powtot_MW'
            expr = 'self.' + name + '.split()'
            try:
                str_powtot_MW = eval(expr)
            except:
                message = 'error in getting EC_launcher config parameters ' + name
                print(message)
                self.services.error(message)
                raise Exception(message)
            print(name, ' = ', str_powtot_MW)
                
            # check that times and parameters have same length
            if len(str_time_points) != len(str_powtot_MW):
                message = 'error: ', launcher_name, \
                ' time points and parameters have different lengths' 
                print(message)
                self.services.error(message)
                raise Exception(message)
                
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
        print('genray.init() called')

        services = self.services
        global programming, times_parameters_list
        
    # Get global configuration parameters
        cur_state_file = self.get_global_param(services, 'CURRENT_STATE')
        cur_eqdsk_file = self.get_global_param(services, 'CURRENT_EQDSK')
        cql_file = self.get_global_param(services, 'CURRENT_CQL', optional = True)

    # Get component-specific configuration parameters. Note: Not all of these are
    # used in 'init' but if any are missing we get an exception now instead of
    # later
        NPROC = self.get_component_param( services, 'NPROC')
        BIN_PATH = self.get_component_param( services, 'BIN_PATH')
        INPUT_FILES = self.get_component_param( services, 'INPUT_FILES')
        OUTPUT_FILES = self.get_component_param( services, 'OUTPUT_FILES')
        RESTART_FILES = self.get_component_param( services, 'RESTART_FILES')
        BIN_PATH = self.get_component_param( services, 'BIN_PATH')
        GENRAY_BIN = self.get_component_param( services, 'GENRAY_BIN')
        RFMODE = self.get_component_param( services, 'RFMODE')
        ISOURCE_STRING = self.get_component_param( services, 'ISOURCE_STRING')
        GENRAYNML = self.get_component_param( services, 'GENRAYNML')
        ADJ_READ = self.get_component_param( services, 'ADJ_READ')
        PS_ADD_NML = self.get_component_param( services, 'PS_ADD_NML')

        # Get [rf_genray_EC] programming configuration parameters, if present
        n_launchers = 0
        programming = False
        try:
            n_launchers = int(self.N_LAUNCHERS_PROGRAMMED)
        except:
            print('\nNo launcher programming in config file')

        if n_launchers > 0: 
            programming = True
            print('Programming for ', n_launchers, ' launchers')
            times_parameters_list = self.time_points_and_parameters(n_launchers)
        else:
            print('\nUsing launcher parmeters in plasma state')


    # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception as e:
          print('Error in call to stage_state()' , e)
          services.error('Error in call to stage_state()')
          raise Exception('Error in call to stage_state()')
        
    # Get input files  
        try:
          services.stage_input_files(INPUT_FILES)
        except Exception as e:
          print('Error in call to stage_input_files()' , e)
          services.error('Error in call to stage_input_files()')
          raise Exception('Error in call to stage_input_files()')
            
    # Copy current plasma state file to generic name -> cur_state.cdf
        try:
            shutil.copyfile(cur_state_file, 'cur_state.cdf')
        except IOError as xxx_todo_changeme:
            (errno, strerror) = xxx_todo_changeme.args
            print('Error copying file %s to %s' % (cur_state_file, 'cur_state.cdf'), strerror)
            services.error('Error copying cur_state_file -> cur_state.cdf')
            raise Exception('Error copying cur_state_file -> cur_state.cdf')

    # Copy current eqdsk file to generic name -> eqdsk
        try:
            shutil.copyfile(cur_eqdsk_file, 'eqdsk')
        except IOError as xxx_todo_changeme1:
            (errno, strerror) = xxx_todo_changeme1.args
            print('Error copying file %s to %s' % (cur_eqdsk_file, 'eqdsk'), strerror)
            services.error('Error copying cur_eqdsk_file -> eqdsk')
            raise Exception('Error copying cur_eqdsk_file -> eqdsk')

        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_genray_input')

        rfmode = self.RFMODE
        isource_string = self.ISOURCE_STRING
        genraynml = self.GENRAYNML
        adj_read = self.ADJ_READ
        ps_add_nml = self.PS_ADD_NML

    # Call prepare_input - init
        print('rf_genray: calling prepare_input init')        
        log_file = open('log_prepare_genray_input_init', 'w')
        mode = 'init'
        command = prepare_input_bin + ' ' + mode + ' ' +  rfmode + ' ' +\
        isource_string + ' ' + genraynml + ' ' + adj_read + ' ' + ps_add_nml
        
        print('running = ', command)
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
          event_comment =  command)

        retcode = subprocess.call(command.split(), stdout = log_file,\
                                  stderr = subprocess.STDOUT)
        if (retcode != 0):
            print('Error executing genray init ', prepare_input_bin)
            services.error('Error executing genray init')
            raise Exception('Error executing genray init')

    # If parameters are programmed from config file set parameters in genray.in
        if programming == True:
            # Get t0 from plasma state
            ps = Dataset(cur_state_file, 'r', format = 'NETCDF3_CLASSIC')
            t0 = ps.variables['t0'].getValue()
            ps.close()
            # set parameters those for time = t0
            self.set_genray_EC_parameters("genray.in", t0, times_parameters_list)

    # Copy generic cur_state.cdf -> current plasma state file
        try:
            shutil.copyfile('cur_state.cdf', cur_state_file)
        except IOError as xxx_todo_changeme2:
            (errno, strerror) = xxx_todo_changeme2.args
            print('Error copying file %s to %s' % ('cur_state.cdf', cur_state_file), strerror)
            services.error('Error copying cur_state.cdf -> cur_state_file')
            raise Exception('Error copying cur_state.cdf -> cur_state_file')

    # Update plasma state files in plasma_state work directory
        try:
          services.update_state()
        except Exception as e:
          print('Error in call to update_state()', e)
          services.error('Error in call to update_state()')
          raise Exception('Error in call to update_state()')
     
    # Archive output files
    # N.B.  prepare_genray_input in init mode does not produce a complete set 
    #       of ourput files.  This causes an error in stage_output_files().
    #       To solve this we generate a dummy set of output files here with 
    #       system call 'touch'
        for file in self.OUTPUT_FILES.split():
          print('touching ', file)
          subprocess.call(['touch', file])
    # Now stage them      
        try:
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception as e:
          print('Error in call to stage_output_files()', e)
          services.error('Error in call to stage_output_files()')
          raise Exception('Error in call to stage_output_files()')

        return 0
# ------------------------------------------------------------------------------
#
# RESTART function
# Gets restart files from restart directory
# Loads the global configuration parameters from the config file
#
# ------------------------------------------------------------------------------
        
    def restart(self, timeStamp):
      print('genray.restart() called')
      print('rf_genray_EC_p')

      if (self.services == None) :
         print('Error in genray: step (): No services')
         raise Exception('Error in genray: step (): No services')
      services = self.services
      global programming, times_parameters_list

    # Get restart files listed in config file.        
      try:
            restart_root = services.get_config_param( 'RESTART_ROOT')
            restart_time = services.get_config_param( 'RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
      except Exception as e:
            print('Error in call to get_restart_files()' , e)
            services.error('genray: error in call to get_restart_files()')
            raise Exception('genray: error in call to get_restart_files()')

    # Get global configuration parameters
      try:
            self.plasma_state_file = services.get_config_param( 'CURRENT_STATE')
            self.eqdsk_file = services.get_config_param( 'CURRENT_EQDSK')
      except:
            print('genray restart: error in getting config parameters')
            services.error('genray restart: error in getting config parameters')
            raise Exception('genray restart: error in getting config parameters')

        # Get [rf_genray_EC] programming configuration parameters, if present
      n_launchers = 0
      programming = False
      try:
          n_launchers = int(self.N_LAUNCHERS_PROGRAMMED)
      except:
          print('\nCould not get launcher programming from config file')

      if n_launchers > 0: 
          programming = True
          times_parameters_list = self.time_points_and_parameters(n_launchers)
      else:
          print('\nUsing ECH launcher parmeters in plasma state')


      return 0
        
# ------------------------------------------------------------------------------
#
# STEP function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print('genray.step() called')

        if (self.services == None) :
           print('Error in genray: step (): No services')
           raise Exception('Error in genray: step (): No services')
        services = self.services
        global programming, times_parameters_list

        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_genray_input')
        process_output_bin  = os.path.join(self.BIN_PATH, 'process_genray_output')
        #genray_bin = os.path.join(self.BIN_PATH, 'genray')

    # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception as e:
          print('Error in call to stage_state()' , e)
          services.error('Error in call to stage_state()')
          raise Exception('Error in call to stage_state()')

        cur_state_file = services.get_config_param( 'CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param( 'CURRENT_EQDSK')
        
    # Copy current plasma state file to generic name -> cur_state.cdf
        try:
            shutil.copyfile(cur_state_file, 'cur_state.cdf')
        except IOError as xxx_todo_changeme3:
            (errno, strerror) = xxx_todo_changeme3.args
            print('Error copying file %s to %s' % (cur_state_file, 'cur_state.cdf'), strerror)
            raise

    # Copy current eqdsk file to generic name -> eqdsk
        try:
            shutil.copyfile(cur_eqdsk_file, 'eqdsk')
        except IOError as xxx_todo_changeme4:
            (errno, strerror) = xxx_todo_changeme4.args
            print('Error copying file %s to %s' % (cur_eqdsk_file, 'eqdsk'), strerror)
            services.error('Error copying cur_eqdsk_file -> eqdsk')
            raise Exception('Error copying cur_eqdsk_file -> eqdsk')

        rfmode = self.RFMODE
        isource_string = self.ISOURCE_STRING
        genraynml = self.GENRAYNML
        adj_read = self.ADJ_READ
        ps_add_nml = self.PS_ADD_NML

# Check if RF power is zero (or effectively zero).  If true don't run GENRAY just
# run zero_RF_power fortran code
        print('cur_state_file = ', cur_state_file)
        ps = Dataset(cur_state_file, 'r', format = 'NETCDF3_CLASSIC')
        if rfmode == 'EC':
            power_ec = ps.variables['power_ec'][:]
            ps.close()
            print('Total EC power = ', sum(power_ec))
            if(sum(power_ec) < 0.001):
                zero_RF_EC_power = get_component_param(self, services, 'ZERO_EC_POWER_BIN')
                retcode = subprocess.call([zero_RF_EC_power, cur_state_file])
                if (retcode != 0):
                    print('Error executing zero_RF_EC_power ')
                    self.services.error('Error executing zero_RF_EC_power')
                    raise Exception('Error executing zero_RF_EC_power')
        elif rfmode == 'LH':
            power_lh = ps.variables['power_lh'][:]
            ps.close()
            print('Total LH power = ', sum(power_lh))
            if(sum(power_ec) < 0.001):
                zero_RF_LH_power = get_component_param(self, services, 'ZERO_LH_POWER_BIN')
                retcode = subprocess.call([zero_RF_LH_power, cur_state_file])
                if (retcode != 0):
                    print('Error executing zero_RF_LH_power ')
                    self.services.error('Error executing zero_RF_LH_power')
                    raise Exception('Error executing zero_RF_LH_power')                 
            else:
                message = 'rf_genray.py: Unimplemented rfmode = ' + RFMODE
                print(message)
                services.exception(message)
                raise


            # N.B. zero_RF_power does not produce a complete set of GENRAY output
            #      files.  This causes an error in stage_output_files().  To
            #      solve this we generate a dummy set of output files here with
            #      system call 'touch'
            for file in self.OUTPUT_FILES.split():
                subprocess.call(['touch', file])

    # Or actually run GENRAY
        else:
            # Call prepare_input - step
            print('rf_genray step: calling prepare_input')
        
            log_file = open('log_prepare_genray_input_step', 'w')
            mode = 'step'
            command = prepare_input_bin + ' ' + mode + ' ' +  rfmode + ' ' +\
              isource_string + ' ' + genraynml + ' ' + adj_read + ' ' + ps_add_nml

            print('running = ', command)
            services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                event_comment =  command)
        
            retcode = subprocess.call(command.split(), stdout = log_file,\
                                      stderr = subprocess.STDOUT)
            if (retcode != 0):
                print('Error executing genray: ', prepare_input_bin)
                services.error('Error executing genray prepare_input')
                raise Exception('Error executing genray prepare_input')

        # If parameters are programmed from config file set parameters in genray.in
            if programming == True:
                # Get t0 from plasma state
                ps = Dataset(cur_state_file, 'r', format = 'NETCDF3_CLASSIC')
                t0 = ps.variables['t0'].getValue()
                ps.close()
                # set parameters those for time = t0
                self.set_genray_EC_parameters("genray.in", t0, times_parameters_list)
        
    #     Launch genray - N.B: Path to executable is in config parameter GENRAY_BIN
            print('rf_genray: launching genray')
            cwd = services.get_working_dir()
            task_id = services.launch_task(self.NPROC, cwd, self.GENRAY_BIN, logfile='log.genray')
            retcode = services.wait_task(task_id)
            if (retcode != 0):
                print('Error executing command: ', self.GENRAY_BIN)
                services.error('Error executing genray')
                raise Exception('Error executing genray')

        # Call process_output - step
            print('rf_genray step: calling process_output')
        
            log_file = open('log_process_genray_output', 'w')
            mode = 'step'
            command = process_output_bin + ' ' +  rfmode + ' ' + isource_string

            print('running', command)
            services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                event_comment =  command)
        
            retcode = subprocess.call(command.split(), stdout = log_file,\
                                      stderr = subprocess.STDOUT)                                  
            if (retcode != 0):
                print('Error executing genray init ', process_output_bin)
                services.error('Error executing genray process_output')
                raise Exception('Error executing genray process_output')
        
    # Copy generic genray partial plasma state file -> EC_cur_state_file
        try:
            partial_file = cwd + '/RF_EC_' + cur_state_file
            shutil.copyfile('RF_GENRAY_PARTIAL_STATE', partial_file )
        except IOError as xxx_todo_changeme5:
            (errno, strerror) = xxx_todo_changeme5.args
            print('Error copying file %s to %s' % ('cur_state.cdf', cur_state_file), strerror)
            raise


# Merge partial plasma state containing updated EC data
        try:
            services.merge_current_state(partial_file, logfile='log.update_state')
            print('merged GENRAY plasma state data ', partial_file)
        except Exception as e:
            print('Error in call to merge_current_state(' , partial_file, ')')
            self.services.error('Error in call to merge_current_state')
            raise Exception('Error in call to merge_current_state')
            
    # Archive output files
        try:
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception as e:
          print('Error in call to stage_output_files()', e)
          services.error('Error in call to stage_output_files()')
          raise Exception('Error in call to stage_output_files()')
        
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
      print('genray.checkpoint() called')
      if (self.services == None) :
         print('Error in genray: checkpoint(): No services')
         services.error('Error in genray: checkpoint(): No services')
         raise Exception('Error in genray: checkpoint(): No services')
      services = self.services
      services.save_restart_files(timestamp, self.RESTART_FILES)

      return 0
        
# ------------------------------------------------------------------------------
#
# finalize function
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print('genray.finalize() called')
# ------------------------------------------------------------------------------
#
# "Private"  methods
#
# ------------------------------------------------------------------------------


    # Try to get config parameter - wraps the exception handling for get_config_parameter()
    def get_global_param(self, services, param_name, optional=False):

        try:
            value = services.get_config_param(param_name)
            print(param_name, ' = ', value)
        except Exception:
            if optional: 
                print('optional config parameter ', param_name, ' not found')
                value = None
            else:
                message = 'required config parameter ', param_name, ' not found'
                print(message)
                services.exception(message)
                raise
        
        return value

    # Try to get component specific config parameter - wraps the exception handling
    def get_component_param(self, services, param_name, optional=False):

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

