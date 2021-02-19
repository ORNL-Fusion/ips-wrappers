#! /usr/bin/env python

"""
gyro_monitor_comp.py 3-4-2015

This is a monitor component to track progress of a GRYO run within IPS.  It is completely 
separate from the IPS monitor_comp.py component in that it does not use Plasma State and 
and it runs concurrently with the gk_gyro.py component.

This component is also different from IPS monitor component in that the monitor file cannont
be defined until after the first time step.  This is because the init() of the gk_gyro.py
component does not produce output files.  So we don't know some of the dimensions until
after the first time step, e.g. the number of columns in the out.gyro.error file.  So this
init function just pushes out the portal RUNID and config file to the www directory.  The
actually initialization of the monitor file happens in step() after the first time step.

This script exports 4 SWIM component functions:
init - defines the netcdf monitor file for monitor data then calls step to insert the
       initial data (t = t0).  Also saves some internal data needed to restart the
       monitor component (monitorVars and ps_VarsList).  These are are pickled to file
       'monitor_restart'
restart - unpickles the internal monitor data previously saved in file 'monitor_restart'
          and loads these global variables.
step - writes monitor data for current time step to the monitor netcdf file and saves it
finalize - sets the global netcdf attribute to tell ELViz to stop watching this
           monitor file

change log:
"""


import sys
import os
import subprocess
import shutil
import pickle

from ipsframework import Component

# Import the necessary Numeric and netCDF modules
from netCDF4 import *
from numpy import *


# ------------------------------------------------------------------------------
#
# Global definitions
#
# ------------------------------------------------------------------------------


debug = False
first_step = True
monitor_fileName = 'monitor_file.nc'



# List of requested variables to monitor (if dependencies are satisfied)
# The list just below is the default containing everything.  In the component it can
# be overwritten with the configuration file data (some day maybe).

requestedVars = []


# List of variables to monitor (which dependencies are satisfied)
monitorVars = []

# List of grids needed to plot monitor variables
monitorGrids = []

# List of files in work directory containing variables to be monitored
filesList = ['out.gyro.t', 'out.gyro.error']

# List of non-grid dimensions needed for other variables - e.g. scalar lists
monitorDims = []

# Dictionary of Plasma State dependencies for each variable to be monitored:

monitorDefinition = {}

gyro_work_path = os.path.join('../','gk_gyro_gk_gyro_2')

# check that gyro_work_path exists

if debug:
    monitorDefinition.update( {    # testing
            'dummy':['S', 'arb', ['stuff']], 'q':['P', 'arb',['ns'],['rho'] ]
    } )

print 'monitor_comp_version = ', monitor_comp_version
print 'metaData = ',monitorDefinition['monitor_comp_metaData']

#----------------------------------------------------------------------------------------------
#
# Define some utility functions
#
#----------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# Open an input file and return the lines
def get_lines(filename):
	try:
	  file = open(filename, 'r')
	except Exception, e:
	  message = 'get_lines: could not open file ' + filename
	  print  message, e
	  raise Exception, message

    lines = file.readlines()
    file.close()
    return lines

class monitor(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)


    

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

    def init(self, timeStamp=0):
        
        print 'gyro_monitor_comp.init() called'

        services = self.services

        workdir = services.get_working_dir()
        run_id = services.get_config_param('PORTAL_RUNID')
        monitor_file = 'monitor_file.nc'
    #      print 'monitor file = ', monitor_file

        self.cdfFile = run_id+'_monitor_file.nc'
        services.log('w3 monitor file = ' + self.cdfFile)
        htmlFile = run_id +'.html'

    # Get input files  
        try:
          work_dir_path = self.WORK_DIR_PATH
        except Exception, e:
          message = 'gyro_monitor_comp.init: failed to get WORK_DIR_PATH'
          print  message, e
          services.error(message)
          raise Exception, message

    # Generate initial monitor file
        retcode = self.init_monitor_file(cur_state_file, timeStamp)
        if (retcode != 0):
            services.log('Error executing INIT:  init_monitor_file')
            return 1

    # copy monitor file to w3 directory
        try:
            shutil.copyfile(monitor_file,
                            os.path.join(self.W3_DIR, self.cdfFile))
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % \
                (monitor_file, self.cdfFile, strerror)

        htmlText = self.htmlText.replace('@CDF_FILE@',
                           os.path.join(self.W3_BASEURL, self.cdfFile))
        try:
            f = open(os.path.join(self.W3_DIR, htmlFile), 'w')
            f.write(htmlText)
            f.close()
        except IOError, (errno, strerror):
            print 'Error writing to file %s : %s' % \
                (htmlFile, strerror)
            return
        monitorURL = os.path.join(self.W3_BASEURL , htmlFile)
        services.setMonitorURL(monitorURL)

    # Copy config file to w3 directory
        conf_file = services.get_config_param('SIMULATION_CONFIG_FILE')
        print 'conf_file = ', conf_file
        conf_file_name = os.path.split(conf_file)[1]
        new_file_name = run_id + '_' + conf_file_name
        new_full_path = os.path.join(self.W3_DIR, new_file_name)
        try:
            shutil.copyfile(conf_file, new_full_path)
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % \
                (conf_file, new_full_path, strerror)

        return

# ------------------------------------------------------------------------------
#
# restart function
#
# Load the internal state needed to restart the monitor component.  In particular
# monitorVars and ps_VarsList are unpickeld pickled from a file "monitor_restart".
#
# ------------------------------------------------------------------------------

    def restart(self, timeStamp):
        """
        Function restart loads the internal monitor state data needed for
        restart_MonitorComponent
        """
        print 'monitor_comp.restart() called'

        services = self.services
        global monitorVars, ps_VarsList, monitorDefinition
        
        workdir = services.get_working_dir()
        run_id = services.get_config_param('PORTAL_RUNID')
        monitor_file = 'monitor_file.nc'
    #      print 'monitor file = ', monitor_file

        self.cdfFile = run_id+'_monitor_file.nc'
        services.log('w3 monitor file = ' + self.cdfFile)
        htmlFile = run_id +'.html'
        
    # Get restart files listed in config file.        
        try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
        except Exception, e:
            print 'Error in call to get_restart_files()' , e
            raise

        # copy monitor file to w3 directory
        try:
            shutil.copyfile(monitor_file,
                            os.path.join(self.W3_DIR, self.cdfFile))
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % \
                (monitor_file, self.cdfFile, strerror)

        htmlText = self.htmlText.replace('@CDF_FILE@',
                           os.path.join(self.W3_BASEURL, self.cdfFile))
        try:
            f = open(os.path.join(self.W3_DIR, htmlFile), 'w')
            f.write(htmlText)
            f.close()
        except IOError, (errno, strerror):
            print 'Error writing to file %s : %s' % \
                (htmlFile, strerror)
        monitorURL = os.path.join(self.W3_BASEURL , htmlFile)
        self.services.setMonitorURL(monitorURL)
    
        # Load monitorVars and ps_VarsList from pickle file "monitor_restart".

        pickleDict = {'monitorVars' : monitorVars, 'ps_VarsList': ps_VarsList,\
                     'monitorDefinition':monitorDefinition}
#        pickleDict = {'monitorVars' : monitorVars, 'ps_VarsList': ps_VarsList}
        pickFile = open('monitor_restart', 'r')
        pickleDict = pickle.load(pickFile)
        pickFile.close()
        monitorVars = pickleDict['monitorVars']
        ps_VarsList = pickleDict['ps_VarsList']
        monitorDefinition = pickleDict['monitorDefinition']
        print 'monitorDefinition = ', monitorDefinition
        
        print 'monitor restart finished'
        return 0

# ------------------------------------------------------------------------------
#
# step function
#
# Stages current input files (primarily plasma state). Updates the monitor_file.nc
# from current plasma state.  And copies updated monitor file to w3 directory
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print '\nmonitor_comp.step() called'
        global first_step, monitor_fileName

        services = self.services
        if (self.services == None) :
            print 'Error in monitor_comp: step() : no framework services'
            return 1

	# If this is the first call to step() initialize monitor file
		if first_step == True:
		
			# Check that gyro gyro has run long enough to produce output files, i.e. 
			# they exist.  If not wait until they appear.
			
			put code here

			first_step == False
			self.init_monitor_file(cur_state_file, timeStamp = 0)		
		
        monitor_file = monitor_fileName

    # Call Load new data into monitor file
        retcode = self.update_monitor_file(cur_state_file, timeStamp)
        if (retcode != 0):
            print 'Error executing command: update_monitor_file'
            return 1

    # "Archive" output files in history directory
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    # copy montor file to w3 directory
        try:
            shutil.copyfile(monitor_file,
                            os.path.join(self.W3_DIR, self.cdfFile))
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % \
                (monitor_file, self.W3_DIR, strerror)
        return

# ------------------------------------------------------------------------------
#
# checkpoint function
#
# Saves restart files to restart directory.  Should include: monitor_restart and
# monitor_file.nc
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print 'monitor.checkpoint() called'
        services = self.services
        services.save_restart_files(timestamp, self.RESTART_FILES)
        

# ------------------------------------------------------------------------------
#
# finalize function
#
# Calls monitor executable in "FINALIZE" mode which sets
# 'running' atribute to false so Elvis will stop watching the monitor file.
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print 'monitor finalize finished'
        return 0


#----------------------------------------------------------------------------------------------
#
# init_monitor_file function
#
# Analyze the Plasma State dependencies of the variables requested to be monitored.
# Define the initial monitor netcdf file.  Load the initial data into the monitor file
# By a call to the 'step' function.  Also saves the internal state needed to restart
# the monitor component.  In particular monitorVars and ps_VarsList are pickled to a
# file called monitor_restart.
#----------------------------------------------------------------------------------------------

    def init_monitor_file(self, cur_state_file, timeStamp = 0):
        """Function init_monitor_file generates the initial netcdf file for monitor data
           and saves the internal state needed to restart the monitor component
        """

        print ' '
        print 'gyro_monitor_component: init_monitor_file'
        
    # Get global configuration parameters
        try:
            Global_label = services.get_config_param('SIM_NAME')
		except Exception, e:
		  message = 'gyro_monitor_comp init(): could not get config parameter SIM_NAME'
		  print  message, e
		  raise Exception, message
        print 'Global_label = ', Global_label

        try:
            RUN_ID = services.get_config_param('SIM_NAME')
		except Exception, e:
		  message = 'gyro_monitor_comp init(): could not get config parameter RUN_ID'
		  print  message, e
		  raise Exception, message
        print 'RUN_ID = ', RUN_ID

        try:
            tokamak_id = services.get_config_param('SIM_NAME')
		except Exception, e:
		  message = 'gyro_monitor_comp init(): could not get config parameter tokamak_id'
		  print  message, e
		  raise Exception, message
        print 'tokamak_id = ', tokamak_id

        try:
            shot_number = services.get_config_param('SIM_NAME')
		except Exception, e:
		  message = 'gyro_monitor_comp init(): could not get config parameter shot_number'
		  print  message, e
		  raise Exception, message
        print 'shot_number = ', shot_number

	# Get monitor variables
	
		# Get data from out.gyro.error
		input_lines = getlines(os.path.join(gyro_work_path, 'out.gyro.error'))
		# Find out how many species are included i.e. number of columns in file
		n_species = len(input_lines[0].split())
		for i in range(n_species):
			var_name = 'species_' + str(i+1)
			monitorVars.append(var_name)
			monitorDefinition.update( {var_name: ['S', ' ', [] ]})
	
	#
	## Define monitor file
	#

        # Open the monitor file for output
        monitor_file = Dataset(monitor_fileName, 'w', format = 'NETCDF3_CLASSIC')

        # make global netcdf attribute of monitor component version number
        setattr(monitor_file, 'monitor_comp_version', monitor_comp_version)

        # make global netcdf attribute for Global_label
        setattr(monitor_file, 'Global_label', Global_label)

        # make global netcdf attribute for RunID
        setattr(monitor_file, 'RunID', RunID)

        # make global netcdf attribute for tokamak_id
        setattr(monitor_file, 'tokamak_id', tokamak_id)

        # make global netcdf attribute for shot_number
        setattr(monitor_file, 'shot_number', shot_number)

        # Create unlimited time dimension and define time variable
        monitor_file.createDimension('timeDim', None)
        time = monitor_file.createVariable('time', float, ('timeDim',))
        setattr(time, 'units', '-')

        # Create grid dimensions and variables and load up grids
        # use this to keep track of mon_grid name, grid_dim_name and grid_dim_value
        grid_map = {}
        for grid in monitorGrids:
			pass

        # Create monitor variables in netCDF4
        for var in monitorVars:

            if debug:
                print 'creating variable ', var

            # Generate the dimension tuple
            dims = ('timeDim',)
            varKind = monitorDefinition[var][0]

            if varKind == 'S':
                mon_obj = monitor_file.createVariable(var, float, dims)

            elif varKind == 'P':
                varGridsList = monitorDefinition[var][3]
                # add in dimension names for plotting grids
                for grid in varGridsList:
                    dims = dims + grid_map[grid][1]

                mon_obj = monitor_file.createVariable(var, float, dims )

            elif varKind == 'SL':
                varDimList = monitorDefinition[var][3]
                # add in dimension names
                for dim_name in varDimList:
                    dims = dims + (dim_name,)

                mon_obj = monitor_file.createVariable(var, float, dims )

            elif varKind == '2D':   # See note at top about the order of dims in definition
                varGridsList = monitorDefinition[var][3]
                for grid in varGridsList:  # add in dimension names for plotting grids
                    dims = dims + grid_map[grid][1]

                if debug:
                    print var, ' dimensions = ', dims
                mon_obj = monitor_file.createVariable(var, float, dims )

            varUnits = monitorDefinition[var][1]
            #Generate units atribute for var
            setattr(mon_obj, 'units', varUnits)

        # Finished defining monitor file.  Close it.
        monitor_file.close()

        # insert intitial data
        self.step(timeStamp)

        print 'monitor file initialization finished'
        
    # Save monitorVars and ps_VarsList and monitorDefinition, are pickled to file "monitor_restart".
    
        pickleDict = {'monitorVars' : monitorVars, 'monitorDefinition': monitorDefinition}
        pickFile = open('monitor_restart', 'w')
        pickle.dump(pickleDict, pickFile)
        pickFile.close() 

        return 0

#----------------------------------------------------------------------------------------------
#
# update_monitor_file function
#
# Opens current plasma state file and monitor netcdf file.  Pulls the needed variables out of
# the plasma state file.  Computes the monitor variables from plasma state data (in functinon
# calculate_MonitorVariable).  And loads the data into the monitor netcdf variables.
#
#----------------------------------------------------------------------------------------------

    def update_monitor_file(self, cur_state_file, timeStamp = 0):

        print ' '
        print 'monitor_component: update_monitor_file'

        #Open Plasma State file
        plasma_state = Dataset(cur_state_file, 'r', format = 'NETCDF3_CLASSIC')

        # Get all Plasma State variable objects for dependencies
        ps_variables = {}
        for var in ps_VarsList:
            ps_variables[var] = plasma_state.variables[var]

        if debug:
            print ' '
            print 'ps_variables.keys() = ', ps_variables.keys()

        # Open the monitor file for output
        monitor_file = Dataset(monitor_fileName, 'r+', format = 'NETCDF3_CLASSIC')

        if debug:
            all_mon_Dims = monitor_file.dimensions
            all_mon_VarNames = monitor_file.variables.keys()
            print ' '
            print 'all_mon_Dims = ', all_mon_Dims
            print ' '
            print 'all_mon_VarNames = ', all_mon_VarNames

        # Get time variable object
        time = monitor_file.variables['time']
        n_step =time.shape[0]    # Time step number (initially 0)
        print 'time step number = ', n_step
        time[n_step] = float(timeStamp)

        # Insert data into monitor variables
        for var in monitorVars:
            # Get the netcdf variable object for this variable
            var_obj = monitor_file.variables[var]

            # Get Plasma State variables for this monitor variable's dependencies
            varDepsDict ={}
            varDepsList = monitorDefinition[var][2]
            for dep in varDepsList:
                varDepsDict[dep] = ps_variables[dep]

            if debug:
                print ' '
                print 'var =', var
                print 'varDepsDict.keys = ', varDepsDict.keys()
                print 'varDepsDict = ', varDepsDict
            # calculate the monitor variable
            value = self.calculate_MonitorVariable(var, varDepsDict)

            # Load the data into monitor variables
            varKind = monitorDefinition[var][0]
            if varKind == 'S':
                var_obj[n_step] = value
            if varKind == 'P':
                var_obj[n_step,:] = value
            if varKind == 'SL':
                var_obj[n_step,:] = value
            if varKind == '2D':
                var_obj[n_step,:] = value



        # Finished writing to monitor file.  Close it.
        monitor_file.close()

        # Close plasma_state
        plasma_state.close()

        print 'update_monitor_file finished'
        return 0


#----------------------------------------------------------------------------------------------
#
# calculate_MonitorVariable function
#
# For given montor variable it calculates the variable value from the netcdf variable objects
# in the varDepsDict dictionary on which the monitor variable depends.
#
#----------------------------------------------------------------------------------------------

    def calculate_MonitorVariable(self, var, varDepsDict):

        if debug:
            print 'calculate_MonitorVariable, var = ', var
            for dep_name in varDepsDict.keys():
                print dep_name, '= ', varDepsDict[dep_name][:]
                print dep_name, ' shape = ', varDepsDict[dep_name].shape
        
        if debug:
            print ' '
            print var, ' = ', value
    
        return value
