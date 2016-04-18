#! /usr/bin/env python

# version 0.4 4/27/08 (Batchelor)

# ------------------------------------------------------------------------------
#
# EPA component script to drive model_EPA_mdescr_mdescr.f90 executable.
#
# !      The executable requires3 commandline arguments:
# !      1) current state file!
# !      2) mode = one of "INIT", "STEP", "FINALIZE"
# !      3) timeStamp = initial time for "INIT", or = time at end of time stamp for "STEP"
#
# ------------------------------------------------------------------------------

import sys
import os
import subprocess
import getopt
import shutil
import string
from netCDF4 import *
from  component import Component

parameterList = ['Te_0', 'Te_edge', 'alpha_Te_1', 'alpha_Te_2', 'ne_0', 'ne_edge',\
    'alpha_ne_1', 'alpha_ne_2',\
    'frac_ni', 'frac_Ti', 'fracmin_T', 'fracmin_n', \
    'Te_ratio', 'alpha_Te', 'ne_ratio', 'alpha_ne',\
    'T_min_0', 'T_min_ratio', 'alpha_Tmin',\
    'fracmin', 'power_ic']

class model_EPA_mdescr(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

# ------------------------------------------------------------------------------
#
# init function
#
# model_EPA_mdescr init function allocates plasma profiles and initializes rf sources
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp):
        print 'model_EPA_mdescr.init() called'
        print 'adjustable model parameters = ', parameterList

        services = self.services

# Copy current and prior state over to working directory
        self.services.stage_plasma_state()

        self.services.stage_input_files(self.INPUT_FILES)
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        bin = os.path.join(self.BIN_PATH, 'model_EPA_mdescr')

        print 'Executing ', [bin, cur_state_file, 'INIT', timeStamp]
        retcode = subprocess.call([bin, cur_state_file, 'INIT', timeStamp])
        if (retcode != 0):
            message = 'generic_ps_init: Error executing' + bin
            print message
            services.exception(message)
            raise

# Update (original) plasma state
        services.update_plasma_state()

        return

# ------------------------------------------------------------------------------
#
# STEP function
#
# Time dependence of parameters is implmemented here by modifying namelist file
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'model_EPA_mdescr.step() called'
        global parameterList
        services = self.services


# Copy current and prior state over to working directory
        services.stage_plasma_state()
        cur_state_file = services.get_config_param('CURRENT_STATE')
        ps = Dataset(cur_state_file, 'r+', format = 'NETCDF3_CLASSIC')
        t0 = ps.variables['t0'].getValue()

# Time evolution of parameters
        # get lines from namelist file
        inputLines = self.get_lines('model_EPA_mdescr_input.nml')

        evolution_models = {'linear_DT': self.linear_DT}
        print ' '
        print 'evolution_models = ', evolution_models.keys()
        
        # Look in config file for parameters to evolve, get the evolution model and its  
        # arguments
        params_to_change = False
        for param in parameterList:
            model_name = self.try_get_component_param(services, param + '_DT_model', \
                optional = True)
            if model_name != None:
            	model_name = model_name.strip()
            print 'model_name = ', model_name
            print model_name == 'linear_DT'
            params_to_change = True
            if model_name != None:
                if model_name == 'linear_DT':
                    print 'time evolution model = ', model_name
                    DT_param = self.try_get_component_param(services, param + '_DT_param')
                    print param + '_DT_param = ', DT_param
                    paramValue = read_var_from_nml_lines(self, inputLines, param, separator = ',')
                    print 'value for ', param, ' = ', paramValue
                    newValue = linear_DT(self, float(paramValue), timestamp, t0, float(DT_param))
                    print 'new value for ', param, ' = ', newValue

                    # modify that parameter in namelist file
                    lines = self.edit_nml_file(inputLines, param, newValue, separator = ',')
        
        # write modified namelist file        
        if params_to_change:
            self.put_lines('model_EPA_mdescr_input.nml', lines)

# Call model_EPA_mdescr
        bin = os.path.join(self.BIN_PATH, 'model_EPA_mdescr')
        print 'Executing ', [bin, cur_state_file, 'STEP', timeStamp]
        retcode = subprocess.call([bin, cur_state_file, 'STEP', timeStamp])
        if (retcode != 0):
            message = 'generic_ps_init: Error executing' + bin
            print message
            services.exception(message)
            raise

# Update plasma state
        services.update_plasma_state()

# "Archive" output files in history directory
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

# ------------------------------------------------------------------------------
#
# finalize function
#
# Does nothing
# ------------------------------------------------------------------------------



    def finalize(self, timestamp=0.0):
        print 'model_EPA_mdescr finalize() called'

# ------------------------------------------------------------------------------
#
# "Private"  methods
#
# ------------------------------------------------------------------------------

    # Linear time advance f(timestamp) = f(t0) + (timestamp - t0)*DT
    def linear_DT(self, f, timestamp, t0, DT):
        return f + (timestamp - t0)*DT

    # Try to get config parameter - wraps the exception handling for get_config_parameter()
    def try_get_config_param(self, services, param_name, optional=False):

        try:
            value = services.get_config_param(param_name)
            print param_name, ' = ', value
        except Exception:
            if optional:
                print 'config parameter ', param_name, ' not found'
                value = None
            else:
                message = 'required config parameter ', param_name, ' not found'
                print message
                services.exception(message)
                raise

        return value

    # Try to get component specific config parameter - wraps the exception handling
    def try_get_component_param(self, services, param_name, optional=False):

        if hasattr(self, param_name):
            value = getattr(self, param_name)
            print param_name, ' = ', value
        elif optional:
            print 'optional config parameter ', param_name, ' not found'
            value = None
        else:
            message = 'required component config parameter ', param_name, ' not found'
            print message
            services.exception(message)
            raise

        return value


    # A utility to do simple editing of txt files on a line by line basis.
    #---------------------------------------------------------------------------------------
    # Open an input file and return the lines
    def get_lines(self, filename):
        file = open(filename, 'r')
        lines = file.readlines()
        file.close()
        return lines

    #---------------------------------------------------------------------------------------
    # Open an output file and write lines into it
    def put_lines(self, filename, lines):
        file = open(filename, 'w')
        file.writelines(lines)
        file.close()


    #---------------------------------------------------------------------------------------
    # Editing utilities
    #---------------------------------------------------------------------------------------

    def edit_nml_file(self, lines, var, values, separator = ','):

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
        print 'var_lines = ', var_lines     
            
        # Insert line with new values at lines[i]
        lines[var_line_number] = lines[var_line_number].split('=')[0] + ' = '
        for val in values:
            lines[var_line_number] = lines[var_line_number] + str(val) + ', '
        lines[var_line_number] = lines[var_line_number] + '\n'
        print 'New ', lines[var_line_number]
        
        return lines[:var_line_number + 1] + lines[var_line_number + var_lines:]

    def read_var_from_nml_lines(self, lines, var, separator = ','):
    # This routine is very limited, for now it only reads scalar real numbers

        # Find the line in the namelist containing 'var = ' 
        var_line_number = -1
        for i in range(len(lines)):
            line = lines[i]
            if '=' in line:
                split_line = line.split('=')
                #print 'split_line = ', split_line
                if (split_line[0].strip()).lower() == var.lower():
                    var_line_number = i

        if var_line_number == -1:
            message = 'read_var_from_nml_lines: Could not find variable ', var, ' in namelist lines'
            print message
            raise Exception(message)
    
        RHS = lines[var_line_number].split('=')[1]
        
    #     Get rid of newline if there is one
        if RHS[-1] == '\n':
            lines[var_line_number] = lines[var_line_number][:-1]
    
        RHS_list = RHS.split(',')
    #    print 'RHS = ', RHS_list[0].split()
        value = float(RHS_list[0].split()[0])
        print 'value = ', value
        return value
