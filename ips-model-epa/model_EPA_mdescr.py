#! /usr/bin/env python

"""
model_EPA_mdescr.py version 1.1 3/31/2017 (Batchelor)

EPA component script to drive model_EPA_mdescr_mdescr.f90 executable.  See comment header
in model_EPA_mdescr_mdescr.f90 for details.  

Note: Time evolution for the models is implemented in this script.  The fortran,
      model_EPA_mdescr_mdescr.f90, generates profiles and communicates with plasma state.

The executable, model_EPA_mdescr_mdescr.f90, requires 3 commandline arguments:
1) current state file!
2) mode = one of "INIT", "STEP", "FINALIZE"
3) timeStamp = initial time for "INIT", or = time at end of time stamp for "STEP"

For applications that only require an initial plasma state with no time evolution (e.g.
iteration of TORIC and CQL3D with fixed thermal profiles and equilibrium) I have added
an optional config parameter to the EPA section of the simulation config file, INIT_ONLY.
If INIT_ONLY == true, then when the STEP function is called it just returns.  I did this
so that it can be used this way with the generic drivers which automatically run the
EPA component in the time loop if one is present in the config PORTS.

Details of how time evolution models are specified:
The names of the parameters to be evolved, the name of the time dependance models and the
parameters of the evolution models are specified in the [[EPA]] section of the simulation
configuration file.  The list of changeable parameters is in parameterList below.

Specification of a parameter to be evolved requires a line in the configuration file of 
the form:

<parameter>_DT_model = <model name>    (e.g. Te_0_DT_model = ramp_initial_to_final)

There must also be a line giving the parameters of the evolution model of the form:

<parameter>_DT_params = <space separated list of values of the parameters of the evolution model>

(e.g. Te_0_DT_params = <time to start ramp> <time to end ramp> <final temperature>
N.B. We expect to get the intial value of the parameter from the initial input 
     namelist file model_EPA_mdescr_input.nml)

Eventually I should make these keyword based, but for now you just have to know what the
ordered parameter list is for you evolution model.

The list of evolution models implemented so far is very short, one -> ramp_initial_to_final

"""

# ------------------------------------------------------------------------------
# Working Notes
# 3-31-2017
# Added INIT_ONLY
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
    'fracmin', 'power_ic', 'power_lh']

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
        
# "Archive" output files in history directory
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

# Copy initial namelist file so original parameters will be available for time evolution
        try:
            shutil.copyfile('model_EPA_mdescr_input.nml', 'initial_input.nml')
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % ('machine.inp' + suffix, 'machine.inp', strerror)
            logMsg = 'Error copying machine.inp_<suffix> -> machine.inp'
            services.exception(logMsg)
            raise

        return

# ------------------------------------------------------------------------------
#
# STEP function
#
# Time dependence of parameters is implmemented here by modifying namelist file
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        services = self.services
        init_only = self.try_get_component_param(services, 'INIT_ONLY', optional = True)
        if init_only in ['TRUE', 'True', 'true']: return

        print 'model_EPA_mdescr.step() called'
        global parameterList


# Copy current and prior state over to working directory
        services.stage_plasma_state()
        cur_state_file = services.get_config_param('CURRENT_STATE')
        ps = Dataset(cur_state_file, 'r+', format = 'NETCDF3_CLASSIC')
        tinit = ps.variables['tinit'].getValue()

# Time evolution of parameters
        # get lines from namelist file
        inputLines = self.get_lines('model_EPA_mdescr_input.nml')
        initial_nml_Lines = self.get_lines('initial_input.nml')

        evolution_models = {'linear_DT': self.linear_DT,\
                            'ramp_initial_to_final': self.ramp_initial_to_final}
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
                params_to_change = True
        
                if model_name == 'ramp_initial_to_final':
                    print 'model_EPA_mdescr: ramp_initial_to_final'
                    DT_paramsList = self.try_get_component_param(services, param + '_DT_params').split()
                    t_initial = float(DT_paramsList[0])
                    t_final = float(DT_paramsList[1])
                    
                    # Get initial value of parameter from the initial namelist file
                    Value_init = self.read_var_from_nml_lines(initial_nml_Lines, param, separator = ',')
                    print 'intial '+param, ' = ', Value_init
                    
                    #Value_init = float(DT_paramsList[2])
                    Value_final = float(DT_paramsList[2])
                    print 't_initial = ',t_initial, ' t_final = ', t_final,\
                    '  Value_init =  ', Value_init, '  Value_final =  ', Value_final
                    newValue = self.ramp_initial_to_final(float(timeStamp), t_initial,\
                               t_final, Value_init, Value_final)
                    print 't = ', float(timeStamp), ' ', param, ' = ', newValue

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

    # f = f0 for t < t0, f= f1 for t >  t1, linear ramp in between
    def ramp_initial_to_final(self, t, t0, t1, f0, f1):
        if (t1 - t0) < 0.:
            message = 'invalid beginning/end times for ramp  t0 = ', t0, ' t1 = ', t1
            print message
            services.exception(message)
            raise
        
        if t <= t0: return f0
        if t <= t1: return f0 + (f1 - f0)*(t - t0)/(t1 - t0)
        if t > t1:  return f1
        
    # Linear time advance f(timestamp) = f(t0) + (timestamp - t0)*DT
    def linear_DT(self, f0, t, t0, DT):
        return f0 + (t - t0)*DT

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
