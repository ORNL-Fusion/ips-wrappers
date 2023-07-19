#! /usr/bin/env python

"""
simple_file_editing_functions.py 3/13/2018 (Batchelor)
Utilities to read, modify, and write text files with lines in the form of assignment
statements, i.e. of the form <name> = <value>
For now it only deals with single line assignments.  The line can have multiple values
(i.e. a vector) but the line has to be parsed to separate the multiple values.
"""

# Working notes:
# Batchelor 4/21/2020: Copied futurized component from dbb4 branch
# Batchelor 2/1/2021: Added read_var_from_nml_lines
# Batchelor 8/12/2022: Added read_var_from_fortran_star_file
# Batchelor 11/27/2022: Extended read_var_from_fortran_star_file to allow var = value
#                       on a single line

import sys

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
    print('Got to put_lines')
    file = open(filename, 'w')
    file.writelines(lines)
    print('Did file.writelines')
    file.close()

#---------------------------------------------------------------------------------------
# Editing utilities
#---------------------------------------------------------------------------------------

def lines_to_variable_dict(lines):
# Parses lines of the form "variable_name = variable value" into dictionary {name:value}
    variable_dict = {}
    for line in lines:
        if (line.strip() not in ['!', '#']) and (len(line) > 0):  # ignore comments
            if '=' in line:
                name, val = line.split("=")
                name = name.strip()
                val = val.strip()
                variable_dict[name] = val
    return variable_dict

def lines_to_list(lines):
# Collects lines into list, one line per list entry
    line_list = []
    for line in lines:
# Get rid of newline if there is one
        if line[-1] == '\n':
            line = line[:-1]
        line_list.append(line)
    return line_list

def list_variables_in_fortran_star_file(lines):
# Assumes lines come in pairs: First a line with variable name, second a line with variable 
# value.

    variable_list = []
    npairs = int(len(lines)/2)
    
    # Check if an even number of lines
    if (2*npairs != len(lines)):
    	err_mess = 'list_variables_in_fortran_star_file: odd number of lines = ' + str(len(lines))
    	print(err_mess)
    	raise Exception(err_mess)
    	
    for i in range(npairs):
        name = lines[2*i].strip()
        variable_list.append(name)
    return variable_list

def variable_dict_to_lines(variable_dict):
    lines = []
    for var in list(variable_dict.keys()):
        lines.append(var + ' = ' + str(variable_dict[var]) + '\n')
    return lines


def input_file_to_variable_dict(filename):   
    lines = get_lines(filename)
    variable_dict = lines_to_variable_dict(lines)
    return variable_dict


def variable_dict_to_output_file(variable_dict, filename):
    lines = variable_dict_to_lines(variable_dict)
    put_lines(filename, lines)

def add_variables_to_output_file(variable_dict, filename):
    lines = get_lines(filename)
    for name in list(variable_dict.keys()):
        new_line = name + ' = ' + str(variable_dict[name]) + '\n'
        lines.append(new_line)
    put_lines(filename, lines)

            
def modify_variables_in_file(change_dict, filename):
    lines = get_lines(filename)

    # Find the line in the file containing 'var = ' and change value
    var_line_number = -1
    for i in range(len(lines)):
        line = lines[i]
 
        # ignore blank lines
        if len(line.strip()) == 0 : continue
        # ignore comments, could be fooled by '=' in comment
        if (line.strip()[0] in ['!', '#']): continue

        if '=' in line:
            name, val = line.split("=")
            name = name.strip()
            if name in list(change_dict.keys()):
                lines[i] = name + ' = ' + str(change_dict.pop(name)) + '\n'
                
    if len(change_dict) != 0:
        for name in list(change_dict.keys()):
            print('variable to be modified ', name, ' not found in file ', filename)
            raise
    put_lines(filename, lines)

#---------------------------------------------------------------------------------------
# Convert a dictionary variable containing whitespace delimited multiple variables to
# a list of floats.  Useful when lines_to_variable_dict() hits a vector on a single line.
#---------------------------------------------------------------------------------------

def dict_variable_to_list_of_floats(variable_dict, variable):
    chr_list =  variable_dict[variable].split()
    num_list = [float(x) for x in chr_list]
    return num_list

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
        print(message)
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
        next_line_no = var_line_number + var_lines
        next_line = lines[next_line_no]
        if not next_line.strip() : # Check for blank line and skip
            var_lines += 1
            #print('var_line_number = ', var_line_number, 'next_line = ', next_line)
            continue
        if '=' in next_line:   # Could get fooled by = in a quoted string
            test = True
            eq_index = next_line.find('=') # so check if quote before =
            single_quote_index = next_line.find("'")
            if single_quote_index > -1 and single_quote_index < eq_index:
                test = False
            double_quote_index = next_line.find('"')
            if double_quote_index > -1 and double_quote_index < eq_index:
                test = False
        elif next_line[-2] == '/':  # At end of line means end of group. [-1] is newline
            test = True
        else:
            var_lines += 1
            
    # Insert line with new values at lines[i]
    lines[var_line_number] = lines[var_line_number].split('=')[0] + ' = '
    for val in values:
        lines[var_line_number] = lines[var_line_number] + str(val) + ', '
    lines[var_line_number] = lines[var_line_number] + '\n'
    print('New ', lines[var_line_number])
        
    return lines[:var_line_number + 1] + lines[var_line_number + var_lines:]

#---------------------------------------------------------------------------------------
# Read a variable from a fortran namelist file
#---------------------------------------------------------------------------------------

def read_var_from_nml_lines(lines, var, separator = ','):

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
        print(message)
        raise Exception(message)

    RHS = lines[var_line_number].split('=')[1]

#     Get rid of newline if there is one
    if RHS[-1] == '\n':
        lines[var_line_number] = lines[var_line_number][:-1]

    RHS_list = RHS.split(',')
#    print 'RHS = ', RHS_list[0].split()
    value = float(RHS_list[0].split()[0])
    print('value = ', value)
    return value


#---------------------------------------------------------------------------------------
# Find and list the variables defined an nml file 
#---------------------------------------------------------------------------------------

def list_variables_in_nml_file(lines):

    var_list = []
    for i in range(len(lines)):
        line = lines[i]
        if '=' in line:
            split_line = line.split('=')
            var_list.append(split_line[0])
    return var_list

#---------------------------------------------------------------------------------------
# Read a variable from a file of simple format. 
# For example a fortran file written with list directed format (unit,*)
# The data must be either in the form: <variable name> = <data record/data values>,
# or pairs of lines: one line with variable name: <variable name> or <variable name> = 
# followed by one line with data values.  
# If there is only one value it returns a scalar, otherwise it returns a list.
#
# Returns values of var_type which must be in ['string', 'int', 'float']. The user
# must know the names and type of variables to be read and must specify the type in the
# call.
#---------------------------------------------------------------------------------------

def read_var_from_fortran_star_file(lines, var, var_type, separator = ''):

    if var_type.strip() not in ['string', 'int', 'float']:
        message = 'unrecognized variable type ',type, ' for ', var.strip()
        print(message)
        raise Exception(message)
    
    # Find the line containing 'var' also 'var = ' is allowed
    var_line_number = -1
    for i in range(len(lines)):
#        line = lines[i].strip()
        line = lines[i].strip('\n').strip() # get rid of new_line
        
        # If just var with no '=' this better be the first of a pair of lines
        if '=' not in line:
            if line.split()[0] == var.strip():
                var_line_number = i
                if separator == '' :
                    value = lines[i+1].split()
                else:
                    value = lines[i+1].split(separator)

        else:
        # If there is an '=' this could be the first of a pair or a var/value line
            split_line = line.split('=')
            LHS = split_line[0].strip()
            RHS = split_line[1]            
            if LHS == var.strip():
                var_line_number = i
                if len(RHS) == 0: # Nothing after the '=' means first of a line pair
                    if separator == '' :
                        value = lines[i+1].split()
                    else:
                        value = lines[i+1].split(separator)
                else:            
                     if separator == '' :
                        value = RHS.split()
                     else:
                        value = RHS.split(separator)

    if var_line_number == -1:
        message = 'read_var_from_fortran_star_file: Could not find variable ',\
                   var, ' in file'
        print(message)
        raise Exception(message)
    
#    If var_type == 'string' do nothing. Otherwise cast as int or float.

    if var_type == 'int' :
        for i in range(len(value)) :
            value[i] = int(value[i])

    if var_type == 'float' :
        for i in range(len(value)) :
            value[i] = float(value[i])          
      
    
    if len(value) == 1: return value[0]
    else: return value

#_________________________________________________________________________________________________

if __name__ == '__main__':

    input_file = 'genray.in'
    lines = get_lines(input_file)        
    var = 'powtot'

#     val = read_var_from_nml_lines(lines, var, separator = ',')    
#     print(var, ' = ', val)

    powtot_t=['1.0d0','2.0d0','3.0d0','4.0d0','5.0d0']

    lines = edit_nml_file(lines, 'powtot', powtot_t, separator=',')
    
    put_lines("new_file", lines)
#    new_lines = get_lines("new_file")
#     val = read_var_from_nml_lines(new_lines, var, separator = ',')    
#     print(var, ' = ', val)

#     input_file = 'scan_summary.ITER'
#     lines = get_lines(input_file)
#         
#     var = 'scan_trace_time'
#     val = read_var_from_fortran_star_file(lines, var, 'float', separator = '')
#     print(var, ' = ', val)
#         
#     var = 'scan_id'
#     val = read_var_from_fortran_star_file(lines, var, 'string', separator = '')
#     print(var, ' = ', val)
#     
#     var = 'dim_v_vector'
#     val = read_var_from_fortran_star_file(lines, var, 'int', separator = '')
#     print(var, ' = ', val)
#     
#     var = 'scan_date_v'
#     val = read_var_from_fortran_star_file(lines, var, 'int', separator = '')
#     print(var, ' = ', val)
#     
#     var = 'ray_stop_flag_run'
#     val = read_var_from_fortran_star_file(lines, var, 'string', separator = '')
#     print(var, ' = ', val)
# 
#     var = 'end_ray_vec_run'
#     val = read_var_from_fortran_star_file(lines, var, 'float', separator = '')
#     print(var, ' = ', val)


