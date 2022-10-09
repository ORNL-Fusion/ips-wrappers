#! /usr/bin/env python

"""
simple_assignment_file_edit.py 3/13/2018 (Batchelor)
Utilities to read, modify, and write text files with lines in the form of assignment
statements, i.e. of the form <name> = <value>
For now it only deals with single line assignments.
"""

# Working notes:
# Batchelor 4/21/2020: Copied futurized component from dbb4 branch
# Batchelor 2/1/2021: Added read_var_from_nml_lines
#

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
    file = open(filename, 'w')
    file.writelines(lines)
    file.close()

#---------------------------------------------------------------------------------------
# Editing utilities
#---------------------------------------------------------------------------------------

def lines_to_variable_dict(lines):
    variable_dict = {}
    for line in lines:
        if (line.strip() not in ['!', '#']) and (len(line) > 0):  # ignore comments
            if '=' in line:
                name, val = line.split("=")
                name = name.strip()
                val = val.strip()
                variable_dict[name] = val
    return variable_dict

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

#_________________________________________________________________________________________________

if __name__ == '__main__':

    VD = {'x': 1.0, 'y': 2.000001}
    lines = variable_dict_to_lines(VD)
    print('lines = ', lines)

    variable_dict = lines_to_variable_dict(lines)
    print('variable_dict = ', variable_dict)
    
    variable_dict_to_output_file(VD, 'out_file')
    
    read_var_from_nml_lines(lines, 'x', separator = ',')
    
    VD_2 = input_file_to_variable_dict('little_dict2')
    print(VD_2)

    change_dict = {'x': [113.0, 142], 'y': 22.000001, 'switch' : "'on'"}
    modify_variables_in_file(change_dict, 'little_dict2')

    change_dict = {'Q': 600}
    add_variables_to_output_file(change_dict, 'test_file_2')

    read_var_from_nml_lines(lines, x, separator = ',')
