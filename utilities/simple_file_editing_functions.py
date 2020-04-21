#! /usr/bin/env python

"""
simple_assignment_file_edit.py 3/13/2018 (Batchelor)
Utilities to read, modify, and write text files with lines in the form of assignment
statements, i.e. of the form <name> = <value>
For now it only deals with single line assignments.
"""

# Working notes:
# Batchelor 4/21/2020: Copied futurized component from dbb4 branch
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

#_________________________________________________________________________________________________

if __name__ == '__main__':

    VD = {'x': 1.0, 'y': 2.000001}
    lines = variable_dict_to_lines(VD)
    print('lines = ', lines)

    variable_dict = lines_to_variable_dict(lines)
    print('variable_dict = ', variable_dict)
    
    variable_dict_to_output_file(VD, 'out_file')

    VD_2 = input_file_to_variable_dict('little_dict2')
    print(VD_2)

    change_dict = {'x': [113.0, 142], 'y': 22.000001, 'switch' : "'on'"}
    modify_variables_in_file(change_dict, 'little_dict2')

    change_dict = {'Q': 600}
    add_variables_to_output_file(change_dict, 'test_file_2')
