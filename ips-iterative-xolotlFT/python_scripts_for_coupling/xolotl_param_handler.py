#!/usr/bin/env python
#  Example file for a Xolotl class object.

# Check if the value can be represented as a float. If the string can be vaildly
# cast as a float value, then it must be a float otherwise an exception is
# thrown.
# Inputs:
# value      Value to test.
# Returns:
# boolean    True if a floating value false otherwise.
def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

# Check if the value can be represented as an int. If the string can be vaildly
# cast as a in value, then it must be a float otherwise an exception is thrown.
# Inputs:
# value    Value to test.
# Returns:
# boolean    True if a integer value false otherwise.
def is_int(value):
    try:
        int(value)
        return True
    except ValueError:
        return False

# Determine if the value is a list. Elements in the list are separated by a ' '
# character. If any ' ' is in the string, it must be a list.
# Inputs:
# value    Value to test.
# Returns:
# boolean    True if a list false otherwise.
def is_list(value):
    return ' ' in value

def write(file_handle, key, value):
    
    import math
    #remove space between key and value (Xolotl doesn't like it)
    #{0:s} = {1:g}\n --> {0:s}={1:g}\n
    if isinstance(value, float):
        file_handle.write('{0:s}={1:g}\n'.format(key, value))
    elif isinstance(value, int):
        file_handle.write('{0:s}={1:d}\n'.format(key, value))
    elif isinstance(value, str):
        file_handle.write('{0:s}={1:s}\n'.format(key, value))
    elif isinstance(value, list):
        file_handle.write('{0:s}='.format(key))
        for v in value:
            if isinstance(v, int):
                file_handle.write('{0:d} '.format(v))
            elif isinstance(v, float):
                file_handle.write('{0:g} '.format(v))
            else:
                file_handle.write('{0:s} '.format(v))
        file_handle.write('\n')
    elif isinstance(value, dict):
        file_handle.write('{0:s}='.format(key))
        for k, v in value.iteritems():
            if v != None:
                if isinstance(v, int):                  
                    file_handle.write('{0:s} {1:d} '.format(k, v))
                elif isinstance(v, float):
                    #inf != INFINITY in petscArg
                    if math.isinf(v):
                        v='INFINITY'
                        file_handle.write('{0:s} {1:s} '.format(k, v))                    
                    else:
                        file_handle.write('{0:s} {1:g} '.format(k, v))
                else:
                    file_handle.write('{0:s} {1:s} '.format(k, v))
            else:
                file_handle.write('{0:s} '.format(k))
        file_handle.write('\n')

# This class represents the parameters of a xolotl input file.
class xolotl_params:
    
    #  Dictionary to hold parameters.
    parameters = {}

    # Read an input file and populate the parameters dictionary.
    # Inputs:
    # file_name    Path to the file to read.
    def read(self, file_name):
        file_handle = open(file_name, 'r')
        lines = file_handle.readlines()
        file_handle.close()

        # Loop over every line to parse the input.
        for line in lines:
            # Remove the string newline.
            line = line[:-1]
            
            # An equal sign separate the key and the value.
            key, value = line.split('=')

            # Check the value for different types. Currently its a string but
            # we want it work in a more basic data type.
            if is_int(value):
                self.parameters[key] = int(value)
            elif is_float(value):
                self.parameters[key] = float(value)
            else:
                # Since it is not a simple integer, parse the string to figure
                # out if it represents a more complex type.
                self.parse_string(key, value)

    # Handle the cases that are not simple integers or floats.
    # Inputs:
    # key       For the current string value.
    # value     String to parse.
    def parse_string(self, key, value):
        if (key == 'petscArgs'):
            # petscArgs is best represented as a dictionary. Parse the value for
            # this key as a special case. parse_petscargs expects a list of
            # strings so split it before the call.
            args = {}
            self.parse_petscargs(args, value.split(' '))
            self.parameters[key] = args
            return
        elif is_list(value):
            # If this is a list, cast each element of the list to its basic
            #type.
            values = []
            for v in value.split(' '):
                if is_int(v):
                    values.append(int(v))
                elif is_float(v):
                    values.append(float(v))
                else:
                    values.append(v)
            self.parameters[key] = values
        else:
            # The value must be a string. Store it directly.
            self.parameters[key] = value

    # Recursive method to parse all the elements of the petscArgs key. This
    # method assumes that the a petsc argument only takes 0 or 1 value.
    # Inputs:
    # args       List of arguments already parsed.
    # values     List of strings to classify.
    def parse_petscargs(self, args, values):
        if len(values) == 0:
            # There are no more values left to classify.
            return
        
        # At least one value remains assume it is an argument. Then check the
        # next value to see if it is a string or the next argument. If the next
        # value is a float, integer, or plain string, associate that value with
        # key. If the next value is a new flag, the current argument takes no
        # value.
        arg = values.pop(0)
        if len(values) > 0:
            if is_int(values[0]):
                args[arg] = int(values.pop(0))
            elif is_float(values[0]):
                args[arg] = float(values.pop(0))
            elif values[0][0] != '-':
                args[arg] = values.pop(0)
            else:
                args[arg] = None
        else:
            args[arg] = None

        # Recursively classify the remaining values.
        self.parse_petscargs(args, values)

    # Write a new input file. This loops through all the values in the
    # dictionary and formats them into an input file.
    # Inputs:
    # file_name    Path to the file to write.
    def write(self, file_name):
        file_handle = open(file_name, 'w')
        
        for key, value in self.parameters.iteritems():
            write(file_handle, key, value)

        file_handle.close()

        #strip end of line space
        file_handle = open(file_name, 'r')

        for line in file_handle:
            line.strip()

        file_handle.close()

if __name__ == '__main__':
    xp = xolotl_params()
    xp.read('params.txt')
    
    # examples on changing values.
    xp.parameters['petscArgs']['-ts_adapt_dt_max'] = 0.24
    xp.parameters['vizHandler'] = 1.0
    xp.parameters['startTemp'] = 773.5

    xp.write('test.txt')

