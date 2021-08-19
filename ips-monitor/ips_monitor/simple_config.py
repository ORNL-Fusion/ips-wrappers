import sys

"""
simple_config.py

A class that implements a very simple configuration file parser.  It doesn't have
recursive section levels nor does it support variable interpolation.

It recognizes three kinds of lines in the config file:

1) Lines starting with # (possibly preceeded by white space) are ignored as comments

2) Lines of the form [<string>] are regarded as beginning a section named <string>.

3) Lines of the from <variable_name> = <variable_value_string> are variable definitions.

   Coninuation - If the last character of the line is \ the next line is appended to the
   variable_value_string as a continuation line
   
   If the first and last characters of variable_value_string are
   [ and ] respectively, the value is interpreted as a comma delimited list.  The
   value is then converted to a Python list of strings.
   
   It is possible include lists of strings which might contain white space. In a list, if 
   the first character following the initial [, and the last character before the final ]  is
   ' or ", the line is parsed as a ' or " delimited list of strings.  So you can include ' or "
   in a string by using the other one as the delimiter.  As of now it only knows how to deal
   with one space between successive quoted strings

An instance of simple_config has the following attributes
1) sections[] -> a list of the names of sections that it contains, if any
2) variables[] -> a list of names of top level variables, if any
3) Top level variable attributes (maybe none) <variable_name> = <variable_value_string>
4) Instances of class simple_config_section, one for each section contained (maybe none)

An instance of simple_config_section has the following attributes
1) variables[] -> a list of names of variable attributes contained
2) Variable attributes (maybe none) <variable_name> = <variable_value_string>

Notes:

Variable values are strings.  If conversion to numbers or something else is needed it's up
to the user to do the conversion.

Indentation is not significant.  White space on the left is deleted before parsing.

For contined lines a space is added between the previous line and the contuation line so
that list items don't get welded together if there is no space.  Therefore:
-> Don't break a oontinued line in the middle of a word.

"""

__all__ = ['simple_config']

debug = False

class simple_config_section:
    def __init__(self):
        self.variables = []

class simple_config:
    def __init__(self, path):
        print('Created %s' % (self.__class__))
        self.config_file = path
        if debug:
            print("init: config_file = ", self.config_file)
        
        self.sections = []
        self.variables = []

# Read the file
        file = open(self.config_file, 'r')
        self.lines = file.readlines()
        file.close()
        
        if debug :
            print("   lines =")
            for line in self.lines:
                print(line)

# Parse the lines
        section_name = 'self'
        continuation_line = False

        for line in self.lines:
            ln = line.lstrip()
            
            #Check if this is a blank line.  If so skip
            if len(ln) == 0:
                continue
            
            # Check if this is a continuation line, if so add to previous line
            if continuation_line:
                ln = old_line + ' ' + ln
                
                if debug:
                    print('continued line = \n', ln)
                    
            # Check if this line is to be continued. If so save line up to '\' and get the next one. 
            # Do the parsing after all continuations have been added
            if ln[-2] == '\\' :
                old_line = ln[:-2]
                continuation_line = True
                continue
            else:
                continuation_line = False

            # Now parse
            if ln[0] == '#':    # skip commnets
                continue
            
            # Check if start of new section. Create a new instance of simple_config_section
            elif ln[0] == '[' and ln[-2] == ']':
                section_name = ln[1:-2]
                self.sections.append(section_name)
                section = vars(self)[section_name] = simple_config_section()
                                
                if debug:
                    print('new section -> ', section_name)
            
            # This has got to be a variable definition
            else:
                var_line = ln.partition('=')
                
                if var_line[1] != '=':  # Doesn't have form of variable def
                    print('simple_config: line not in form of a variable def')
                    print('var_line = ', var_line)
                    raise Exception('simple_config: line not in form of a variable def')
                    
                var_name = var_line[0].strip()
                var_value = var_line[2].strip()
                
                if debug:
                    print('len(var_value) = ', len(var_value), '  var_value = ', var_value)
                
               # Check if the value is a list.  In which case convert to a Python list.
               # It can't be a list unless its length >= 2
                if len(var_value) >= 2:                 
                    
                    if var_value[0] == '[' and var_value[-1] == ']':
                        inner_string = var_value[1:-1]
                        
                        var_value = inner_string.split(',') 
                        # For good measure strip white space from list elements
                        var_value = [x.strip() for x in var_value]
                
                if section_name == 'self':  # This is a top level Global variable
                    vars(self)[var_name] = var_value
                    self.variables.append(var_name)
                else:
                    section.variables.append(var_name)
                    vars(section)[var_name] = var_value
                    print('section varibles: ', var_name, ' = ', var_value)
                    

#---------------------------------------------------------------------------------------
#
# main
#
#---------------------------------------------------------------------------------------

if __name__ == '__main__':

    print("Open a config file")
    file_name = 'dbb.config'
    config = simple_config(file_name)
    
    print('\n dir(config) = ', dir(config))
    
    print('\n config.variables = ', config.variables)
    
    print('\n config.sections = ', config.sections)
    
    print('\n config.x = ', config.x)
    
    print('\n config.y = ', config.y)
    
    print('\n dir(config.section_1) = ', dir(config.section_1))
        
    print('\n config.section_1.variables = ', config.section_1.variables)
        
    print('\n config.section_1.x = ', config.section_1.x)
    
    print('\n config.section_1.z = ', config.section_1.z)

    print('\n config.section_2.p = ', config.section_2.p)
    for x in config.section_2.p:
        print(x)

    print('\n config.section_2.q = ', config.section_2.q)
    for x in config.section_2.q:
        print(x)
    
    
    
