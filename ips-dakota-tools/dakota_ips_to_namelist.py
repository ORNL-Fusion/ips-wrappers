#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper to bridge the DAKOTA params to a namelist input file.
#
#-------------------------------------------------------------------------------

from component import Component

#-------------------------------------------------------------------------------
#
#  Helper function to detect the end of a namelist input file. This finds if the
#  '/' character is found withing a string. However, this if this character is
#  is located after a '!' character or between '\'' characters ignore it.
#
#-------------------------------------------------------------------------------
def contains_end(line):

#  Check and remove every character after the !.
    line = line.split('!')[0]

    if (line == ''):
        return False, ''

    if (line[0] == '/'):
        return True, ''

    in_single_quote = line[0] == '\''
    in_double_quote = line[0] == '\"'
    
    for i, c in enumerate(line[1:]):

#  Do not close the quote when a single quote is encounterd in the following
#  situations.
#
#  '\''
#  "'"
        if ((c == '\'') and line[i - 1] != '\\' and not in_double_quote):
            in_single_quote = not in_single_quote
        
#  Do not close the quote when a single quote is encounterd in the following
#  situations.
#
#  "\""
#  '"'
        elif ((c == '\"') and line[i - 1] != '\\' and not in_single_quote):
            in_double_quote = not in_double_quote
        elif ((c == '/') and not in_double_quote and not in_single_quote):
            return True, line[:i]
        elif ((c == '&') and not in_double_quote and not in_single_quote):
            if (i + 4 <= len(line)) and (line[i:i + 4] == '&end'):
                return TRUE, line[:i]

    return False, ''

#-------------------------------------------------------------------------------
#
#  DAKOTA to fortran namelist Component Constructor
#
#-------------------------------------------------------------------------------
class dakota_ips_to_namelist(Component):
    def __init__(self, services, config):
        print('dakota_ips_to_namelist: Construct')
        Component.__init__(self, services, config)
    
#-------------------------------------------------------------------------------
#
#  DAKOTA to fortran namelist Component init method. This method reads the
#  parameter file from datoka and adds to an existing namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('dakota_ips_to_namelist: init')
        
        self.services.stage_plasma_state()
        
#  Parse the environment variables for the prefix.
        params = {}
    
        for key, value in self.__dict__.iteritems():
            if (self.PREFIX in key):
                params[key.replace(self.PREFIX, '')] = value
    
#  Append the new parameters to the end of the namelist input file before the
#  end. Mark the start of the dakota parameter with a comment. To avoid
#  littering the namelist file. NAMELIST_FILE should reference a file in the
#  plasma state.
        namelist_file = open(self.NAMELIST_FILE, 'r')
        namelist_file_lines = namelist_file.readlines()
        namelist_file.close()
        
#  Reopen for writting.s
        namelist_file = open(self.NAMELIST_FILE, 'w')
        
        todem_found = False
        for line in namelist_file_lines:
#  Name list input files can have strings containing path separators. Check for
#  an equals sign to avoid these.
            end_found, short_line = contains_end(line)
            if (end_found and not todem_found):
                namelist_file.write(short_line)
                namelist_file.write('\n!  DAKOTA params\n')
                for key, value in params.iteritems():
                    namelist_file.write('%s = %s\n'%(key, value))
                namelist_file.write('/')
            elif ('!  DAKOTA params\n' in line):
                todem_found = True
                namelist_file.write(line)
                for key, value in params.iteritems():
                    namelist_file.write('%s = %s\n'%(key, value))
            else:
                namelist_file.write(line)
                    
        namelist_file.close()

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  DAKOTA to fortran namelist Component step method. Not used.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('dakota_ips_to_namelist: step')
    
#-------------------------------------------------------------------------------
#
#  DAKOTA to fortran namelist Component finalize method. This cleans up afterwards. Not
#  used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('dakota_ips_to_namelist: finalize')
