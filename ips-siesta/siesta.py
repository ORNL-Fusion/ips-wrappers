#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for SIESTA component. This wapper only takes a SIESTA input file
#  and runs SIESTA.
#
#-------------------------------------------------------------------------------

from component import Component
import os

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
#  SIESTA Component Constructor
#
#-------------------------------------------------------------------------------
class siesta(Component):
    def __init__(self, services, config):
        print('siesta: Construct')
        Component.__init__(self, services, config)
        
        self.siesta_exe = self.SIESTA_EXE

#-------------------------------------------------------------------------------
#
#  SIESTA Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        print('siesta: init')
        self.services.stage_plasma_state()

        current_siesta_namelist = self.services.get_config_param('CURRENT_SIESTA_NAMELIST')
    
        params = keywords['name_list_params']
    
#  Need to set the run parameters.
    
#  Start by reading in all the lines and closing the file.
        siesta_namelist = open(current_siesta_namelist, 'r')
        siesta_namelist_lines = siesta_namelist.readlines()
        siesta_namelist.close()
        
        siesta_ladd_pert = False
        siesta_lresistive = False
        siesta_lrestart = False
        
        #  Reopen for writting.
        siesta_namelist = open(current_siesta_namelist, 'w')
        
        for line in siesta_namelist_lines:
#  Name list input files can have strings containing path separators. Check for
#  an equals sign to avoid these.
            end_found, short_line = contains_end(line)
            if (end_found):
                siesta_namelist.write(short_line)
                siesta_namelist.write('\n!  SIESTA params\n')
                for key, value in params.iteritems():
                    siesta_namelist.write('%s = %s\n'%(key, value))
                siesta_namelist.write('/\n')
                break
            elif ('!  SIESTA params\n' in line):
                siesta_namelist.write(line)
                for key, value in params.iteritems():
                    siesta_namelist.write('%s = %s\n'%(key, value))
                siesta_namelist.write('/\n')
                break
            else:
                siesta_namelist.write(line)
        
        siesta_namelist.close()

#-------------------------------------------------------------------------------
#
#  SIESTA Component step method. This runs siesta.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('siesta: step')
        
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.SIESTA_EXE,
                                            logfile = 'siesta.log')
    
        if (self.services.wait_task(task_id)):
            self.services.error('siesta: step failed.')

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  SIESTA Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('siesta: finalize')
