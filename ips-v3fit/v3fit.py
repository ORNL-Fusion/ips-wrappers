#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for V3FIT component. This wapper only takes a V3FIT input file
#  and runs V3FIT.
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
#  V3FIT Component Constructor
#
#-------------------------------------------------------------------------------
class v3fit(Component):
    def __init__(self, services, config):
        print('v3fit: Construct')
        Component.__init__(self, services, config)
        
        self.v3fit_exe = self.V3FIT_EXE

#-------------------------------------------------------------------------------
#
#  V3FIT Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        print('v3fit: init')
        self.services.stage_plasma_state()
        
        self.current_v3fit_namelist = self.services.get_config_param('CURRENT_V3FIT_NAMELIST')

        params = keywords['name_list_params']

#  Append the new parameters to the end of the namelist input file before the
#  end. Mark the start of the dakota parameter with a comment. To avoid
#  littering the namelist file. NAMELIST_FILE should reference a file in the
#  plasma state.

        v3fit_namelist = open(self.current_v3fit_namelist, 'r')
        v3fit_namelist_lines = v3fit_namelist.readlines()
        v3fit_namelist.close()
        
        #  Reopen for writting.
        v3fit_namelist = open(self.current_v3fit_namelist, 'w')

        for line in v3fit_namelist_lines:
#  Name list input files can have strings containing path separators. Check for
#  an equals sign to avoid these.
            end_found, short_line = contains_end(line)
            if (end_found):
                v3fit_namelist.write(short_line)
                v3fit_namelist.write('\n!  V3FIT params\n')
                for key, value in params.iteritems():
                    v3fit_namelist.write('%s = %s\n'%(key, value))
                v3fit_namelist.write('/\n')
                break
            elif ('!  V3FIT params\n' in line):
                v3fit_namelist.write(line)
                for key, value in params.iteritems():
                    v3fit_namelist.write('%s = %s\n'%(key, value))
                v3fit_namelist.write('/\n')
                break
            else:
                v3fit_namelist.write(line)

        v3fit_namelist.close()

#-------------------------------------------------------------------------------
#
#  V3FIT Component step method. This runs V3FIT.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('v3fit: step')

        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            self.V3FIT_EXE,
                                            self.current_v3fit_namelist,
                                            logfile = 'v3fit.log')
    
        if (self.services.wait_task(task_id)):
            self.services.error('v3fit: step failed.')

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  V3FIT Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('v3fit: finalize')
