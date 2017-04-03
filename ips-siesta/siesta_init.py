#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC init component. This wapper only takes a SIESTA input
#  file and runs SIESTA.
#
#-------------------------------------------------------------------------------

from component import Component
import shutil

#-------------------------------------------------------------------------------
#
#  SIESTA init Component Constructor
#
#-------------------------------------------------------------------------------
class vmec_init(Component):
    def __init__(self, services, config):
        print('siesta_init: Construct')
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  SIESTA init Component init method. This method prepairs the namelist input
#  file and creates a dummy out put file. This allows staging the plasma state
#  files.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('siesta_init: init')

        self.services.stage_input_files(self.INPUT_FILES)
 
#  Currently SIESTA has the namelist input file hard coded. This could change in the future.  
        current_siesta_namelist = self.services.get_config_param('CURRENT_SIESTA_NAMELIST')
#        shutil.copyfile(self.INPUT_FILES, current_siesta_namelist)
        
#  Create a dummy restart file so the plasma state has something to update to.
        open(self.services.get_config_param('CURRENT_SIESTA_RESTART_FILE'), 'a').close()
        
#  Need to set and reset some values in the siest namelist input file. Primarily
#  the vmec wout file name input needs to be reset to the to the internal name..
        current_vmec_wout_file = self.services.get_config_param('CURRENT_VMEC_WOUT_FILE')

#  Start by reading in all the lines and closing the file.
        siesta_namelist = open(current_siesta_namelist, 'r')
        siesta_namelist_lines = siesta_namelist.readlines()
        siesta_namelist.close()
        
#  Reopen the file for writing.
        siesta_namelist = open(current_siesta_namelist, 'w')
        for line in siesta_namelist_lines:
            if ('WOUT_FILE' in line):
                siesta_namelist.write('WOUT_FILE = \'%s\'\n'%(current_vmec_wout_file))
                siesta_vmec_wout_set = True
            
#  The siesta namelist input files contains strings that may have path
#  separators in them. Check for an equal sign to avoid these cases.
            elif (('/' in line) and ('=' not in line)):
                if (not siesta_vmec_wout_set):
                    siesta_namelist.write('WOUT_FILE = \'%s\'\n'%(current_vmec_wout_file))
                siesta_namelist.write('/\n')
            else:
                siesta_namelist.write(line)
    
        siesta_namelist.close()
    
        self.services.update_plasma_state()


#-------------------------------------------------------------------------------
#
#  SIESTA init Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('siesta_init: step')

#-------------------------------------------------------------------------------
#
#  SIESTA init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('siesta_init: finalize')
