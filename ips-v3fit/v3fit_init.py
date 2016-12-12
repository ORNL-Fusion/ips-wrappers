#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for V3FIT init component. This wapper only takes a V3FIT input
#  file and runs V3FIT.
#
#-------------------------------------------------------------------------------

from component import Component
import shutil

#-------------------------------------------------------------------------------
#
#  V3FIT init Component Constructor
#
#-------------------------------------------------------------------------------
class v3fit_init(Component):
    def __init__(self, services, config):
        print('v3fit_init: Construct')
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  V3FIT init Component init method. This method prepairs the namelist input
#  file and creates a dummy out put file. This allows staging the plasma state
#  files. In the v3fit namelist input file configure the v3fit namelist input
#  with the task, internal vmec input name and optional name of the wout file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('v3fit_init: init')
    
        self.services.stage_input_files(self.INPUT_FILES)

        current_v3fit_namelist = self.services.get_config_param('CURRENT_V3FIT_NAMELIST')
        shutil.copyfile(self.INPUT_FILES, current_v3fit_namelist)

#  Create a dummy result file so the plasma state has something to update to.
        open(self.services.get_config_param('CURRENT_V3FIT_RESULT_FILE'), 'a').close()

#  Need to set and reset some values in the v3fit namelist input file. Primarily
#  the vmec namelist input needs to be reset to the to the internal name and the
#  v3fit task must be set. Optionally a wout file may be specified.

#  All V3FIT namelist input files need a VMEC file name.
        current_vmec_namelist = self.services.get_config_param('CURRENT_VMEC_NAMELIST')

#  Start by reading in all the lines and closing the file.
        v3fit_namelist = open(current_v3fit_namelist, 'r')
        v3fit_namelist_lines = v3fit_namelist.readlines()
        v3fit_namelist.close()
        
        v3fit_task_set = False
        v3fit_vmec_nli_set = False
        v3fit_vmec_wout_set = False
        
#  Reopen the file for writing.
        v3fit_namelist = open(current_v3fit_namelist, 'w')
        for line in v3fit_namelist_lines:
            if ('my_task' in line):
                v3fit_namelist.write('my_task = \'%s\'\n'%(self.V3FIT_TASK))
                v3fit_task_set = True
            elif ('vmec_nli_filename' in line):
                v3fit_namelist.write('vmec_nli_filename = \'%s\'\n'%(current_vmec_namelist))
                v3fit_vmec_nli_set = True
            elif ('vmec_wout_input' in line):
                v3fit_namelist.write('vmec_wout_input = \'%s\'\n'%(self.V3FIT_WOUT_FILE))
                v3fit_vmec_wout_set = True
            elif ('/' in line):
                if (not v3fit_task_set):
                    v3fit_namelist.write('my_task = \'%s\'\n'%(self.V3FIT_TASK))
                if (not v3fit_vmec_nli_set):
                    v3fit_namelist.write('vmec_nli_filename = \'%s\'\n'%(current_vmec_namelist))
                if (not v3fit_vmec_wout_set):
                    v3fit_namelist.write('vmec_wout_input = \'%s\'\n'%(self.V3FIT_WOUT_FILE))
                v3fit_namelist.write('/\n')
            else:
                v3fit_namelist.write(line)

        v3fit_namelist.close()

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  V3FIT init Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('v3fit_init: step')

#-------------------------------------------------------------------------------
#
#  V3FIT init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('v3fit_init: finalize')
