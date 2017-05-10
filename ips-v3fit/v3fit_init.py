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
        current_vmec_wout_file = self.services.get_config_param('CURRENT_VMEC_WOUT_FILE')
        
#  Start by reading in all the lines and closing the file.
        v3fit_namelist = open(current_v3fit_namelist, 'r')
        v3fit_namelist_lines = v3fit_namelist.readlines()
        v3fit_namelist.close()
        
        v3fit_task_set = False
        v3fit_vmec_nli_set = False
        v3fit_vmec_wout_set = False
        v3fit_siesta_nli_set = False
        v3fit_siesta_restart_set = False
        
#  Reopen the file for writing.
        v3fit_namelist = open(current_v3fit_namelist, 'w')
        for line in v3fit_namelist_lines:
            if ('my_task' in line):
                v3fit_namelist.write('my_task = \'%s\'\n'%(self.V3FIT_TASK))
                v3fit_task_set = True
            if ('model_eq_type' in line):
                v3fit_namelist.write('model_eq_type = \'%s\'\n'%(self.MODEL_EQ_TYPE))
                v3fit_model_set = True
            elif ('vmec_nli_filename' in line):
                v3fit_namelist.write('vmec_nli_filename = \'%s\'\n'%(current_vmec_namelist))
                v3fit_vmec_nli_set = True
            elif ('vmec_wout_input' in line):
                v3fit_namelist.write('vmec_wout_input = \'%s\'\n'%(current_vmec_wout_file))
                v3fit_vmec_wout_set = True
            elif (self.MODEL_EQ_TYPE == 'siesta' and 'siesta_nli_filename' in line):
                v3fit_namelist.write('siesta_nli_filename = \'%s\'\n'%(self.services.get_config_param('CURRENT_SIESTA_NAMELIST')))
                v3fit_siesta_nli_set = True
            elif (self.MODEL_EQ_TYPE == 'siesta' and 'siesta_restart_filename' in line):
                v3fit_namelist.write('siesta_restart_filename = \'%s\'\n'%(self.services.get_config_param('CURRENT_SIESTA_RESTART_FILE')))
                v3fit_siesta_restart_set = True

#  The v3fit namelist input files contains strings that may have path separators
#  in them. Check for an equal sign to avoid these cases.
            elif (('/' in line) and ('=' not in line)):
                if (not v3fit_task_set):
                    v3fit_namelist.write('my_task = \'%s\'\n'%(self.V3FIT_TASK))
                if (not v3fit_model_set):
                    v3fit_namelist.write('model_eq_type = \'%s\'\n'%(self.MODEL_EQ_TYPE))
                if (not v3fit_vmec_nli_set):
                    v3fit_namelist.write('vmec_nli_filename = \'%s\'\n'%(current_vmec_namelist))
                if (not v3fit_vmec_wout_set):
                    v3fit_namelist.write('vmec_wout_input = \'%s\'\n'%(self.current_vmec_wout_file))
                if (self.MODEL_EQ_TYPE == 'siesta' and not v3fit_siesta_nli_set):
                    v3fit_namelist.write('siesta_nli_filename = \'%s\'\n'%(self.services.get_config_param('CURRENT_SIESTA_NAMELIST')))
                if (self.MODEL_EQ_TYPE == 'siesta' and not v3fit_siesta_restart_set):
                    v3fit_namelist.write('siesta_restart_filename = \'%s\'\n'%(self.services.get_config_param('CURRENT_SIESTA_RESTART_FILE')))
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
