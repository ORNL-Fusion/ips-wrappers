#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC component. This wapper only takes a VMEC input file and
#  runs VMEC.
#
#-------------------------------------------------------------------------------

from component import Component
import os
from omfit.classes.omfit_namelist import OMFITnamelist
from utilities import ZipState
import json

#-------------------------------------------------------------------------------
#
#  VMEC Component Constructor
#
#-------------------------------------------------------------------------------
class vmec(Component):
    def __init__(self, services, config):
        print('vmec: Construct')
        Component.__init__(self, services, config)
        self.task_queue = {}

#-------------------------------------------------------------------------------
#
#  VMEC Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        print('vmec: init')
        
        self.current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        self.current_wout_file = 'wout_{}.nc'.format(self.current_vmec_namelist.replace('input.','',1))
        self.current_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        current_vmec_flags = self.services.get_config_param('CURRENT_VMEC_FLAGS')
        
#  Stage plasma state.
        self.services.stage_plasma_state()

#  Unzip files from the plasma state. Use mode a so files can be read and
#  written to.
        self.zip_ref = ZipState.ZipState(self.current_state, 'a')
        self.zip_ref.extract(self.current_vmec_namelist)
        self.zip_ref.extract('flags.json')
        
#  Load the json file.
        with open('flags.json', 'r') as flags_file:
            self.flags = json.load(flags_file)

        if len(keywords) > 0:
            self.flags['state'] = 'needs_update'
        
#  Update parameters in the namelist.
            namelist = OMFITnamelist(self.current_vmec_namelist)

            for key, value in keywords.iteritems():
                if '(' in key :
                    key, indices = key.split('(')
                    indices, extra = indices.split(')')
                    indices = [[int(i) - 1] for i in indices.split(',')]
                    namelist['indata'][key][indices] = value
                else:
                    namelist['indata'][key] = value
    
            namelist.save()

#-------------------------------------------------------------------------------
#
#  VMEC Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('vmec: step')

        if self.flags['state'] == 'needs_update':
            self.task_queue['vmec'] = self.services.launch_task(self.NPROC,
                                                                self.services.get_working_dir(),
                                                                self.VMEC_EXE,
                                                                self.current_vmec_namelist,
                                                                logfile = 'vmec.log')

#  Update flags.
            self.flags['state'] = 'updated'
            with open('flags.json', 'w') as flags_file:
                json.dump(self.flags, flags_file)

#  Wait for VMEC to finish.
            if (self.services.wait_task(self.task_queue['vmec'], True) or not os.path.exists(self.current_wout_file)):
                self.services.error('vmec: step failed.')
            del self.task_queue['vmec']

#  Add the wout file to the plasma state.
            self.zip_ref.write([self.current_vmec_namelist, self.current_wout_file, 'flags.json'])
            self.zip_ref.close()

            self.services.update_plasma_state()
            self.services.stage_output_files(timeStamp, self.OUTPUT_FILES,
                                             keep_old_files=False)
        else:
#  Update flags.
            self.flags['state'] = 'unchanged'
            with open('flags.json', 'w') as flags_file:
                json.dump(self.flags, flags_file)
            
            self.zip_ref.write(['flags.json'])
            self.zip_ref.close()
            
            self.services.update_plasma_state()
            self.services.stage_output_files(timeStamp, self.OUTPUT_FILES,
                                             keep_old_files=False)

#-------------------------------------------------------------------------------
#
#  VMEC Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('vmec: finalize')
        self.services.wait_tasklist(self.task_queue.values(), True)
