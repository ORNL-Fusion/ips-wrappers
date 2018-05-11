#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for SIESTA component. This driver only runs the SIESTA component.
#
#-------------------------------------------------------------------------------

from component import Component
import os
import shutil
from utilities import ZipState

#-------------------------------------------------------------------------------
#
#  SIESTA Driver Constructor
#
#-------------------------------------------------------------------------------
class siesta_driver(Component):
    def __init__(self, services, config):
        print('siesta_driver: Construct')
        Component.__init__(self, services, config)
        self.async_queue = {}
        self.ports = {}
        self.model_workers = {}
        self.first_run = True
    
#-------------------------------------------------------------------------------
#
#  SIESTA Driver init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        print('siesta_driver: init')
    
#  Separate out the siesta and vmec keywords.
        siesta_keywords = {}
        vmec_keywords = {}
        for key, value in keywords.iteritems():
            if 'vmec__' in key:
                vmec_keywords[key.replace('vmec__','',1)] = value
            if 'siesta__' in key:
                siesta_keywords[key.replace('siesta__','',1)] = value

#  Get config filenames.
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        self.current_siesta_state = self.services.get_config_param('CURRENT_SIESTA_STATE')

#  We need to pass the inputs to the VMEC child workflow.
        self.services.stage_plasma_state()
        
        zip_ref = ZipState.ZipState(self.current_siesta_state, 'a')
        zip_ref.extract(current_vmec_state)

#  If this is the first call, set up the VMEC sub workflow.
        if self.first_run:

#  Get keys for the sub workflow.
            keys = {}
            keys['VMEC_CONFIG'] = self.services.get_config_param('VMEC_CONFIG')
            keys['PWD'] = self.services.get_config_param('PWD')
            
            self.model_workers['vmec'] = {'sim_name': None, 'sub_input_dir': 'vmec_inputs'}
            if os.path.exists(self.model_workers['vmec']['sub_input_dir']):
                shutil.rmtree(self.model_workers['vmec']['sub_input_dir'])
            os.mkdir(self.model_workers['vmec']['sub_input_dir'])
            shutil.copy2(current_vmec_state, self.model_workers['vmec']['sub_input_dir'])

            (self.model_workers['vmec']['sim_name'],
             self.ports['vmec_init'],
             self.ports['vmec_driver']) = self.services.create_sub_workflow('vmec', keys['VMEC_CONFIG'], {
                                                                            'pwd'              : keys['PWD'],
                                                                            'LOG_FILE'         : 'log.vmec.warning',
                                                                            'SIM_NAME'         : 'vmec',
                                                                            'USER_INPUT_FILES' : current_vmec_state
                                                                            }, self.model_workers['vmec']['sub_input_dir'])
                                                                            
            self.first_run = False

        self.async_queue['vmec_init:init'] = self.services.call_nonblocking(self.ports['vmec_init'], 'init', timeStamp)

        self.ports['siesta'] = self.services.get_port('SIESTA')

#  Initalize and run VMEC. Replace values in the siesta state.
        self.services.wait_call(self.async_queue['vmec_init:init'], True)
        self.async_queue['vmec_driver:init'] = self.services.call_nonblocking(self.ports['vmec_driver'], 'init',
                                                                              timeStamp, **vmec_keywords)
        del self.async_queue['vmec_init:init']

        self.services.wait_call(self.async_queue['vmec_driver:init'], True)
        self.async_queue['vmec_driver:step'] = self.services.call_nonblocking(self.ports['vmec_driver'], 'step',
                                                                              timeStamp)
        del self.async_queue['vmec_driver:init']

#  After VMEC has run update the VMEC state.
        self.services.wait_call(self.async_queue['vmec_driver:step'], True)

        self.services.stage_subflow_output_files()
        zip_ref.write(current_vmec_state)
        zip_ref.close()
        self.services.update_plasma_state()
        
        self.async_queue['siesta:init'] = self.services.call_nonblocking(self.ports['siesta'], 'init',
                                                                         timeStamp, **siesta_keywords)
        del self.async_queue['vmec_driver:step']

#-------------------------------------------------------------------------------
#
#  SIESTA Driver step method. This runs the vmec component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('siesta_driver: step')

        self.services.wait_call(self.async_queue['siesta:init'], True)
        self.async_queue['siesta:step'] = self.services.call_nonblocking(self.ports['siesta'], 'step',
                                                                         timeStamp)
        del self.async_queue['siesta:init']
    
#-------------------------------------------------------------------------------
#
#  SIESTA Driver finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('siesta_driver: finalize')
        self.services.wait_call_list(self.async_queue.values(), True)
        self.async_queue = {}

        for key, value in self.ports.iteritems():
            self.async_queue['{}:finalize'.format(key)] = self.services.call_nonblocking(value, 'finalize',
                                                                                         timeStamp)

        self.services.wait_call_list(self.async_queue.values(), True)
        self.async_queue = {}
