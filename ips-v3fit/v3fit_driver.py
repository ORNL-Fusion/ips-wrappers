#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for V3FIT component. This driver only uses the VMEC component to
#  stage in plasma state files. The v3fit components generate the actual wout
#  file.
#
#-------------------------------------------------------------------------------

from component import Component
import os
import shutil
from utilities import ZipState

#-------------------------------------------------------------------------------
#
#  V3FIT Driver Constructor
#
#-------------------------------------------------------------------------------
class v3fit_driver(Component):
    def __init__(self, services, config):
        print('v3fit_driver: Construct')
        Component.__init__(self, services, config)
        self.async_queue = {}
        self.ports = {}
        self.model_workers = {}
        self.first_run = True

#-------------------------------------------------------------------------------
#
#  V3FIT Driver init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        print('v3fit_driver: init')
        
#  Separate out the siesta, vmec and v3fit keywords.
        eq_keywords = {}
        v3fit_keywords = {}
        for key, value in keywords.iteritems():
            if 'vmec__' in key or 'siesta__' in key:
                eq_keywords[key] = value
            if 'v3fit__' in key:
                v3fit_keywords[key.replace('v3fit__','',1)] = value

#  Get config filenames.
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        current_siesta_state = self.services.get_config_param('CURRENT_SIESTA_STATE')
        self.current_v3fit_state = self.services.get_config_param('CURRENT_V3FIT_STATE')

#  We need to pass the inputs to the SIESTA or VMEC child workflow.
        self.services.stage_plasma_state()
            
        zip_ref = ZipState.ZipState(self.current_v3fit_state, 'a')

#  If this is the first call, set up the VMEC or SIESTA sub workflow.
        if self.first_run:
            if current_siesta_state in zip_ref:
                zip_ref.extract(current_siesta_state)
                    
#  Get keys for the SIESTA sub workflow.
                keys = {}
                keys['SIESTA_CONFIG'] = self.services.get_config_param('SIESTA_CONFIG')
                keys['PWD'] = self.services.get_config_param('PWD')
                    
                self.model_workers['eq'] = {'sim_name': None, 'input_dir': 'siesta_inputs'}
                if os.path.exists(self.model_workers['eq']['input_dir']):
                    shutil.rmtree(self.model_workers['eq']['input_dir'])
                os.mkdir(self.model_workers['eq']['input_dir'])
                shutil.copy2(current_siesta_state, self.model_workers['eq']['input_dir'])
                    
                (self.model_workers['eq']['sim_name'],
                 self.ports['eq_init'],
                 self.ports['eq_driver']) = self.services.create_sub_workflow('siesta', keys['SIESTA_CONFIG'], {
                                                                              'pwd'              : keys['PWD'],
                                                                              'LOG_FILE'         : 'log.siesta.warning',
                                                                              'SIM_NAME'         : 'siesta',
                                                                              'USER_INPUT_FILES' : current_siesta_state
                                                                              }, self.model_workers['eq']['input_dir'])
                
            else:
                zip_ref.extract(current_vmec_state)

#  Get keys for the VMEC sub workflow.
                keys = {}
                keys['VMEC_CONFIG'] = self.services.get_config_param('VMEC_CONFIG')
                keys['PWD'] = self.services.get_config_param('PWD')

                self.model_workers['eq'] = {'sim_name': None, 'input_dir': 'vmec_inputs'}
                if os.path.exists(self.model_workers['eq']['input_dir']):
                    shutil.rmtree(self.model_workers['eq']['input_dir'])
                os.mkdir(self.model_workers['eq']['input_dir'])
                shutil.copy2(current_vmec_state, self.model_workers['eq']['input_dir'])
                
                (self.model_workers['eq']['sim_name'],
                 self.ports['eq_init'],
                 self.ports['eq_driver']) = self.services.create_sub_workflow('vmec', keys['VMEC_CONFIG'], {
                                                                              'pwd'              : keys['PWD'],
                                                                              'LOG_FILE'         : 'log.vmec.warning',
                                                                              'SIM_NAME'         : 'vmec',
                                                                              'USER_INPUT_FILES' : current_vmec_state
                                                                              }, self.model_workers['eq']['input_dir'])
                    
            self.first_run = False

        self.async_queue['eq_init:init'] = self.services.call_nonblocking(self.ports['eq_init'], 'init',
                                                                          timeStamp)
        
        self.ports['v3fit'] = self.services.get_port('V3FIT')
            
#  Initialize and run the equilibrium. Replace values in the V3FIT state.
        self.services.wait_call(self.async_queue['eq_init:init'], True)
        self.async_queue['eq_driver:init'] = self.services.call_nonblocking(self.ports['eq_driver'], 'init',
                                                                            timeStamp, **eq_keywords)
        del self.async_queue['eq_init:init']
        
        self.services.wait_call(self.async_queue['eq_driver:init'], True)
        self.async_queue['eq_driver:step'] = self.services.call_nonblocking(self.ports['eq_driver'], 'step',
                                                                            timeStamp)
        del self.async_queue['eq_driver:init']
            
#  After the equilibrium has run update the state.
        self.services.wait_call(self.async_queue['eq_driver:step'], True)
        
        self.services.stage_subflow_output_files()
        if current_siesta_state in zip_ref:
            zip_ref.write(current_siesta_state)
        else:
            zip_ref.write(current_vmec_state)
        zip_ref.close()
        self.services.update_plasma_state()

        self.async_queue['v3fit:init'] = self.services.call_nonblocking(self.ports['v3fit'], 'init',
                                                                        timeStamp, **v3fit_keywords)
        del self.async_queue['eq_driver:step']

#-------------------------------------------------------------------------------
#
#  V3FIT Driver step method. This runs the v3fit component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0, **keywords):
        print('v3fit_driver: step')

        self.services.wait_call(self.async_queue['v3fit:init'], True)
        self.async_queue['v3fit:step'] = self.services.call_nonblocking(self.ports['v3fit'], 'step',
                                                                        timeStamp, **keywords)
    
        del self.async_queue['v3fit:init']

        self.services.wait_call_list(self.async_queue.values(), True)
        self.async_queue = {}

#  Prepare the output files for a super work flow.
        self.services.stage_plasma_state()
        
#  The super flow may need to rename the output file. Check is the current state
#  matches if output file. If it does not rename the plasma state so it can be
#  staged.
        if not os.path.exists(self.OUTPUT_FILES):
            os.rename(self.current_v3fit_state, self.OUTPUT_FILES)

#-------------------------------------------------------------------------------
#
#  V3FIT Driver finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('v3fit_driver: finalize')
        self.services.wait_call_list(self.async_queue.values(), True)
        self.async_queue = {}
        
        for key, value in self.ports.iteritems():
            self.async_queue['{}:finalize'.format(key)] = self.services.call_nonblocking(value, 'finalize',
                                                                                         timeStamp)
    
        self.services.wait_call_list(self.async_queue.values(), True)
        self.async_queue = {}
