#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for V3FIT component. This driver only uses the VMEC component to
#  stage in state files. The v3fit components generate the actual wout file.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
import os
import shutil
from ips_component_utilities import ZipState
from ips_component_utilities import ScreenWriter

#-------------------------------------------------------------------------------
#
#  V3FIT Driver Constructor
#
#-------------------------------------------------------------------------------
class v3fit_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        self.eq_worker = {'sim_name': None, 'init': None, 'driver': None}

#-------------------------------------------------------------------------------
#
#  V3FIT Driver init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'v3fit_driver: init')
        
#  Separate out the siesta, vmec and v3fit keywords.
        eq_keywords = {}
        v3fit_keywords = {}
        for key, value in keywords.items():
            if 'vmec__' in key or 'siesta__' in key:
                eq_keywords[key] = value
            if 'v3fit__' in key:
                v3fit_keywords[key.replace('v3fit__','',1)] = value

#  Get config filenames.
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        current_siesta_state = self.services.get_config_param('CURRENT_SIESTA_STATE')
        self.current_v3fit_state = self.services.get_config_param('CURRENT_V3FIT_STATE')

#  We need to pass the inputs to the SIESTA or VMEC child workflow.
        self.services.stage_state()
            
        zip_ref = ZipState.ZipState(self.current_v3fit_state, 'a')

#  If this is the first call, set up the VMEC or SIESTA sub workflow.
        if timeStamp == 0.0:
            if os.path.exists('eq_input_dir'):
                shutil.rmtree('eq_input_dir')
            os.mkdir('eq_input_dir')

            self.v3fit_port = self.services.get_port('V3FIT')
            
            if current_siesta_state in zip_ref:
#  Get keys for the SIESTA sub workflow.
                keys = {'PWD'              : self.services.get_config_param('PWD'),
                        'USER_INPUT_FILES' : current_siesta_state,
                        'SIM_NAME'         : '{}_siesta'.format(self.services.get_config_param('SIM_NAME')),
                        'LOG_FILE'         : 'log.{}_siesta.warning'.format(self.services.get_config_param('SIM_NAME')),
                        'OUTPUT_LEVEL'     : self.services.get_config_param('OUTPUT_LEVEL')
                       }

                try:
                    self.services.get_config_param('VMEC_NAMELIST_INPUT')
                except:
                    pass
                else:
                    keys['VMEC_NAMELIST_INPUT'] = self.services.get_config_param('VMEC_NAMELIST_INPUT')

                try:
                    self.services.get_config_param('SIESTA_NAMELIST_INPUT')
                except:
                    pass
                else:
                    keys['SIESTA_NAMELIST_INPUT'] = self.services.get_config_param('SIESTA_NAMELIST_INPUT')

                siesta_config = self.services.get_config_param('SIESTA_CONFIG')
                    
                (self.eq_worker['sim_name'],
                 self.eq_worker['init'],
                 self.eq_worker['driver']) = self.services.create_sub_workflow('siesta', siesta_config,
                                                                               keys, 'eq_input_dir')
                
            else:
#  Get keys for the VMEC sub workflow.
                keys = {'PWD'              : self.services.get_config_param('PWD'),
                        'USER_INPUT_FILES' : current_vmec_state,
                        'SIM_NAME'         : '{}_vmec'.format(self.services.get_config_param('SIM_NAME')),
                        'LOG_FILE'         : 'log.{}_vmec.warning'.format(self.services.get_config_param('SIM_NAME')),
                        'OUTPUT_LEVEL'     : self.services.get_config_param('OUTPUT_LEVEL')
                       }

                try:
                    self.services.get_config_param('VMEC_NAMELIST_INPUT')
                except:
                    pass
                else:
                    keys['VMEC_NAMELIST_INPUT'] = self.services.get_config_param('VMEC_NAMELIST_INPUT')

                vmec_config = self.services.get_config_param('VMEC_CONFIG')
                
                (self.eq_worker['sim_name'],
                 self.eq_worker['init'],
                 self.eq_worker['driver']) = self.services.create_sub_workflow('vmec', vmec_config,
                                                                               keys, 'eq_input_dir')

#  Copy new subworkflow inputs to the input directory.
        if current_siesta_state in zip_ref:
            zip_ref.extract(current_siesta_state)
            shutil.copy2(current_siesta_state, 'eq_input_dir')
        else:
            zip_ref.extract(current_vmec_state)
            shutil.copy2(current_vmec_state, 'eq_input_dir')

#  Initialize and run the equilibrium. Replace values in the V3FIT state.
        self.services.call(self.eq_worker['init'], 'init', timeStamp)
        self.services.call(self.eq_worker['driver'], 'init', timeStamp, **eq_keywords)
        self.services.call(self.eq_worker['driver'], 'step', timeStamp)
            
#  After the equilibrium has run update the state.
        self.services.stage_subflow_output_files()
        if current_siesta_state in zip_ref:
            zip_ref.write(current_siesta_state)
        else:
            zip_ref.write(current_vmec_state)
        zip_ref.close()
        self.services.update_state()

#  Initialize V3FIT.
        self.wait = self.services.call_nonblocking(self.v3fit_port, 'init',
                                                   timeStamp, **v3fit_keywords)

#-------------------------------------------------------------------------------
#
#  V3FIT Driver step method. This runs the v3fit component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'v3fit_driver: step')

#  Run V3FIT.
        self.services.wait_call(self.wait, True)
        self.services.call(self.v3fit_port, 'step', timeStamp, **keywords)

#  Prepare the output files for a super work flow. Need to remove any old output
#  files first before staging the state.
        if os.path.exists(self.OUTPUT_FILES):
            os.remove(self.OUTPUT_FILES)
        self.services.stage_state()
        
#  The super flow may need to rename the output file. Check if the current state
#  matches the output file. If it does not rename the state so it can be staged.
        if not os.path.exists(self.OUTPUT_FILES):
            os.rename(self.current_v3fit_state, self.OUTPUT_FILES)

#-------------------------------------------------------------------------------
#
#  V3FIT Driver finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'v3fit_driver: finalize')
        
        self.wait = [
                     self.services.call_nonblocking(self.eq_worker['init'], 'finalize', timeStamp),
                     self.services.call_nonblocking(self.eq_worker['driver'], 'finalize', timeStamp),
                     self.services.call_nonblocking(self.v3fit_port, 'finalize', timeStamp)
                    ]
            
        self.services.wait_call_list(self.wait, True)
