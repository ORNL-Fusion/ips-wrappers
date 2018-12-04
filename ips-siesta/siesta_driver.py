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
from utilities import ScreenWriter

#-------------------------------------------------------------------------------
#
#  SIESTA Driver Constructor
#
#-------------------------------------------------------------------------------
class siesta_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        self.vmec_worker = {}
    
#-------------------------------------------------------------------------------
#
#  SIESTA Driver init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'siesta_driver: init')
    
#  Separate out the siesta and vmec keywords.
        siesta_keywords = {}
        vmec_keywords = {}
        for key, value in keywords.iteritems():
            if 'vmec__' in key:
                vmec_keywords[key] = value
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
        if timeStamp == 0.0:
#  Get the siesta port
            self.siesta_port = self.services.get_port('SIESTA')

#  Get keys for the sub workflow.
            keys = {'PWD'              : self.services.get_config_param('PWD'),
                    'SIM_NAME'         : '{}_vmec'.format(self.services.get_config_param('SIM_NAME')),
                    'USER_INPUT_FILES' : current_vmec_state,
                    'LOG_FILE'         : 'log.vmec.warning',
                    'OUTPUT_LEVEL'     : self.services.get_config_param('OUTPUT_LEVEL')
                   }
            
            if os.path.exists('vmec_input_dir'):
                shutil.rmtree('vmec_input_dir')
            os.mkdir('vmec_input_dir')
            
            vmec_config = self.services.get_config_param('VMEC_CONFIG')
            
            self.vmec_worker = {'sim_name': None, 'init': None, 'driver': None}
            (self.vmec_worker['sim_name'],
             self.vmec_worker['init'],
             self.vmec_worker['driver']) = self.services.create_sub_workflow('vmec', vmec_config,
                                                                             keys, 'vmec_input_dir')

        shutil.copy2(current_vmec_state, 'vmec_input_dir')

#  Initalize and run VMEC. Replace values in the siesta state.
        self.services.call(self.vmec_worker['init'], 'init', timeStamp)
        self.services.call(self.vmec_worker['driver'], 'init', timeStamp, **vmec_keywords)
        self.services.call(self.vmec_worker['driver'], 'step', timeStamp)

#  After VMEC has run update the VMEC state.
        self.services.stage_subflow_output_files()
        zip_ref.write(current_vmec_state)
        zip_ref.close()
        self.services.update_plasma_state()

#  Initialize SIESTA.
        self.wait = self.services.call_nonblocking(self.siesta_port, 'init',
                                                   timeStamp, **siesta_keywords)

#-------------------------------------------------------------------------------
#
#  SIESTA Driver step method. This runs the vmec component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'siesta_driver: step')

#  Run SIESTA.
        self.services.wait_call(self.wait, True)
        self.services.call(self.siesta_port, 'step', timeStamp)

#  Prepare the output files for a super work flow. Need to remove any old output
#  files first before staging the plasma state.
        if os.path.exists(self.OUTPUT_FILES):
            os.remove(self.OUTPUT_FILES)
        self.services.stage_plasma_state()

#  The super flow may need to rename the output file. Check is the current state
#  matches if output file. If it does not rename the plasma state so it can be
#  staged.
        if not os.path.exists(self.OUTPUT_FILES):
            os.rename(self.current_siesta_state, self.OUTPUT_FILES)

#-------------------------------------------------------------------------------
#
#  SIESTA Driver finalize method. This cleans up afterwards.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'siesta_driver: finalize')

        self.wait = [
                     self.services.call_nonblocking(self.vmec_worker['init'], 'finalize', timeStamp),
                     self.services.call_nonblocking(self.vmec_worker['driver'], 'finalize', timeStamp),
                     self.services.call_nonblocking(self.siesta_port, 'finalize', timeStamp)
                    ]

        self.services.wait_call_list(self.wait, True)
