#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for VMEC component. This driver only runs the VMEC component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
import os

#-------------------------------------------------------------------------------
#
#  VMEC Driver Constructor
#
#-------------------------------------------------------------------------------
class vmec_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  VMEC Driver init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'vmec_driver: init')

#  Separate out the vmec keywords.
        vmec_keywords = {}
        for key, value in keywords.items():
            if 'vmec__' in key:
                vmec_keywords[key.replace('vmec__','',1)] = value
    
#  Initialize vmec.
        self.vmec_port = self.services.get_port('VMEC')
        self.wait = self.services.call_nonblocking(self.vmec_port, 'init',
                                                   timeStamp, **vmec_keywords)
    
#-------------------------------------------------------------------------------
#
#  VMEC Driver step method. This runs the vmec component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'vmec_driver: step')

#  Run vmec.
        self.services.wait_call(self.wait, True)
        self.services.call(self.vmec_port, 'step', timeStamp, **keywords)
    
#  Prepare the output files for a super work flow. Need to remove any old output
#  files first before staging the state.
        if os.path.exists(self.OUTPUT_FILES):
            os.remove(self.OUTPUT_FILES)
        self.services.stage_state()

#  The super flow may need to rename the output file. Check is the current state
#  matches if output file. If it does not rename the state so it can be staged.
        if not os.path.exists(self.OUTPUT_FILES):
            os.rename(self.services.get_config_param('CURRENT_VMEC_STATE'),
                      self.OUTPUT_FILES)
    
#-------------------------------------------------------------------------------
#
#  VMEC Driver finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'vmec_driver: finalize')
        self.services.call(self.vmec_port, 'finalize', timeStamp)
