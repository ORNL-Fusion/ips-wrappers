#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive VMEC Runner Driver component.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from ips_component_utilities import ZipState
from ips_component_utilities import ScreenWriter
import os

#-------------------------------------------------------------------------------
#
#  Massive VMEC Runner Driver Component Constructor
#
#-------------------------------------------------------------------------------
class massive_serial_runner_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  Massive VMEC Runner Driver Compoenet init method. This method prepairs the
#  state.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_driver: init')
        
        if timeStamp == 0.0:
            self.current_state = self.services.get_config_param('CURRENT_MVR_STATE')
            self.massive_vmec_runner_port = self.services.get_port('MVR')
        
        self.wait = self.services.call_nonblocking(self.massive_vmec_runner_port,
                                                   'init', timeStamp)

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Driver Compoenet step method. This runs the vmec
#  component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_driver: step')

#  Run massive serial runner.
        self.services.wait_call(self.wait, True)
        self.services.call(self.massive_vmec_runner_port, 'step', timeStamp)
        
#  Prepare the output files for a super workflow. Need to remove any old output
#  files first before staging the state.
        if os.path.exists(self.OUTPUT_FILES):
            os.remove(self.OUTPUT_FILES)
        self.services.stage_state()

#  The super flow may need to rename the output file. Check if the current state
#  matches the output file. If it does not rename the state so it can be staged.
        if not os.path.exists(self.OUTPUT_FILES):
            os.rename(self.current_state)

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Driver Compoenet finalize method. This cleans up after.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_driver: finalize')
        self.services.call(self.massive_vmec_runner_port, 'finalize', timeStamp)
