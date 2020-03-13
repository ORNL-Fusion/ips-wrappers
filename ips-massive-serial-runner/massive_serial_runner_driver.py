#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive Serial Runner driver component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
import os

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component Constructor
#
#-------------------------------------------------------------------------------
class massive_serial_runner_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init method. This method prepairs the namelist input
#  file.
#
#------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_driver: init')

#  Initialize the massive serial runner.
        self.massive_serial_runner_port = self.services.get_port('MSR')
        self.wait = self.services.call_nonblocking(self.massive_serial_runner_port,
                                                   'init', timeStamp)

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner step method. This runs the vmec component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_driver: step'

#  Run massive serial runner.
        self.services.wait_call(self.wait, True)
        self.services.call(self.massive_serial_runner_port, 'step', timeStamp)

#  Prepare the output files for a super work flow. Need to remove any old output
#  files first before staging the state.
        if os.path.exists(self.OUTPUT_FILES):
            os.remove(self.OUTPUT_FILES)
        self.services.stage_state()

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner step method. This cleans up afterwards.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_driver: finalize')
        self.services.call(self.massive_serial_runner_port, 'finalize', timeStamp)
