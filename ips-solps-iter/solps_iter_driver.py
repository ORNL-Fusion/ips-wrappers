#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for SOLPS-ITER Driver component. This wapper runs b2-eierne then
#  computes synthetic signals from the output.
#
#-------------------------------------------------------------------------------

from component import Component

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Driver Component Constructor
#
#-------------------------------------------------------------------------------
class solps_iter_driver(Component):
    def __init__(self, services, config):
        print('solps_iter_driver: Construct')
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Driver Component init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('solps_iter_driver: init')
        
        solps_init = self.services.get_port('SOLPS_INIT')
        self.services.call(solps_init, 'init', timeStamp)

        self.solps = self.services.get_port('SOLPS')
        self.services.call(self.solps, 'init', timeStamp)
    
#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Driver step Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('solps_iter_driver: step')

        self.services.call(self.solps, 'step', timeStamp, task=self.SOLPS_TASK)

#-------------------------------------------------------------------------------
#
#  SOLPS-ITER Driver Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('solps_iter_driver: finalize')
