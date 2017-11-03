#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for V3FIT component. This driver only uses the VMEC component to
#  stage in plasma state files. The v3fit components generate the actual wout
#  file.
#
#-------------------------------------------------------------------------------

from component import Component

#-------------------------------------------------------------------------------
#
#  V3FIT Driver Constructor
#
#-------------------------------------------------------------------------------
class v3fit_driver(Component):
    def __init__(self, services, config):
        print('v3fit_driver: Construct')
        Component.__init__(self, services, config)
            
#-------------------------------------------------------------------------------
#
#  V3FIT Driver init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('v3fit_driver: init')
        
        vmec_init = self.services.get_port('VMEC_INIT')
        self.services.call(vmec_init, 'init', timeStamp)

        v3fit_init = self.services.get_port('V3FIT_INIT')
        self.services.call(v3fit_init, 'init', timeStamp)

        self.v3fit_comp = self.services.get_port('V3FIT')
        self.services.call(self.v3fit_comp, 'init', timeStamp)
    
#-------------------------------------------------------------------------------
#
#  V3FIT Driver step method. This runs the vmec component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('v3fit_driver: step')
        self.services.call(self.v3fit_comp, 'step', timeStamp)

#-------------------------------------------------------------------------------
#
#  V3FIT Driver finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('v3fit_driver: finalize')
        self.services.call(self.v3fit_comp, 'finalize', timeStamp)
