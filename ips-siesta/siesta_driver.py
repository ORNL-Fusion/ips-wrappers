#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for SIESTA component. This driver only runs the SIESTA component.
#
#-------------------------------------------------------------------------------

from component import Component

#-------------------------------------------------------------------------------
#
#  SIESTA Driver Constructor
#
#-------------------------------------------------------------------------------
class siesta_driver(Component):
    def __init__(self, services, config):
        print('siesta_driver: Construct')
        Component.__init__(self, services, config)
            
#-------------------------------------------------------------------------------
#
#  SIESTA Driver init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('siesta_driver: init')
        
        vmec_init = self.services.get_port('VMEC_INIT')
        self.services.call(vmec_init, 'init', timeStamp)
        
        siesta_init = self.services.get_port('SIESTA_INIT')
        self.services.call(siesta_init, 'init', timeStamp)
        
        self.vmec_comp = self.services.get_port('VMEC')
        self.services.call(self.vmec_comp, 'init', timeStamp)
    
#-------------------------------------------------------------------------------
#
#  SIESTA Driver step method. This runs the vmec component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('siesta_driver: step')
        
        self.services.call(self.vmec_comp, 'step', timeStamp)
        
        self.siesta_comp = self.services.get_port('SIESTA')
        self.services.call(self.siesta_comp, 'init', timeStamp)
        self.services.call(self.siesta_comp, 'step', timeStamp)
    
#-------------------------------------------------------------------------------
#
#  SIESTA Driver finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('siesta_driver: finalize')
        self.services.call(self.vmec_comp, 'finalize', timeStamp)
        self.services.call(self.siesta_comp, 'finalize', timeStamp)
