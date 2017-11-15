#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for VMEC component. This driver only runs the VMEC component.
#
#-------------------------------------------------------------------------------

from component import Component

#-------------------------------------------------------------------------------
#
#  VMEC Driver Constructor
#
#-------------------------------------------------------------------------------
class vmec_driver(Component):
    def __init__(self, services, config):
        print('vmec_driver: Construct')
        Component.__init__(self, services, config)
            
#-------------------------------------------------------------------------------
#
#  VMEC Driver init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('vmec_driver: init')
        
        vmec_init = self.services.get_port('AINIT')
        self.services.call(vmec_init, 'init', timeStamp)

        self.vmec_comp = self.services.get_port('VMEC')
        self.services.call(self.vmec_comp, 'init', timeStamp)
    
#-------------------------------------------------------------------------------
#
#  VMEC Driver step method. This runs the vmec component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('vmec_driver: step')
        self.services.call(self.vmec_comp, 'step', timeStamp)

#-------------------------------------------------------------------------------
#
#  VMEC Driver finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('vmec_driver: finalize')
        self.services.call(self.vmec_comp, 'finalize', timeStamp)
