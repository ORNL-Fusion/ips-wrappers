#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for V3POST component. This the VMEC component then use the
#  generated wout file to run v3post.
#
#-------------------------------------------------------------------------------

from component import Component

#-------------------------------------------------------------------------------
#
#  V3POST Driver Constructor
#
#-------------------------------------------------------------------------------
class v3post_driver(Component):
    def __init__(self, services, config):
        print('v3post_driver: Construct')
        Component.__init__(self, services, config)
            
#-------------------------------------------------------------------------------
#
#  V3POST Driver init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('v3post_driver: init')
        
        vmec_init = self.services.get_port('VMEC_INIT')
        self.services.call(vmec_init, 'init', timeStamp)

        v3fit_init = self.services.get_port('V3FIT_INIT')
        self.services.call(v3fit_init, 'init', timeStamp)

        self.vmec_comp = self.services.get_port('VMEC')
        self.services.call(self.vmec_comp, 'init', timeStamp)

        self.v3fit_comp = self.services.get_port('V3FIT')
    
#-------------------------------------------------------------------------------
#
#  V3POST Driver step method. This runs the vmec component followed by v3fit.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('v3post_driver: step')
        
        self.services.call(self.vmec_comp, 'step', timeStamp)
        
#  VMEC generated a new wout file. Call the V3FIT init method to state the
#  updated file.
        self.services.call(self.v3fit_comp, 'init', timeStamp)
        self.services.call(self.v3fit_comp, 'step', timeStamp)

#-------------------------------------------------------------------------------
#
#  V3POST Driver finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('v3post_driver: finalize')
        self.services.call(self.vmec_comp, 'finalize', timeStamp)
        self.services.call(self.v3fit_comp, 'finalize', timeStamp)
