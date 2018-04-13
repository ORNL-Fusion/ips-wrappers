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
        self.async_queue = {}
        self.ports = {}

#-------------------------------------------------------------------------------
#
#  VMEC Driver init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        print('vmec_driver: init')
        
        self.ports['vmec'] = self.services.get_port('VMEC')
        self.async_queue['vmec:init'] = self.services.call_nonblocking(self.ports['vmec'], 'init',
                                                                       timeStamp, **keywords)
    
#-------------------------------------------------------------------------------
#
#  VMEC Driver step method. This runs the vmec component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('vmec_driver: step')
        
        self.services.wait_call(self.async_queue['vmec:init'], True)
        self.async_queue['vmec:step'] = self.services.call_nonblocking(self.ports['vmec'],
                                                                       'step', timeStamp)
        del self.async_queue['vmec:init']
        
        self.services.wait_call(self.async_queue['vmec:step'], True)
        del self.async_queue['vmec:step']
    
#-------------------------------------------------------------------------------
#
#  VMEC Driver finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('vmec_driver: finalize')
        
        self.services.wait_call_list(self.async_queue.values(), True)
        self.async_queue = {}
        
        for portname, port in self.ports.iteritems():
            self.async_queue[portname + ':finalize'] = self.services.call_nonblocking(port, 'finalize',
                                                                                      timeStamp)
        
        self.services.wait_call_list(self.async_queue.values(), True)
        self.async_queue = {}

