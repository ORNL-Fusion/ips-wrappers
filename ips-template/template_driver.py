
#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for TEMPLATE Driver component. This template shows an example of
#  how to set up an ips component driver.
#
#-------------------------------------------------------------------------------

from component import Component

#-------------------------------------------------------------------------------
#
#  TEMPLATE Driver Class
#
#-------------------------------------------------------------------------------
class template_driver(Component):

#-------------------------------------------------------------------------------
#
#  template_driver Component Constructor
#
#-------------------------------------------------------------------------------
    def __init__(self, services, config):
        print('template_driver: Construct')
        Component.__init__(self, services, config)
        self.running_components = {}
    
#-------------------------------------------------------------------------------
#
#  template_driver Component init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('template_driver: init')

#  Initialize the overall workflow. This is only getting called once so there is
#  no need to save the port.
        template_init_port = self.services.get_port('INIT')
        self.services.call(template_init_port, 'init', timeStamp)

#  Initialize the component for the first time.
        self.template = self.services.get_port('TEMPLATE')
        self.running_components['template:init'] = self.services.call_nonblocking(self.template, 'init', timeStamp)

#-------------------------------------------------------------------------------
#
#  template_driver Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('template_driver: step')

#  Wait until all the dependent components are finished.
        self.services.wait_call_list([self.running_components['template:init']], True)
        self.running_components['template:step'] = self.services.call_nonblocking(self.template, 'step', timeStamp)
#  Clear the waiting condition since we are no longer waiting for the init
#  method of template to run.
        del self.running_components['template:init']

#-------------------------------------------------------------------------------
#
#  template_driver Component finalize method.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('template_driver: finalize')
#  Wait until all the components are finished.
        self.services.wait_call_list(self.running_components.values(), True)
