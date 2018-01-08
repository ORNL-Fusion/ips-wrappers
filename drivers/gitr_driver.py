#! /usr/bin/env python

from  component import Component
import pickle

class gitr_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        try:
            self.worker_comp = self.services.get_port('WORKER')
        except Exception:
            self.services.exception('Error accessing worker component')
            raise
        self.services.call(self.worker_comp, 'init', 0.0)
        return

    def step(self, timeStamp=0.0):
        print 'gitr_driver: beginning step call' 
        self.services.call(self.worker_comp, 'step', 0.0)
        print 'GITRDriver: finished worker call' 
        return

    def finalize(self, timeStamp=0.0):
        return

