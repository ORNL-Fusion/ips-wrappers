#! /usr/bin/env python

from  component import Component

class hello_worker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        return

    def step(self, timeStamp=0.0):
        print 'Hello from hello_worker'
        return
    
    def finalize(self, timeStamp=0.0):
        return
    
