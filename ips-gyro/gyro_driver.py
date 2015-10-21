#! /usr/bin/env python

from  component import Component

class gyro_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        return

    def step(self, timeStamp=0.0):
        print 'gyro_driver: beginning step call' 
        try:
            gyro_comp = self.services.get_port('GYROKINETIC')
        except Exception:
            self.services.exception('Error accessing GYROKINETIC component')
            raise

        self.services.call(gyro_comp, 'init', 0.0)
        print 'gyro_driver: finished GYROKINETIC init call' 

        self.services.call(gyro_comp, 'step', 0.0)
        print 'gyro_driver: finished GYROKINETIC step call' 
        return

        self.services.call(gyro_comp, 'finalize', 0.0)
        print 'gyro_driver: finished GYROKINETIC finalize call' 
        return

    def finalize(self, timeStamp=0.0):
        return

