#! /usr/bin/env python

from ipsframework import Component

class rustbca_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeStamp=0.0):
        self.services.stage_input_files(self.INPUT_FILES)
        self.services.update_state()
        return

    def step(self, timeStamp=0.0):
        print('rustbca_driver: beginning step call') 
        try:
            worker_comp = self.services.get_port('WORKER')
        except Exception:
            self.services.exception('Error accessing worker component')
            raise
        self.services.call(worker_comp, 'step', 0.0)
        print('RustBCA Driver: finished worker call') 
        return

    def finalize(self, timeStamp=0.0):
        return

