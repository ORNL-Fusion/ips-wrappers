#! /usr/bin/env python

from ipsframework import Component

class rustbca_worker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeStamp=0.0):
        return

    def step(self, timeStamp=0.0):
        # Here the worker will update the input.toml
        # of RustBCA at every time instance
        self.services.stage_state()
        print('Running rustBCA worker step...')
        task_id = self.services.launch_task(
            self.NPROC, 
            self.services.get_working_dir(), 
            self.EXECUTABLE,
            self.INPUT_FILES,
            )
        self.services.wait_task(task_id)
        print('done!')
        return
    
    def finalize(self, timeStamp=0.0):
        return
    
