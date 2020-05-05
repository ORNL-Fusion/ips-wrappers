#! /usr/bin/env python

from  component import Component
import fileinput
import hpic

class hpicWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeStamp=0.0):
        return

    def step(self, timeStamp=0.0):
        self.services.stage_plasma_state()
        
        print('Hello from hpic worker')
        
        dict = {'solpsTarg':'solpsTarg.txt',
                'solpsInds':[int(i) for i in self.SOLPS_INDS.split()],
                'hpicPaths':[str(i) for i in self.HPIC_PATHS.split()]}
        
        hpic.plot_hpic_iead(solps_path=dict['solpsTarg'], \
            HpicDataFolder = dict['hpicPaths'], \
            solpsIndex = dict['solpsInds'])        
        
        self.services.update_state()

        return
    
    def finalize(self, timeStamp=0.0):
        return
