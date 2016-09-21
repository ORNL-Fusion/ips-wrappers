#! /usr/bin/env python

import sys
import os
import subprocess
import getopt
import shutil
import math
from component import Component

class epaDriver(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0):
        return

    def step(self, timestamp=0):
        services = self.services

        #services.setWorkingDirectory(self)
        epaComp = services.get_port('EPA')
        fpComp = services.get_port('FOKKER_PLANCK')
        if(epaComp == None):
            print 'Error accessing Equilibrium component'
            sys.exit(1)
        if(fpComp == None):
            print 'Error accessing FP component'
            sys.exit(1)

        timeloop = services.get_time_loop()
        tlist_str = ['%.9f'%t for t in timeloop]
        print 'epa_fp_driver: tlist_str',tlist_str

        services.call(epaComp,'init')
        services.call(fpComp,'init')

        for t in tlist_str:
            print 'Current time = ', t
# call EPA and FP component
            services.call(epaComp,'step', t)
            services.call(fpComp,'step', t)

        services.call(epaComp, 'finalize')
        services.call(fpComp, 'finalize')

#      #services.setWorkingDirectory(self)

    def finalize(self, timestamp=0.0):
        pass
