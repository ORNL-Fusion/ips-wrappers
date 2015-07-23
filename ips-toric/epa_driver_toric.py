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
        if(epaComp == None):
            print 'Error accessing Equilibrium component'
            sys.exit(1)
        rfComp = services.get_port('RF_IC')
        if(rfComp == None):
            print 'Error accessing RF component'
            sys.exit(1)

# Get timeloop for simulation
        timeloop = services.get_time_loop()
        tlist_str = ['%.2f'%t for t in timeloop]

# Call init for each component
        services.call(epaComp,'init')
        services.call(rfComp,'init')

        for t in tlist_str:
            print 'Current time = ', t
# call EPA component
            services.call(epaComp,'step', t)
            services.call(rfComp,'step', t)

        services.call(epaComp, 'finalize')
        services.call(rfComp, 'finalize')

        #services.setWorkingDirectory(self)

    def finalize(self, timestamp=0.0):
        pass
