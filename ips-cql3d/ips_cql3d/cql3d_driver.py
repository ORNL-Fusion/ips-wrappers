#! /usr/bin/env python

# ssf - some comments explaining each step would be helpful for users
#       I'll add what I *think* is going on.  Please correct me if I am wrong.

import sys
import os
import subprocess
import getopt
import shutil
import math
from ipsframework import Component

class cql3dDriver(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0):
        # ssf - insert additional init stuff here
        return

    def step(self, timestamp=0):
        services = self.services

        # ssf - set working directory and cd into it
        #services.setWorkingDirectory(self)

        # ssf - get references to the components to run in simulation
        fpComp = services.get_port('FOKKER_PLANCK')

        if(fpComp == None):
            print 'Error accessing FOKKER_PLANCK physics components'
            sys.exit(1)

        # ssf - get timeloop for simulation
        timeloop = services.get_time_loop()
        tlist_str = ['%.2f'%t for t in timeloop]

        # ssf - call init for each component
        services.call(fpComp,'init')

        # ssf - iterate through the timeloop
        for t in tlist_str:
            print 'Current time = ', t

            # ssf - call step for each component
            services.call(fpComp,'step', t)

            # ssf - post step processing: stage plasma state, stage output
            services.stage_state()
            services.stage_output_files(t, self.OUTPUT_FILES)

        # ssf - post simulation: call finalize on each component
        services.call(fpComp, 'finalize')

    def finalize(self, timestamp=0.0):
        # ssf - insert additional post simulation processing here
        pass
