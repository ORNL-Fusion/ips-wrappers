#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 SOLPS5 omponent 
 -----------------------------------------------------------------------
"""

import sys,os,shutil
import subprocess
from numpy import *

from component import Component

class solps5(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):

        services = self.services
        services.stage_plasma_state()
        pwd = services.get_working_dir()

        # Change to SOLPS5 dir

        # Change to csh

        # Source env.solps5.edison

            # Perhaps ...
            # http://stackoverflow.com/questions/3503719/emulating-bash-source-in-python
            # http://stackoverflow.com/questions/16744334/python-sourcing-a-csh-and-passing-setenv-to-a-new-subprocess

        # Create link to solps run directory

        # Change to run directory

        # Stage input files

        return

    def step(self, timeStamp=0.0):

        services = self.services
        SetupScript = self.env
        
        #--- entry

        #--- excutable

        #--- stage plasma state files

        #--- get plasma state file names

        #--- stage input files


        #--- generate genray input

        #--- run genray

        #--- get genray output

        #--- update plasma state files

        #--- archive output files

        return


    def restart(self, timeStamp=0.0):

        return
    

    def checkpoint(self, timestamp=0.0):

        return

    def finalize(self, timeStamp=0.0):

        return
    
