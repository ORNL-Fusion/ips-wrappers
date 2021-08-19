#! /usr/bin/env python

import sys
import os
import subprocess
import getopt
import shutil
import string
from ipsframework import Component

class cql3d(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        self.firstTime  = True
        self.curTime = -1.0
        self.prevTime = -1.0
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0):
        print 'cql3d.init() called'
        self.services.stage_input_files(self.INPUT_FILES)
        self.OUTPUT_FILES += ' '+ self.services.get_config_param('RUN_ID')+'.ps'
        return

    def step(self, timeStamp):
        print 'cql3d.step() calld'

        if (self.services == None) :
            print 'Error in aorsa:;step () : init() function not called before step().'
            sys.exit(1)
        services = self.services

        prepare_input = os.path.join(self.BIN_PATH, 'prepare_cql3d_input')
        process_output  = os.path.join(self.BIN_PATH, 'process_cql3d_output')
        cql3d_bin = os.path.join(self.BIN_PATH, 'xcql3d')

# Copy current and prior state over to working directory
        services.stage_state()

# Call prepare_cql3d-input
        cur_state_file = services.get_config_param('CURRENT_STATE')
        prior_state_file = services.get_config_param('PRIOR_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

# Copy files (Rename) to meet prepare_cql3d_input expectations
        fileList = [(cur_state_file,'cur_state.cdf'),
                    (prior_state_file, 'prior_state'),
                    (cur_eqdsk_file, 'eqdskin')]
        for (src, dist) in fileList:
            if (src != dist):
                print 'Copying %s to %s' %(src, dist)
                try:
                    shutil.copyfile(src, dist)
                except Exception, e:
                    print 'fp_cql3d.py 1 error: Error copying state files' , e
#              return 1
        if (self.firstTime):
            inputFileName = 'cqlinput.1r.1'
            self.firstTime = False
            self.curTime = float(timeStamp)
            self.prevTime = self.curTime
            nsteps=1
#         deltaT = self.curTime
            deltaT = 1.e-9
        else:
            inputFileName = 'cqlinput.1r.2'
            self.prevTime = self.curTime
            self.curTime = float(timeStamp)
            nsteps=50
            deltaT = (self.curTime - self.prevTime)/nsteps

        shutil.copyfile(inputFileName, 'cqlinput')
        cmd = ([prepare_input, '%i' % (nsteps), '%e' % (deltaT), services.get_config_param('RUN_ID')])
        print cmd
        retcode = subprocess.call(cmd)
        if (retcode != 0):
            print 'fp_cql3d.py error: Error executing ', prepare_input
            return 1

#BH: Copy old cqlinput to cqlinput.0.  Copy cqlinput_new to cqlinput
#      shutil.copyfile(os.path.join(workdir,inputFileName), os.path.join(workdir, 'cqlinput.0'))
        shutil.copyfile('cqlinput_new', 'cqlinput')

# Call cql3d (in parallel)
#      retcode = services.launchJob(cql3d_bin, self.NPROC)
        retcode = subprocess.call([cql3d_bin])
        if (retcode != 0):
            print 'Error executing command: ', cql3d_bin
            sys.exit(1)

        fd = open('distrfunc', 'r')
        outBuf = ''
        buf = fd.readline()
        while (buf != ''):
            if ' RFFILE  =' in buf:
                fd.readline()
                fd.readline()
            else:
                outBuf += buf
            buf = fd.readline()
        fd.close()

        fd = open('distrfunc', 'w')
        fd.write(outBuf)
        fd.close()

        fname = services.get_config_param('RUN_ID')+'.nc'
        shutil.copy(fname, 'mnemonic.nc')

# Call process_output and copy files over
        retcode = subprocess.call([process_output])
        if (retcode != 0):
            print 'Error executing',  process_output
            sys.exit(1)

# Update (original) plasma state
        for (src, dist) in fileList:
            if (src != dist):
                print 'Copying %s to %s' %(dist, src)
                try:
                    shutil.copyfile(dist, src)
                except Exception, e:
                    print 'fp_cql3d.py 1 error: Error copying state files' , e
#              return 1

        services.update_state()

# "Archive" output files in history directory
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    def finalize(self, timestamp=0.0):
        print 'cql3d.finalize() called'
