#! /usr/bin/env python

import sys
import os
import subprocess
from component import Component
import time
import shutil
from stat import *

class Genray_xmhd(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        self.done = False
        self.last_timestamp = -1.0
        self.event_arrived = False
        self.last_event = None
        self.last_dump_time = None

    def init(self, timestamp):
        self.services.subscribe('PHYS_NIMROD_EVENTS', "process_nimrod_events")
	self.services.info("Subscribed GENRAY to NIMROD events.")
        self.services.stage_input_files(self.INPUT_FILES)
        self.services.info("Put GENRAY input files in working directory for init step.")

#       eqdsk_dir = self.services.get_config_param('EQDSKSTART')
# KLUGE - this is hardcoded.  Need to figure out how to make it work in general.
# The above lines are one attempt at it but I don't have the dynamic Python 
# variables figured out very well yet.
        work_dir = self.services.get_working_dir()
        eqdsk_file = os.path.join(work_dir, '../MHD__Nimrod_2/fake_eqdsk')
	self.services.info("Got path for eqdsk file from NIMROD.")

        phys_bin_root = self.services.get_config_param('PHYS_BIN_ROOT')
	self.services.info("Got PHYS_BIN_ROOT for GENRAY init call.")
#        genray_bin = os.path.join(phys_bin_root, 'genray', 'bin', 'xgenray')
# KLUGE - this is a serial version with my modifications.  Need a parallel version
# that has my modifications in PHYS_BIN_ROOT.  Doesn't need to be the same as RWH's.
        genray_bin = os.path.join(phys_bin_root, 'genray', 'bin', 'xgenray_tj')
	self.services.info("Got GENRAY binary path for init call.")

        target_file = os.path.join(work_dir, 'fake_eqdsk')
        self.services.info("Defined path for getting fake_eqdsk file for GENRAY init call.")
	
        shutil.copy(eqdsk_file, target_file)
        self.services.debug('Received new eqdsk file: %s', eqdsk_file)

        self.services.info("Calling GENRAY for initial eqdsk processing.")
        task_id = self.services.launch_task(self.NPROC, work_dir, genray_bin,
                                       logfile='genray_step.log')
        ret_code = self.services.wait_task(task_id)
        genray_fname = 'genray.nc'

# Signal the nimrod component to compute new eqdsk file
        event_data = {}
        event_data['GENRAY_FILE'] = os.path.join(work_dir, genray_fname)
        self.services.publish('PHYS_GENRAY_EVENTS', 'GENRAY_FILE_AVAILABLE', event_data)
	self.services.info("Published GENRAY event and finished GENRAY init.")
        return

    def step(self, timestamp):
        self.services.info('Starting step() - timestamp = %s', str(timestamp))
        work_dir = self.services.get_working_dir()
	print "Got work_dir for GENRAY step call"
        phys_bin_root = self.services.get_config_param('PHYS_BIN_ROOT')
	print "Got PHYS_BIN_ROOT for GENRAY step call"
#KLUGE see also sim-nimrod-genray.conf line 60
#        genray_bin = os.path.join(phys_bin_root, 'genray', 'bin', 'xgenray_par')
#        genray_bin = os.path.join(phys_bin_root, 'genray', 'bin', 'xgenray')
        genray_bin = os.path.join(phys_bin_root, 'genray', 'bin', 'xgenray_tj')
	print "Got GENRAY binary path for step call"
        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_genray_input')
	print "Got GENRAY's prepare_genray_input binary path"
        process_output_bin = os.path.join(self.BIN_PATH, 'process_genray_output')
        print "Got GENRAY's process_genray_output binary path"
	
        self.services.stage_plasma_state()

        start = self.last_timestamp
        stop = timestamp

        counter = 0
        while True:
            # Wait for nimrod to generate an eqdsk file
            self.event_arrived = False
            while True:
                self.services.process_events()
                if self.event_arrived:
                    break
                time.sleep(2)

            if self.done:
                break
            #
            # Copy new sources to work directory (using the same name every time)
            eqdsk_file = self.last_event['EQDSK_FILE']
            target_file = os.path.join(work_dir, 'fake_eqdsk')
            shutil.copy(eqdsk_file, target_file)
            self.services.debug('Received new eqdsk file: %s', eqdsk_file)
#            try:
#                print "About to call prepare_genray_input"
#		subprocess.call([prepare_input_bin,'step','EC','0',
#		                'genray.dat','disabled','disabled'])
#            except:
#                self.services.exception('Error invoking %s', prepare_input_bin)

            print "Calling GENRAY, I think"
            task_id = self.services.launch_task(self.NPROC, work_dir, genray_bin,
                                           logfile='genray_step.log')
            ret_code = self.services.wait_task(task_id)
            # PLace holder for the actual genray invocation
            genray_fname = 'genray.nc'
#            genray_fname = 'genray.nc.' + str(counter)
#            counter += 1
#            subprocess.call(['touch', genray_fname])

#            try:
#                print "About to call process_genray_output"
#		subprocess.call([process_output_bin])
#            except Exception:
#                self.services.exception('Error invoking %s', process_output_bin)
            #
            # Signal the nimrod component to compute new eqdsk file
            #
            event_data = {}
            event_data['GENRAY_FILE'] = os.path.join(work_dir, genray_fname)
            self.services.publish('PHYS_GENRAY_EVENTS', 'GENRAY_FILE_AVAILABLE', event_data)

        return

    def process_nimrod_events(self, topicName, theEvent):
        self.event_arrived = True
        self.last_event = theEvent.getBody()
        try:
            self.done =  self.last_event['NIMROD_DONE']
        except KeyError:
            pass
        return

    def finalize(self, timestamp):
        pass
