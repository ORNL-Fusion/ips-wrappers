#! /usr/bin/env python

import sys
import os
import subprocess
from ipsframework import Component
import time
import shutil
from stat import *

class ng_genray(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        self.Nimrod_done = False
        self.last_timestamp = -1.0
        self.last_event = None

    def init(self, timestamp):
        # Figure out the working directory for GENRAY component.
        work_dir = self.services.get_working_dir()

        # Define PHYS_BIN_ROOT for the GENRAY component.
        phys_bin_root = self.services.get_config_param('PHYS_BIN_ROOT')

        # Specify the GENRAY binary path.
        genray_bin = os.path.join(phys_bin_root, 'genray', 'bin', 'xgenray')

        # Specify the prepare_genray_input and process_genray_output paths.
        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_genray_input')
        process_output_bin = os.path.join(self.BIN_PATH, 'process_genray_output')

        # Subscribe GENRAY to the PHYS_NIMSET_EVENTS channel so it can receive the
        # announcement of new EQDSK files.
        self.services.subscribe('PHYS_NIMSET_EVENTS', "process_nimset_events")
        newComment='ng_genray subscribed to the PHYS_NIMSET_EVENTS channel.'
	self.services.info(newComment)
        self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)

        # Now wait for the announcement of a new EQDSK on PHYS_NIMSET_EVENTS.
        self.eqdsk_arrived = False
        while True:
            # Check to see if anything happened.
            self.services.process_events()          
            if self.eqdsk_arrived:
                # New EQDSK available.  Copy it to GENRAY's work directory, and
		# break out of the check-to-see-if-anything-happened loop.
                genray_file = self.last_event['EQDSK_FILE']
                target_file = os.path.join(work_dir, 'g000001.0000')
                self.services.info('Copying nimset output to GENRAY working directory.')
                shutil.copy(genray_file, target_file)
                break
            # Otherwise, nothing new happened; wait two seconds and check again.
            time.sleep(2)
	    
        # If we get here, it means that we have a new EQDSK.
        newComment='ng_genray received a new EQDSK file path from the PHYS_NIMSET_EVENTS channel.'
        self.services.info(newComment)
        self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)

        # Move the current plasma state files to the GENRAY working directory.
        newComment='ng_genray is getting the current Plasma State data.'
        self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
        self.services.info(newComment)
        self.services.stage_state()

        # Move the input files (presently only the genray.template file) to the GENRAY working
        # directory so they will be available for use.
        newComment='ng_genray is moving the template GENRAY input file to its working directory.'
        self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
        self.services.stage_input_files(self.INPUT_FILES)

        # Now call prepare_genray_input to make the GENRAY input file with the
        # correct NIMROD profiles.
        newComment='ng_genray is calling prepare_genray_input.'
        self.services.send_portal_event(event_comment=newComment)
        try:
            self.services.info(newComment)
            subprocess.call([prepare_input_bin,'step',self.RFMODE, \
                            self.ISOURCE_STRING,self.GENRAYNML,self.ADJ_READ, \
                            self.PS_ADD_NML])
            self.services.info("Finished prepare_genray_input.")
        except:
            self.services.exception('Error invoking %s.', prepare_input_bin)
        self.services.send_portal_event('GENRAY_EVENT',event_comment='Finished prepare_genray_input.')

        # Now call GENRAY itself.
        newComment='ng_genray is launching GENRAY.'
        self.services.info(newComment)
        self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
        task_id = self.services.launch_task(self.NPROC, work_dir, genray_bin, \
                                      logfile='genray_init.log')
        ret_code = self.services.wait_task(task_id)
        newComment='The initial GENRAY call finished executing.'
        self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
        self.services.info(newComment)

        # Provenance
        shutil.copy('genray.nc','genray.nc.0')
        shutil.copy('g000001.0000','g000001.0000.0')

        # Now call process_genray_output to merge the results with the Plasma State.
        newComment='ng_genray is calling process_genray_output.'
        self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
        try:
            self.services.info(newComment)
            subprocess.call([process_output_bin,self.RFMODE,self.ISOURCE_STRING])
            self.services.info("Finished process_genray_output.")
        except:
            self.services.exception('Error invoking %s.', process_output_bin)
        self.services.send_portal_event('GENRAY_EVENT',event_comment='Finished process_genray_output.')

        # Now move the GENRAY output files to the results directory.
        newComment='ng_genray is storing output files from GENRAY.'
        self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
        self.services.stage_output_files(timestamp,'genray.nc.0')
        self.services.info("Moved genray.nc file to results directory.")

        # Merge plasma state.
        shutil.copy('cur_state.cdf','cur_state.cdf.0')
        self.services.send_portal_event('GENRAY_EVENT',event_comment='ng_genray is updating the Plasma State.')
        self.services.update_state() 

        # Now publish the new genray.nc filepath to PHYS_GENRAY_EVENTS.
        event_data={}
        event_data['GENRAY_FILE'] = os.path.join(work_dir, 'genray.nc')
        self.services.publish('PHYS_GENRAY_EVENTS', 'GENRAY_FILE_AVAILABLE', event_data)
        self.services.info("Published file path for new genray.nc to PHYS_GENRAY_EVENTS.")

        newComment='The INIT call for the RF component (ng_genray) completed successfully.'
        self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
        self.services.info(newComment)
        return

    def step(self, timestamp):
        self.services.info("Starting the ng_genray step call.")

        # Get the working directory.
        work_dir = self.services.get_working_dir()
	self.services.info("Got work_dir for ng_genray step call.")

        # Figure out expected file sizes so we don't open partially written files.
        base_eqdsk_path=os.path.join(work_dir, 'g000001.0000.0')
	eqdsk_size=os.path.getsize(base_eqdsk_path)
	base_ncfile_path=os.path.join(work_dir, 'genray.nc.0')
	ncfile_size=os.path.getsize(base_ncfile_path)

        # Get PHYS_BIN_ROOT.
        phys_bin_root = self.services.get_config_param('PHYS_BIN_ROOT')
	self.services.info("Got PHYS_BIN_ROOT for ng_genray step call.")

        # Get binary path.
        genray_bin = os.path.join(phys_bin_root, 'genray', 'bin', 'xgenray')
	self.services.info("Got GENRAY binary path for ng_genray step call.")

        # Get prepare_genray_input and process_genray_output binary paths.
        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_genray_input')
	self.services.info("Got prepare_genray_input binary path.")
        process_output_bin = os.path.join(self.BIN_PATH, 'process_genray_output')
        self.services.info("Got process_genray_output binary path.")

        start = self.last_timestamp
        stop = timestamp

        self.Nimrod_done = False
        counter = 1
        while not self.Nimrod_done:
            # Perpetually check to see if new EQDSK files are available.
            while True:
                # Default assumption - no new EQDSK.
                self.eqdsk_arrived = False
                # Now check to see if the default assumption is wrong.
                self.services.process_events()
                if self.eqdsk_arrived:
                    # New EQDSK.  Copy it to GENRAY's work directory and break out of the
		    # continuously-check-for-new-EQDSKs loop.
                    genray_file = self.last_event['EQDSK_FILE']
                    target_file = os.path.join(work_dir, 'g000001.0000')
                    self.services.info('Copying nimset output to GENRAY working directory.')
                    shutil.copy(genray_file, target_file)
                    break
                # Otherwise, no new EQDSK; check to see if NIMROD is done.
                if self.Nimrod_done:
                    # Then NIMROD has finished and we won't need to process any more EQDSK
                    # files, so break out of the continuously-check-for-new-EQDSKs loop.
                    break
		# Otherwise, we're still running and can expect more EQDSKs momentarily.
		# Wait two seconds and then try again.
		time.sleep(2)
		
            # If we get here, it is because NIMROD is done or we have a new EQDSK to process.
	    # If NIMROD is done, break out and end the STEP call.
            if self.Nimrod_done:
                break

            # Otherwise, we have a new file to process.
            newComment='ng_genray received a new EQDSK file path from the PHYS_NIMSET_EVENTS channel.'
            self.services.info(newComment)
            self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)

            # Make sure we are done writing it.
	    neweqdsk_size=os.path.getsize(target_file)
	    doneWriting = False
	    while not doneWriting:
	        if (neweqdsk_size<eqdsk_size):
	            time.sleep(.1)
		    neweqdsk_size=os.path.getsize(target_file)
		else:
		    doneWriting = True
		    break
		    
            # Move the current plasma state files to the GENRAY working directory.
            newComment='ng_genray is getting the current Plasma State data.'
            self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
            self.services.info(newComment)
            self.services.stage_state()

            # Now call prepare_genray_input to make the GENRAY input file with the
            # correct NIMROD profiles.
            newComment='ng_genray is calling prepare_genray_input.'
            self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
            try:
                self.services.info(newComment)
                subprocess.call([prepare_input_bin,'step',self.RFMODE, \
                                self.ISOURCE_STRING,self.GENRAYNML,self.ADJ_READ, \
                                self.PS_ADD_NML])
                self.services.info("Finished prepare_genray_input.")
            except:
                self.services.exception('Error invoking %s.', prepare_input_bin)
            self.services.send_portal_event('GENRAY_EVENT',event_comment='Finished prepare_genray_input.')

            newComment='ng_genray is launching GENRAY.'
            self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
            # Now call GENRAY itself to make the genray.nc file that NIMROD reads.
            task_id = self.services.launch_task(self.NPROC, work_dir, genray_bin, \
                                           logfile='genray_step.log'+str(counter))
            ret_code = self.services.wait_task(task_id)
            self.services.send_portal_event('GENRAY_EVENT',event_comment='Finished GENRAY.')
 
            # Now call process_genray_output to merge the results with the Plasma State.
            newComment='ng_genray is calling process_genray_output.'
            self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
            try:
                self.services.info(newComment)
                subprocess.call([process_output_bin,self.RFMODE, \
                                self.ISOURCE_STRING])
                self.services.info("Finished process_genray_output.")
            except:
                self.services.exception('Error invoking %s.', process_output_bin)
            self.services.send_portal_event('GENRAY_EVENT',event_comment='Finished process_genray_output.')

            # Provenance of various files.
            shutil.copy('g000001.0000','g000001.0000.'+str(counter))
            shutil.copy('genray.nc','genray.nc.'+str(counter))
        
            # Now move the GENRAY output files to the results directory.
            newComment="ng_genray is storing GENRAY's output file genray.nc."
            self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
            self.services.info(newComment)
            self.services.stage_output_files(timestamp,'genray.nc.'+str(counter))

            # Update the Plasma State.
	    newComment='ng_genray is updating the Plasma State.'
            self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
            shutil.copy('cur_state.cdf','cur_state.cdf.'+str(counter))
            self.services.update_state()

            # Publish the announcement of a new Plasma State to the monitor channel so 
	    # that the monitor file can be updated.
            new_state_file = 'cur_state.cdf.'+str(counter)
            new_state_file_path = os.path.join(work_dir, new_state_file)
            event_data = {}
            event_data['EVENT_TYPE'] = 'UPDATE_MONITOR'
            event_data['STATE_FILE'] = new_state_file_path
            event_data['TIME_STAMP'] = counter
            self.services.publish('SIMULATION_MONITOR', 'UPDATE_MONITOR', event_data)
            newComment='ng_genray requested that the monitor component publish an update.'
            self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)

            # Update timestamp.  This is a dumb way to do it.  Works for now though.
            # Eventually, should use NIMROD dumpfile tag, I guess.
            self.services.update_time_stamp(new_time_stamp = counter)

            # Increment counter.
            counter = counter + 1

            # Make sure we are done writing genray.nc.
	    newncfile_path=os.path.join(work_dir,'genray.nc')
	    newncfile_size=os.path.getsize(newncfile_path)
	    doneWriting = False
	    while not doneWriting:
	        if (newncfile_size<ncfile_size):
	            time.sleep(.1)
		    newncfile_size=os.path.getsize(newncfile_path)
		else:
		    doneWriting = True
		    break

            # Publish the availability of a new genray.nc to PHYS_GENRAY_EVENTS.
            event_data = {}
            event_data['GENRAY_FILE'] = os.path.join(work_dir, 'genray.nc')
            self.services.publish('PHYS_GENRAY_EVENTS', 'GENRAY_FILE_AVAILABLE', event_data)
            newComment='ng_genray published an updated genray.nc file path to PHYS_GENRAY_EVENTS.'
            self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)

            # Now check to see if NIMROD has finished while we were working on all that.
            self.services.process_events()
            # I don't think we need these, since if self.Nimrod_done becomes true, we just kick
	    # out of the loop anyway.
	    # if self.Nimrod_done:
	        # break

        # NIMROD is done, so signal the monitor component that we're done too.
        event_data = {}
        event_data['EVENT_TYPE'] = 'END_SIMULATION'
        self.services.publish('SIMULATION_MONITOR', 'UPDATE_MONITOR', event_data)
	
        newComment='The STEP call for the RF component (ng_genray) completed successfully.'
        self.services.send_portal_event('GENRAY_EVENT',event_comment=newComment)
        self.services.info(newComment)
        return

    def process_nimset_events(self, topicName, theEvent):
        self.last_event = theEvent.getBody()
        if self.last_event.has_key('EQDSK_FILE'):
            # A new eqdsk file is available.
            self.eqdsk_arrived = True
        if self.last_event.has_key('NIMROD_DONE'):
            # NIMROD is finished executing.
            self.Nimrod_done = True
            self.eqdsk_arrived = False
            return
        else:
            pass
        return

    def finalize(self, timestamp):
        pass
