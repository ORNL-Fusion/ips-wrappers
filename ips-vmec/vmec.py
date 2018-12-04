#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC component. This wapper only takes a VMEC input file and
#  runs VMEC.
#
#-------------------------------------------------------------------------------

from component import Component
import os
from omfit.classes.omfit_namelist import OMFITnamelist
from utilities import ZipState
from utilities import ScreenWriter

#-------------------------------------------------------------------------------
#
#  VMEC Component Constructor
#
#-------------------------------------------------------------------------------
class vmec(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  VMEC Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'vmec: init')
        
        self.current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        self.current_wout_file = 'wout_{}.nc'.format(self.current_vmec_namelist.replace('input.','',1))
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        
#  Stage plasma state.
        self.services.stage_plasma_state()

#  Unzip files from the plasma state. Use mode a so files can be read and
#  written to.
        self.zip_ref = ZipState.ZipState(current_vmec_state, 'a')
        self.zip_ref.extract(self.current_vmec_namelist)

        if len(keywords) > 0:
            self.zip_ref.set_state(state='needs_update')
        
#  Update parameters in the namelist.
            namelist = OMFITnamelist(self.current_vmec_namelist)

            for key, value in keywords.iteritems():
                namelist['indata'][key] = value
    
            namelist.save()

#-------------------------------------------------------------------------------
#
#  VMEC Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'vmec: step')

        flags = self.zip_ref.get_state()

        if 'state' in flags and flags['state'] == 'needs_update':
            task_wait = self.services.launch_task(self.NPROC,
                                                  self.services.get_working_dir(),
                                                  self.VMEC_EXE,
                                                  self.current_vmec_namelist,
                                                  logfile = 'vmec.log')

#  Update flags.
            self.zip_ref.set_state(state='updated')

#  Wait for VMEC to finish.
            if (self.services.wait_task(task_wait) and not os.path.exists(self.current_wout_file)):
                self.services.error('vmec: step failed.')

#  Add the wout file to the plasma state.
            self.zip_ref.write([self.current_vmec_namelist, self.current_wout_file])

        else:
#  Update flags.
            self.zip_ref.set_state(state='unchanged')

        self.zip_ref.close()
            
        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  VMEC Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'vmec: finalize')
