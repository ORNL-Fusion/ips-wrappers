#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for SIESTA component. This wapper only takes a SIESTA input file
#  and runs SIESTA.
#
#-------------------------------------------------------------------------------

from component import Component
import os
from omfit.classes.omfit_namelist import OMFITnamelist
from utilities import ZipState
from utilities import ScreenWriter
from utilities import NamelistItem

#-------------------------------------------------------------------------------
#
#  SIESTA Component Constructor
#
#-------------------------------------------------------------------------------
class siesta(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  SIESTA Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'siesta: init')
        self.services.stage_state()

        self.current_siesta_namelist = self.services.get_config_param('SIESTA_NAMELIST_INPUT')
        self.current_siesta_state = self.services.get_config_param('CURRENT_SIESTA_STATE')
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        current_wout_file = 'wout_{}.nc'.format(current_vmec_namelist.replace('input.','',1))
    
#  Stage state.
        self.services.stage_state()
    
#  Unzip files from the state. Use mode a so files can be read and written to.
        self.zip_ref = ZipState.ZipState(self.current_siesta_state, 'a')
        self.zip_ref.extract(self.current_siesta_namelist)
        self.zip_ref.extract(current_vmec_state)

        with ZipState.ZipState(current_vmec_state, 'r') as vmec_zip_ref:
            vmec_zip_ref.extract(current_wout_file)
            flags = vmec_zip_ref.get_state()
            if 'state' in flags and flags['state'] == 'updated':
                self.zip_ref.set_state(state='needs_update')

#  Update parameters in the namelist.
        self.set_namelist(wout_file=current_wout_file, **keywords)

#-------------------------------------------------------------------------------
#
#  SIESTA Component step method. This runs siesta.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'siesta: step')
        
        flags = self.zip_ref.get_state()
        
        if 'state' in flags and flags['state'] == 'needs_update':
            self.set_namelist(ladd_pert=True, lrestart=False)

            self.task_wait = self.services.launch_task(self.NPROC,
                                                       self.services.get_working_dir(),
                                                       self.SIESTA_EXE,
                                                       logfile = 'siesta1.log')

#  Update flags.
            self.zip_ref.set_state(state='updated')
            
#  Restart siesta to ensure convergence.
            self.services.wait_task(self.task_wait)
            
            self.set_namelist(ladd_pert=False, lrestart=True)
            self.task_wait = self.services.launch_task(self.NPROC,
                                                       self.services.get_working_dir(),
                                                       self.SIESTA_EXE,
                                                       logfile = 'siesta2.log')
                
#  Wait for SIESTA to finish.
            if (self.services.wait_task(self.task_wait) or not os.path.exists(self.restart_file)):
                self.services.error('siesta: step failed.')
                
#  Add the restart file to the state.
            self.zip_ref.write([self.current_siesta_namelist, self.restart_file])
        
        else:
#  Update flags.
            self.zip_ref.set_state(state='unchanged')
        
        self.zip_ref.close()
        self.services.update_state()
                
#-------------------------------------------------------------------------------
#
#  SIESTA Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'siesta: finalize')

#-------------------------------------------------------------------------------
#
#  SIESTA Component set_namelist method. This sets the namelist input file from
#  the keywords.
#
#-------------------------------------------------------------------------------
    def set_namelist(self, **keywords):
#  Update parameters in the namelist.
        namelist = OMFITnamelist(self.current_siesta_namelist,
                                 collect_arrays={
                                 'mres'     : {'default' : 0,   'shape' : (20,), 'offset' : (1,)},
                                 'HelPert'  : {'default' : 0.0, 'shape' : (20,), 'offset' : (1,)},
                                 'HelPertA' : {'default' : 0.0, 'shape' : (20,), 'offset' : (1,)}
                                 })
        
        self.restart_file = 'siesta_{}.nc'.format(namelist['siesta_info']['restart_ext'])
        
        if 'wout_file' in keywords:
            namelist['siesta_info']['wout_file'] = keywords['wout_file']
            self.update = True
            del keywords['wout_file']
        
        if len(keywords) > 0:
            self.zip_ref.set_state(state='needs_update')
        
            for key, value in keywords.iteritems():
                NamelistItem.set(namelist['siesta_info'], key, value)

        namelist.save()
