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

#-------------------------------------------------------------------------------
#
#  SIESTA Component Constructor
#
#-------------------------------------------------------------------------------
class siesta(Component):
    def __init__(self, services, config):
        print('siesta: Construct')
        Component.__init__(self, services, config)
        self.task_queue = {}

#-------------------------------------------------------------------------------
#
#  SIESTA Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        print('siesta: init')
        self.services.stage_plasma_state()

        self.current_siesta_namelist = self.services.get_config_param('SIESTA_NAMELIST_INPUT')
        self.current_siesta_state = self.services.get_config_param('CURRENT_SIESTA_STATE')
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        current_wout_file = 'wout_{}.nc'.format(current_vmec_namelist.replace('input.','',1))
    
#  Stage plasma state.
        self.services.stage_plasma_state()
    
#  Unzip files from the plasma state. Use mode a so files can be read and
#  written to.
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
        print('siesta: step')
        
        flags = self.zip_ref.get_state()
        
        if 'state' in flags and flags['state'] == 'needs_update':
            self.set_namelist(ladd_pert=True, lrestart=False)

            self.task_queue['siesta'] = self.services.launch_task(self.NPROC,
                                                                  self.services.get_working_dir(),
                                                                  self.SIESTA_EXE,
                                                                  logfile = 'siesta1.log')

#  Update flags.
            self.zip_ref.set_state(state='updated')
            
#  Restart siesta to ensure convergence.
            self.services.wait_task(self.task_queue['siesta'], True)
            del self.task_queue['siesta']
            self.set_namelist(ladd_pert=False, lrestart=True)
            self.task_queue['siesta'] = self.services.launch_task(self.NPROC,
                                                                  self.services.get_working_dir(),
                                                                  self.SIESTA_EXE,
                                                                  logfile = 'siesta2.log')
                
#  Wait for SIESTA to finish.
            if (self.services.wait_task(self.task_queue['siesta'], True) or not os.path.exists(self.restart_file)):
                self.services.error('siesta: step failed.')
            del self.task_queue['siesta']
                
#  Add the restart file to the plasma state.
            self.zip_ref.write([self.current_siesta_namelist, self.restart_file])
            self.zip_ref.close()
    
            self.services.update_plasma_state()
        
        else:
#  Update flags.
            self.zip_ref.set_state(state='unchanged')
            self.zip_ref.close()
            
            self.services.update_plasma_state()
                
#-------------------------------------------------------------------------------
#
#  SIESTA Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('siesta: finalize')
        self.services.wait_tasklist(self.task_queue.values(), True)

#-------------------------------------------------------------------------------
#
#  SIESTA Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def set_namelist(self, **keywords):
#  Update parameters in the namelist.
        namelist = OMFITnamelist(self.current_siesta_namelist)
        
        self.restart_file = 'siesta_{}.nc'.format(namelist['siesta_info']['restart_ext'])
        
        if 'wout_file' in keywords and namelist['siesta_info']['wout_file'] != keywords['wout_file']:
            namelist['siesta_info']['wout_file'] = keywords['wout_file']
            self.update = True
            del keywords['wout_file']
        
        if len(keywords) > 0:
            self.zip_ref.set_state(state='needs_update')
        
            for key, value in keywords.iteritems():
                if '(' in key :
                    key, indices = key.split('(')
                    indices, extra = indices.split(')')
                    indices = [[int(i) - 1] for i in indices.split(',')]
                    namelist['siesta_info'][key][indices] = value
                else:
                    namelist['siesta_info'][key] = value

        namelist.save()
