#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for VMEC component. This wapper only takes a VMEC input file and
#  runs VMEC.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
import os
from omfit.classes.omfit_namelist import OMFITnamelist
from utilities import ZipState
from utilities import ScreenWriter
from utilities import NamelistItem

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
        
#  Stage state.
        self.services.stage_state()

#  Unzip files from the state. Use mode a so files can be read and written to.
        self.zip_ref = ZipState.ZipState(current_vmec_state, 'a')
        self.zip_ref.extract(self.current_vmec_namelist)

        if len(keywords) > 0:
            self.zip_ref.set_state(state='needs_update')
        
#  Update parameters in the namelist.
            namelist = OMFITnamelist(self.current_vmec_namelist,
                                     collect_arrays={
                                     'ns_array'    : {'default' : 0, 'shape' : (100,),    'offset' : (1,),     'sparray' : True},
                                     'niter_array' : {'default' : 0, 'shape' : (100,),    'offset' : (1,),     'sparray' : True},
                                     'rbs'         : {'default' : 0, 'shape' : (203,101), 'offset' : (-101,0), 'sparray' : True},
                                     'rbc'         : {'default' : 0, 'shape' : (203,101), 'offset' : (-101,0), 'sparray' : True},
                                     'zbs'         : {'default' : 0, 'shape' : (203,101), 'offset' : (-101,0), 'sparray' : True},
                                     'zbc'         : {'default' : 0, 'shape' : (203,101), 'offset' : (-101,0), 'sparray' : True},
                                     'am'          : {'default' : 0, 'shape' : (21,),     'offset' : (0,),     'sparray' : True},
                                     'ai'          : {'default' : 0, 'shape' : (21,),     'offset' : (0,),     'sparray' : True},
                                     'ac'          : {'default' : 0, 'shape' : (21,),     'offset' : (0,),     'sparray' : True},
                                     'am_aux_s'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                     'am_aux_f'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                     'ai_aux_s'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                     'ai_aux_f'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                     'ac_aux_s'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                     'ac_aux_f'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                     'raxis'       : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                     'zaxis'       : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                     'raxis_cc'    : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                     'raxis_cs'    : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                     'zaxis_cc'    : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                     'zaxis_cs'    : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                     'ftol_array'  : {'default' : 0, 'shape' : (100,),    'offset' : (1,),     'sparray' : True},
                                     'extcur'      : {'default' : 0, 'shape' : (300,),    'offset' : (1,),     'sparray' : True}
                                     })

            for key, value in keywords.items():
                NamelistItem.set(namelist['indata'], key, value)

            namelist.save()

#-------------------------------------------------------------------------------
#
#  VMEC Component step method. This runs vmec.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'vmec: step')

        flags = self.zip_ref.get_state()

        if 'state' in flags and flags['state'] == 'needs_update' or 'force_update' in keywords:
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

#  Add the wout file to the state.
            self.zip_ref.write([self.current_vmec_namelist, self.current_wout_file])

        else:
#  Update flags.
            self.zip_ref.set_state(state='unchanged')

        self.zip_ref.close()
            
        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  VMEC Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'vmec: finalize')
