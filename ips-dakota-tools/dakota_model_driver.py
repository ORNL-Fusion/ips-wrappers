#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper to bridge the DAKOTA model driver. This defines a work flow that
#  can be used for the DAKOTA blackbox model.
#
#-------------------------------------------------------------------------------

from component import Component
import shutil
import ipsutil

#-------------------------------------------------------------------------------
#
#  DAKOTA model driver Component Constructor
#
#-------------------------------------------------------------------------------
class dakota_model_driver(Component):
    def __init__(self, services, config):
        print('dakota_model_driver: Construct')
        Component.__init__(self, services, config)
    
#-------------------------------------------------------------------------------
#
#  DAKOTA model driver Component init method. This method reads the
#  parameter file from datoka and adds to an existing namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('dakota_model_driver: init')

#  Create dummy files for plasma state.
        datoka_init = self.services.get_port('DAKOTA_INIT')
        dakota_init_id = self.services.call_nonblocking(datoka_init, 'init', timeStamp)

        vmec_init = self.services.get_port('VMEC_INIT')
        vmec_init_id = self.services.call_nonblocking(vmec_init, 'init', timeStamp)
        
        v3fit_init = self.services.get_port('V3FIT_INIT')
        v3fit_init_id = self.services.call_nonblocking(v3fit_init, 'init', timeStamp)

#  While waiting get all the other needed ports.
        dakota_to_vmec = self.services.get_port('DAKOTA_TO_VMEC')
        dakota_to_v3fit = self.services.get_port('DAKOTA_TO_V3FIT')
        self.v3fit_comp = self.services.get_port('V3FIT')
        self.v3fit_to_dakota_comp = self.services.get_port('V3FIT_TO_DAKOTA')

#  VMEC component init method needs the vmec_init call and dakota_init calls
#  completed before running.
        self.services.wait_call_list([dakota_init_id, vmec_init_id], True)
        self.dakota_to_vmec_id = self.services.call_nonblocking(dakota_to_vmec, 'init', timeStamp)

#  V3FIT component init method needs the v3fit_init call and dakota_init calls
#  completed before running.
#        self.services.wait_call_list([dakota_init_id, v3fit_init_id], True)
        self.services.wait_call_list([v3fit_init_id], True)

#  We don't need to wait on these calls until the step method is called. Store
#  the call id in the driver object.
        self.dakota_to_v3fit_id = self.services.call_nonblocking(dakota_to_v3fit, 'init', timeStamp)
    
#-------------------------------------------------------------------------------
#
#  DAKOTA model driver Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('dakota_model_driver: step')
    
# Run VMEC.
#        self.services.call(self.vmec_comp, 'step', timeStamp)

#  Wait until plasma state is ready.
        self.services.wait_call_list([self.dakota_to_vmec_id, self.dakota_to_v3fit_id], True)

#  All the rest need to run serially.

#  VMEC generated a new wout file. Call the V3FIT init method to state the
#  updated file. Bug in PARVMEC prevents using the wout file from vmec_comp.step
#  from being used in V3POST.
        self.services.call(self.v3fit_comp, 'init', timeStamp)
        self.services.call(self.v3fit_comp, 'step', timeStamp)

#  Create DAKOTA file with the cost functions computed from V3FIT.
        self.services.call(self.v3fit_to_dakota_comp, 'init', timeStamp)

#-------------------------------------------------------------------------------
#
#  DAKOTA model driver Component finalize method. This cleans up afterwards. Not
#  used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('dakota_model_driver: finalize')

        self.services.stage_plasma_state()
        
        current_dakota_cost_file = self.services.get_config_param('CURRENT_DAKOTA_COST_FILE')
        shutil.copyfile(current_dakota_cost_file, self.OUTPUT_FILES)

        root_directory = self.services.get_config_param('SIM_ROOT')
        current_directory = self.services.get_working_dir()
        ipsutil.copyFiles(current_directory, self.OUTPUT_FILES, root_directory)
