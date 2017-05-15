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
        self.dakota_init_id = self.services.call_nonblocking(datoka_init, 'init', timeStamp)

        vmec_init = self.services.get_port('VMEC_INIT')
        self.vmec_init_id = self.services.call_nonblocking(vmec_init, 'init', timeStamp)

        siesta_init = self.services.get_port('SIESTA_INIT')
        self.siesta_init_id = self.services.call_nonblocking(siesta_init, 'init', timeStamp)
        
        v3fit_init = self.services.get_port('V3FIT_INIT')
        self.v3fit_init_id = self.services.call_nonblocking(v3fit_init, 'init', timeStamp)

#  While waiting get all the other needed ports.
        self.vmec_comp = self.services.get_port('VMEC')
        self.siesta_comp = self.services.get_port('SIESTA')
        self.v3fit_comp = self.services.get_port('V3FIT')
        self.v3fit_to_dakota_comp = self.services.get_port('V3FIT_TO_DAKOTA')

        self.vmec_params = {}
        self.siesta_params = {}
        self.v3fit_params = {}
        
        for key, value in self.__dict__.iteritems():
            if ('dakota_vmec_' in key):
                self.vmec_params[key.replace('dakota_vmec_', '')] = value
            elif ('dakota_siesta_' in key):
                self.siesta_params[key.replace('dakota_siesta_', '')] = value
            elif ('dakota_v3fit_' in key):
                self.v3fit_params[key.replace('dakota_v3fit_', '')] = value

#-------------------------------------------------------------------------------
#
#  DAKOTA model driver Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('dakota_model_driver: step')
    
#  Run VMEC.
        self.services.wait_call_list([self.dakota_init_id, self.vmec_init_id], True)
        self.services.call(self.vmec_comp, 'init', timeStamp, name_list_params=self.vmec_params)
        self.services.call(self.vmec_comp, 'step', timeStamp)
        
#  Run the first pass of SIESTA.
        self.services.wait_call_list([self.siesta_init_id], True)
        self.siesta_params['LADD_PERT'] = 'T'
        self.siesta_params['LRESISTIVE'] = 'T'
        self.siesta_params['LRESTART'] = 'F'
        self.services.call(self.siesta_comp, 'init', timeStamp, name_list_params=self.siesta_params)
        self.services.call(self.siesta_comp, 'step', timeStamp)

#  Run the second pass of SIESTA.
        self.siesta_params['LADD_PERT'] = 'F'
        self.siesta_params['LRESISTIVE'] = 'F'
        self.siesta_params['LRESTART'] = 'T'
        self.services.call(self.siesta_comp, 'init', timeStamp, name_list_params=self.siesta_params)
        self.services.call(self.siesta_comp, 'step', timeStamp)
        
#  Run V3FIT.
        self.services.wait_call_list([self.v3fit_init_id], True)
        self.services.call(self.v3fit_comp, 'init', timeStamp, name_list_params=self.v3fit_params)
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
