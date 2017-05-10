#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper to bridge the V3FIT results to the DAKOTA input component.
#
#-------------------------------------------------------------------------------

from component import Component
from netCDF4 import Dataset
import numpy
import os

#-------------------------------------------------------------------------------
#
#  V3FIT to DAKOTA Component Constructor
#
#-------------------------------------------------------------------------------
class v3fit_to_ips_dakota(Component):
    def __init__(self, services, config):
        print('v3fit_to_ips_dakota: Construct')
        Component.__init__(self, services, config)
    
#-------------------------------------------------------------------------------
#
#  V3FIT to DAKOTA Component init method. This method reads the result file and
#  creates an input file to feed back into DAKOTA.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('v3fit_to_ips_dakota: init')
        
        self.services.stage_plasma_state()
    
        current_v3fit_result_file = self.services.get_config_param('CURRENT_V3FIT_RESULT_FILE')
    
#  Read the values from the result file to generate chi.
        v3fit_data_set = Dataset(current_v3fit_result_file)
    
        nsteps = v3fit_data_set.variables['nsteps'][:]
    
        observed = v3fit_data_set.variables['signal_observed_value'][:]
        modeled = v3fit_data_set.variables['signal_model_value'][:][nsteps - 1,:,0]
        weight = v3fit_data_set.variables['signal_weight'][:]
        sigma = v3fit_data_set.variables['signal_sigma'][:][nsteps - 1,:]

        v3fit_data_set.close()

        chis = numpy.sqrt(weight)*(observed - modeled)/sigma
    
#  Write chi to an dakota input file.
        current_dakota_cost_file = os.path.join(self.services.get_config_param('SIM_ROOT'), 'RESULT')
        dakota_in = open(current_dakota_cost_file, 'w')
    
        for chi in chis:
            dakota_in.write('{:22.15E}\n'.format(chi))

        dakota_in.close()

#-------------------------------------------------------------------------------
#
#  V3FIT to DAKOTA Component step method. Not used.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('v3fit_to_ips_dakota: step')

#-------------------------------------------------------------------------------
#
#  V3FIT to DAKOTA Component finalize method. This cleans up afterwards. Not
#  used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('v3fit_to_ips_dakota: finalize')
