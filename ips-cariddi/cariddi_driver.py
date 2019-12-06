#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for CARIDDI component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
import os
import shutil
import netCDF4

#-------------------------------------------------------------------------------
#
#  CARIDDI Driver Constructor
#
#-------------------------------------------------------------------------------
class cariddi_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        self.eq_worker = {'sim_name': None, 'init': None, 'driver': None}

#-------------------------------------------------------------------------------
#
#  CARIDDI Driver init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi_driver: init')

#  Get config filenames.
        self.current_cariddi_state = self.services.get_config_param('CURRENT_CARIDDI_STATE')

#  We need to pass the inputs to the V3FIT child workflow.
        self.services.stage_state()

        self.zip_ref = ZipState.ZipState(self.current_cariddi_state, 'a')

#  If this is the first call, set up the V3FIT sub workflow.
        if timeStamp == 0.0:
            self.current_v3fit_state = self.services.get_config_param('CURRENT_V3FIT_STATE')
            self.current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
            current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
            self.current_wout_file = 'wout_{}.nc'.format(current_vmec_namelist.replace('input.','',1))
            
            if os.path.exists('eq_input_dir'):
                shutil.rmtree('eq_input_dir')
            os.mkdir('eq_input_dir')

            self.zip_ref.extract(self.current_v3fit_state)
            shutil.copy2(self.current_v3fit_state, 'eq_input_dir')

            self.cariddi_port = self.services.get_port('CARIDDI')

#  Get keys for the V3FIT sub workflow.
            keys = {'PWD'                  : self.services.get_config_param('PWD'),
                    'USER_INPUT_FILES'     : self.current_v3fit_state,
                    'SIM_NAME'             : '{}_v3fit'.format(self.services.get_config_param('SIM_NAME')),
                    'LOG_FILE'             : 'log.{}_v3fit.warning'.format(self.services.get_config_param('SIM_NAME')),
                    'OUTPUT_LEVEL'         : self.services.get_config_param('OUTPUT_LEVEL'),
                    'VMEC_NAMELIST_INPUT'  : self.services.get_config_param('VMEC_NAMELIST_INPUT'),
                    'V3FIT_NAMELIST_INPUT' : self.services.get_config_param('V3FIT_NAMELIST_INPUT')
                   }

            v3fit_config = self.services.get_config_param('V3FIT_CONFIG')

            (self.eq_worker['sim_name'],
             self.eq_worker['init'],
             self.eq_worker['driver']) = self.services.create_sub_workflow('v3fit', v3fit_config,
                                                                           keys, 'eq_input_dir')

            eq_keywords = {'vmec__MGRID_FILE' : self.services.get_config_param('MGRID_FILE')}

#  Initialize the equilibrium.
            self.services.call(self.eq_worker['init'], 'init', timeStamp)
            self.services.call(self.eq_worker['driver'], 'init', timeStamp, **eq_keywords)

            self.cariddi_port = self.services.get_port('CARIDDI')

#-------------------------------------------------------------------------------
#
#  CARIDDI Driver step method. This runs the V3FIT subworkflow and CARIDDI
#  component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi_driver: step')

#  FIXME: Choose some stopping criteria for the loop..
        magnetic_axis = 0
        delta_magnetic_axis = 100

        while delta_magnetic_axis > 1.0E-6:
            self.services.call(self.cariddi_port, 'init', timeStamp)
            self.services.call(self.cariddi_port, 'step', timeStamp)

#  Run the equilibrium. Replace values in the V3FIT state. Since the changes are
#  external to the namelist input, force an update to the V3FIT state.
            if timeStamp > 0.0:
                self.services.call(self.eq_worker['driver'], 'init', timeStamp)

            self.services.call(self.eq_worker['driver'], 'step', timeStamp, force_update=True)
            if timeStamp > 0.0:
                with ZipState.ZipState(self.current_v3fit_state, 'r') as v3fit_ref:
                    v3fit_ref.extract(self.current_vmec_state)
                    with ZipState.ZipState(self.current_vmec_state, 'r') as vmec_ref:
                        vmec_ref.extract(self.current_wout_file)

                        data = netCDF4.Dataset(self.current_wout_file)
                        new_magnetic_axis = data.variables['rmnc'][1,1]
                        delta_magnetic_axis = abs(new_magnetic_axis - magnetic_axis)
                        magnetic_axis = new_magnetic_axis

#  After the equilibrium has run update the state.
            self.services.stage_subflow_output_files()
            self.zip_ref.write(self.current_v3fit_state)
            self.zip_ref.close()
            self.services.update_state()

#  Reopen the zip file after the state has been updated.
            self.zip_ref = ZipState.ZipState(self.current_cariddi_state, 'a')

            self.services.stage_output_files(timeStamp, self.current_cariddi_state)

            timeStamp += 1.0

#-------------------------------------------------------------------------------
#
#  CARIDDI Driver finalize method. This cleans up afterwards.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi_driver: finalize')

        self.wait = [
                     self.services.call_nonblocking(self.eq_worker['init'], 'finalize', timeStamp),
                     self.services.call_nonblocking(self.eq_worker['driver'], 'finalize', timeStamp),
                     self.services.call_nonblocking(self.cariddi_port, 'finalize', timeStamp)
                    ]

        self.services.wait_call_list(self.wait, True)
