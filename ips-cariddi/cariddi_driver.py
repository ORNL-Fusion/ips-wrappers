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
import json

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
            self.current_siesta_state = self.services.get_config_param('CURRENT_SIESTA_STATE')
            self.current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
            self.current_vmec_profile = self.services.get_config_param('CURRENT_VMEC_PROFILE')
            current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
            self.current_wout_file = 'wout_{}.nc'.format(current_vmec_namelist.replace('input.','',1))

            self.time_steps = int(self.services.get_config_param('NUMBER_OF_TIME_STEPS'))

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

            eq_keywords = {'vmec__mgrid_file' : self.services.get_config_param('MGRID_FILE')}

#  Initialize the equilibrium.
            self.services.call(self.eq_worker['init'], 'init', timeStamp)

            self.cariddi_port = self.services.get_port('CARIDDI')

#-------------------------------------------------------------------------------
#
#  CARIDDI Driver step method. This runs the V3FIT subworkflow and CARIDDI
#  component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi_driver: step')

        outter_loop_time_stamp = 0
        inner_loop_time_stamp = 0
        time_stamp = '{}.{}'.format(outter_loop_time_stamp,
                                    inner_loop_time_stamp)
        ScreenWriter.screen_output(self, 'quiet', 'Initialization: Time Stamp = {}'.format(time_stamp))

#  Zero out the mgrid_file and get intial equilibrium profiles for t = 0.0.
        self.services.call(self.cariddi_port, 'init', time_stamp)
        self.services.call(self.cariddi_port, 'step', time_stamp,
                           task_list=['zero_mgrid',
                                      'get_profile'])
        self.get_updated_state()

#  Extract out the updated profiles. Since this is the first ever run, set the
#  the mgrid file path to point to the zeroed out file.
        eq_keywords = self.get_eq_profiles()
        eq_keywords['vmec__mgrid_file'] = self.services.get_config_param('MGRID_FILE')

#  Run the inital equilibrium with the first set of profiles.
        self.services.call(self.eq_worker['driver'], 'init', float(time_stamp), **eq_keywords)
        self.services.call(self.eq_worker['driver'], 'step', float(time_stamp), force_update=True)
        self.get_updated_substate(time_stamp)

#  Get the surface currents and generate then eddy.nc file.
        inner_loop_time_stamp = inner_loop_time_stamp + 1
        time_stamp = '{}.{}'.format(outter_loop_time_stamp,
                                    inner_loop_time_stamp)
        ScreenWriter.screen_output(self, 'quiet', 'Initialization: Time Stamp = {}'.format(time_stamp))

        self.services.call(self.cariddi_port, 'init', time_stamp)
        self.services.call(self.cariddi_port, 'step', time_stamp,
                           task_list=['get_current',
                                      'make_eddy'])

#  Start outer loop. FIXME: Need a time counter for the outter loop.
#--  Outer Loop  ---------------------------------------------------------------
        while outter_loop_time_stamp < self.time_steps - 1:
#  Increment outer loop counter and reset inner loop counter.
            outter_loop_time_stamp = outter_loop_time_stamp + 1
            inner_loop_time_stamp = 0
            time_stamp = '{}.{}'.format(outter_loop_time_stamp,
                                        inner_loop_time_stamp)
            ScreenWriter.screen_output(self, 'quiet', 'Outer Loop: Time Stamp = {}'.format(time_stamp))

#  Get the new profiles.
            self.services.call(self.cariddi_port, 'init', time_stamp)
            self.services.call(self.cariddi_port, 'step', time_stamp,
                               task_list=['get_profile'])
            self.get_updated_state()

#  Update VMEC for t = t + 1 profiles.
            eq_keywords = self.get_eq_profiles()
            self.services.call(self.eq_worker['driver'], 'init', float(time_stamp), **eq_keywords)
            self.services.call(self.eq_worker['driver'], 'step', float(time_stamp), force_update=True)
            self.get_updated_substate(time_stamp)

#  Calculate the surface current from the new wout file and update the mgrid
#  file.
            self.services.call(self.cariddi_port, 'init', time_stamp)
            wait_call = self.services.call_nonblocking(self.cariddi_port, 'step', time_stamp,
                                                       task_list=['get_current',
                                                                  'set_mgrid'])

#  While the surface current is being computed, extract the current woutfile to
#  get the stopping criteria.
            magnetic_axis = self.get_magnetic_axis()

#  Wait for the mgrid file to be finished updating.
            self.services.wait_call(wait_call, True)

            if outter_loop_time_stamp == 2.0:
                return

#  FIXME: Choose some stopping criteria for the loop...
            delta_magnetic_axis = 100

#--  Inner Loop  ---------------------------------------------------------------
            while delta_magnetic_axis > 1.0E-7:
                inner_loop_time_stamp = inner_loop_time_stamp + 1
                time_stamp = '{}.{}'.format(outter_loop_time_stamp,
                                            inner_loop_time_stamp)
                ScreenWriter.screen_output(self, 'quiet',
                                           'Inner Loop: Time Stamp = {}, Magnetic Axis = {}, Delta Axis = {}'.format(time_stamp,
                                                                                                                     magnetic_axis,
                                                                                                                     delta_magnetic_axis))

#  Update equilibrium with new vacuum contribution.
                self.services.call(self.eq_worker['driver'], 'init', float(time_stamp))
                self.services.call(self.eq_worker['driver'], 'step', float(time_stamp), force_update=True)
                self.get_updated_substate(time_stamp)

#  Calculate the surface current from the new wout file and update the mgrid
#  file.
                self.services.call(self.cariddi_port, 'init', time_stamp)
                wait_call = self.services.call_nonblocking(self.cariddi_port, 'step', time_stamp,
                                                           task_list=['get_current',
                                                                      'set_mgrid'])

#  Measure change in magnetic axis position and update magnetic axis.
                new_magnetic_axis = self.get_magnetic_axis()
                delta_magnetic_axis = abs(new_magnetic_axis - magnetic_axis)
                magnetic_axis = new_magnetic_axis

#  Wait for the mgrid file to be finished updating.
                self.services.wait_call(wait_call, True)

            ScreenWriter.screen_output(self, 'quiet',
                                       'Inner Loop End: Time Stamp = {}, Magnetic Axis = {}, Delta Axis = {}'.format(time_stamp,
                                                                                                                     magnetic_axis,
                                                                                                                     delta_magnetic_axis))

#--  End Inner Loop  -----------------------------------------------------------

#  Equilibrium converged to new position. Update eddy.nc file.
            self.services.call(self.cariddi_port, 'init', time_stamp)
            self.services.call(self.cariddi_port, 'step', time_stamp,
                               task_list=['make_eddy'])

        ScreenWriter.screen_output(self, 'quiet', 'Outer Loop End: Time Stamp = {}'.format(time_stamp))

#--  End Outer Loop  -----------------------------------------------------------

#-------------------------------------------------------------------------------
#
#  CARIDDI Driver finalize method. This cleans up afterwards.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi_driver: finalize')

        wait_call = [
                     self.services.call_nonblocking(self.eq_worker['init'], 'finalize', timeStamp),
                     self.services.call_nonblocking(self.eq_worker['driver'], 'finalize', timeStamp),
                     self.services.call_nonblocking(self.cariddi_port, 'finalize', timeStamp)
                    ]

        self.services.wait_call_list(wait_call, True)

#-------------------------------------------------------------------------------
#
#  CARIDDI Driver get magnetic axis. Note this is not an IPS method.
#
#-------------------------------------------------------------------------------
    def get_magnetic_axis(self):
        with ZipState.ZipState(self.current_v3fit_state, 'r') as v3fit_ref:
            if self.current_siesta_state in v3fit_ref:
                v3fit_ref.extract(self.current_siesta_state)
                with ZipState.ZipState(self.current_siesta_state, 'r') as siesta_ref:
                    siesta_ref.extract(self.current_vmec_state)
            else:
                v3fit_ref.extract(self.current_vmec_state)

        with ZipState.ZipState(self.current_vmec_state, 'r') as vmec_ref:
            vmec_ref.extract(self.current_wout_file)

        return netCDF4.Dataset(self.current_wout_file).variables['rmnc'][0,0]

#-------------------------------------------------------------------------------
#
#  CARIDDI Driver get equilibrium profiles. Note this is not an IPS method.
#
#-------------------------------------------------------------------------------
    def get_eq_profiles(self):
        flags = self.zip_ref.get_state()

        if 'state' in flags and flags['state'] != 'unchanged':
            self.zip_ref.extract(self.current_vmec_profile)
            with open(self.current_vmec_profile, 'r') as json_ref:
                return json.load(json_ref)
        else:
            return {}

#-------------------------------------------------------------------------------
#
#  CARIDDI Driver Get the updated state. The current open state must be closed
#  before the updated state can be staged. State the current state for output.
#  Note this is not an IPS method.
#
#-------------------------------------------------------------------------------
    def get_updated_state(self):
        self.zip_ref.close()
        self.services.stage_state()
        self.zip_ref = ZipState.ZipState(self.current_cariddi_state, 'a')

#-------------------------------------------------------------------------------
#
#  CARIDDI Driver Get the updated substate. After the substate is added to the
#  current state update the current state so the other components can get the
#  current substate. The open state must be closed before the state can be
#  updated. Reopen the current state after update.
#
#-------------------------------------------------------------------------------
    def get_updated_substate(self, timeStamp=0.0):
        self.services.stage_subflow_output_files()

        with ZipState.ZipState(self.current_v3fit_state, 'r') as v3fit_ref:
            flags = v3fit_ref.get_state()

        if 'state' in flags and flags['state'] != 'unchanged':
            self.zip_ref.write(self.current_v3fit_state)
            self.zip_ref.set_state(state='updated')
            self.zip_ref.close()
            self.services.update_state()
            self.services.stage_output_files(timeStamp, self.current_cariddi_state)
            self.zip_ref = ZipState.ZipState(self.current_cariddi_state, 'a')
