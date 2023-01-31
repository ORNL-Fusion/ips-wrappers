#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for V3FIT component. This wapper only takes a V3FIT input file
#  and runs V3FIT.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from ips_component_utilities import ZipState
from ips_component_utilities import ScreenWriter
import os

#-------------------------------------------------------------------------------
#
#  CARIDDI Component Constructor
#
#-------------------------------------------------------------------------------
class cariddi(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  CARIDDI Component init method. This method prepairs the input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi: init')

#  Stage state.
        self.services.stage_state()

#  Get config filenames.
        if timeStamp == '0.0':
            self.current_cariddi_matrix = self.services.get_config_param('CARIDDI_MATRIX_FILE')
            self.current_cariddi_geometry = self.services.get_config_param('CARIDDI_GEOMETRY_FILE')
            self.cariddi_matrix_path = self.services.get_config_param('CARIDDI_MATRIX_PATH')
            self.current_cariddi_state = self.services.get_config_param('CURRENT_CARIDDI_STATE')

            self.current_v3fit_state = self.services.get_config_param('CURRENT_V3FIT_STATE')

            self.current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
            current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
            self.current_wout_file = 'wout_{}.nc'.format(current_vmec_namelist.replace('input.','',1))

            self.mgrid_file = self.services.get_config_param('MGRID_FILE')
            self.current_vmec_profile = self.services.get_config_param('CURRENT_VMEC_PROFILE')

            self.zip_ref = ZipState.ZipState(self.current_cariddi_state, 'a')
            if 'A1.nc' in self.zip_ref:
                self.zip_ref.extract('A1.nc')
        else:
            self.zip_ref = ZipState.ZipState(self.current_cariddi_state, 'a')

# Extract input files.
        self.flags = self.zip_ref.get_state()

        if 'state' in self.flags and self.flags['state'] != 'unchanged':
            self.zip_ref.extract(self.current_v3fit_state)

            with ZipState.ZipState(self.current_v3fit_state, 'r') as v3fit_zip_ref:
                v3fit_zip_ref.extract(self.current_vmec_state)

            with ZipState.ZipState(self.current_vmec_state, 'r') as vmec_zip_ref:
                if self.current_wout_file in vmec_zip_ref:
                    vmec_zip_ref.extract(self.current_wout_file)

#-------------------------------------------------------------------------------
#
#  CARIDDI Component step method. This runs V3FIT.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi: step')

        if isinstance(timeStamp, str):
            time_index = timeStamp.split('.')[0]
        else:
            time_index = int(timeStamp)

        if 'task_list' in keywords:
            for task in keywords['task_list']:
#  Initalize eddy file.
                if task == 'init_eddy':
                    ScreenWriter.screen_output(self, 'quiet', 'Initalize eddy current file: Time Stamp = {}'.format(timeStamp))
                    task_wait = self.services.launch_task(self.NPROC,
                                                          self.services.get_working_dir(),
                                                          self.CARIDDI_PRE_EXE,
                                                          '--matrix_file={}'.format(self.current_cariddi_matrix),
                                                          '--matrix_path={}'.format(self.cariddi_matrix_path),
                                                          logfile = 'cariddi_pre_init_eddy_{}.log'.format(timeStamp))
                    self.services.wait_task(task_wait)
                    continue

#  Initalize vector potential file.
                if task == 'init_potential':
                    ScreenWriter.screen_output(self, 'quiet', 'Initalize potential file: Time Stamp = {}'.format(timeStamp))
                    task_wait = self.services.launch_task(self.NPROC,
                                                          self.services.get_working_dir(),
                                                          self.CARIDDI_PRE_EXE,
                                                          '--coordinate_file={}'.format(self.current_cariddi_geometry),
                                                          '--matrix_path={}'.format(self.cariddi_matrix_path),
                                                          logfile = 'cariddi_pre_init_potential_{}.log'.format(timeStamp))
                    self.services.wait_task(task_wait)
                    self.zip_ref.write('A1.nc')
                    self.zip_ref.set_state(state='updated')
                    continue

#  Zero out the inital mgrid fields.
                if task == 'zero_mgrid':
                    ScreenWriter.screen_output(self, 'quiet', 'Zero out mgrid file: Time Stamp = {}'.format(timeStamp))
                    task_wait = self.services.launch_task(self.NPROC,
                                                          self.services.get_working_dir(),
                                                          self.CARIDDI_BIN_EXE,
                                                          '--mgrid_file={}'.format(self.mgrid_file),
                                                          '--zero_mgrid',
                                                          logfile = 'cariddi_bin_zero_mgrid_{}.log'.format(timeStamp))
                    self.services.wait_task(task_wait)
                    continue

#  Get new profiles for vmec for current time step.
                if task == 'get_profile':
                    ScreenWriter.screen_output(self, 'quiet', 'Get equilibrium profiles: Time Stamp = {}'.format(timeStamp))
                    task_wait = self.services.launch_task(self.NPROC,
                                                          self.services.get_working_dir(),
                                                          self.CARIDDI_PRE_EXE,
                                                          '--matrix_path={}'.format(self.cariddi_matrix_path),
                                                          '--time_index={}'.format(time_index),
                                                          '--vmec_input={}'.format(self.current_vmec_profile),
                                                          logfile = 'cariddi_pre_get_profile_{}.log'.format(timeStamp))
                    self.services.wait_task(task_wait)
                    self.zip_ref.write(self.current_vmec_profile)
                    self.zip_ref.set_state(state='updated')
                    continue

#  Update eddy current file.
                if task == 'make_eddy':
                    ScreenWriter.screen_output(self, 'quiet', 'Update eddy.nc file: Time Stamp = {}'.format(timeStamp))
                    task_wait = self.services.launch_task(self.NPROC,
                                                          self.services.get_working_dir(),
                                                          self.CARIDDI_PRE_EXE,
                                                          '--matrix_path={}'.format(self.cariddi_matrix_path),
                                                          '--time_index={}'.format(time_index),
                                                          '--vmec_current=A1.nc',
                                                          logfile = 'cariddi_pre_make_eddy_{}.log'.format(timeStamp))
                    self.services.wait_task(task_wait)
                    continue

#  Save fields.
                if task == 'save_fields':
                    ScreenWriter.screen_output(self, 'quiet', 'Save fields file: Time Stamp = {}'.format(timeStamp))
                    task_wait = self.services.launch_task(self.NPROC,
                                                          self.services.get_working_dir(),
                                                          self.CARIDDI_PRE_EXE,
                                                          '--matrix_path={}'.format(self.cariddi_matrix_path),
                                                          '--mgrid_file={}'.format(self.mgrid_file),
                                                          logfile = 'cariddi_pre_save_fields_{}.log'.format(timeStamp))
                    self.services.wait_task(task_wait)
                    continue

#  Compute virtual current at coupling surface.
                if task == 'get_current':
                    ScreenWriter.screen_output(self, 'quiet', 'Getting surface current: Time Stamp = {}'.format(timeStamp))
                    task_wait = self.services.launch_task(self.NPROC,
                                                          self.services.get_working_dir(),
                                                          self.SURFACE_EXE,
                                                          '-woutf={}'.format(self.current_wout_file),
                                                          '-surff=A1.nc',
                                                          '-para=-1',
                                                          logfile = 'surface_get_current_{}.log'.format(timeStamp))
                    self.services.wait_task(task_wait)
                    self.zip_ref.write('A1.nc')
                    self.zip_ref.set_state(state='updated')
                    continue

#  Set vacuum fields for the eddy current.
                if task == 'set_mgrid':
                    ScreenWriter.screen_output(self, 'quiet', 'Setting Vacuum Eddy Current Fields: Time Stamp = {}'.format(timeStamp))
                    task_wait = self.services.launch_task(self.NPROC,
                                                          self.services.get_working_dir(),
                                                          self.CARIDDI_BIN_EXE,
                                                          '--mgrid_file={}'.format(self.mgrid_file),
                                                          '--matrix_path={}'.format(self.cariddi_matrix_path),
                                                          '--vmec_current=A1.nc',
                                                          '--set_fields',
                                                          logfile = 'cariddi_bin_set_mgrid_{}.log'.format(timeStamp))
                    self.services.wait_task(task_wait)
                    continue

                if task == 'scale_fields':
                    ScreenWriter.screen_output(self, 'quiet', 'Scaling Vacuum Eddy Current Fields: Time Stamp = {}'.format(timeStamp))
                    task_wait = self.services.launch_task(self.NPROC,
                                                          self.services.get_working_dir(),
                                                          self.CARIDDI_BIN_EXE,
                                                          '--mgrid_file={}'.format(self.mgrid_file),
                                                          '--matrix_path={}'.format(self.cariddi_matrix_path),
                                                          '--vmec_current=A1.nc',
                                                          '--scale_fields',
                                                          logfile = 'cariddi_bin_scale_fields_{}.log'.format(timeStamp))
                    self.services.wait_task(task_wait)
                    continue
        else:
            self.zip_ref.set_state(state='unchanged')

        self.zip_ref.close()
        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  CARIDDI Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi: finalize')
