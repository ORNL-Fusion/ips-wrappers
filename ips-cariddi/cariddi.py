#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for V3FIT component. This wapper only takes a V3FIT input file
#  and runs V3FIT.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
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
        if timeStamp == 0.0:
            self.current_cariddi_input = self.services.get_config_param('CARIDDI_INPUT')
            self.cariddi_matrix_path = self.services.get_config_param('CARIDDI_MATRIX_PATH')
            self.current_cariddi_state = self.services.get_config_param('CURRENT_CARIDDI_STATE')

            self.current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
            current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
            self.current_wout_file = 'wout_{}.nc'.format(current_vmec_namelist.replace('input.','',1))

            self.mgrid_file = self.services.get_config_param('MGRID_FILE')

            with ZipState.ZipState(self.current_cariddi_state, 'r') as zip_ref:
                zip_ref.extract(self.current_cariddi_input)

# Extract input files.
        with ZipState.ZipState(self.current_cariddi_state, 'r') as zip_ref:
            zip_ref.extract(self.current_vmec_state)

        with ZipState.ZipState(self.current_vmec_state, 'r') as zip_ref:
            if self.current_wout_file in zip_ref:
                zip_ref.extract(self.current_wout_file)

#-------------------------------------------------------------------------------
#
#  CARIDDI Component step method. This runs V3FIT.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi: step')

        if os.path.exists(self.current_wout_file):
#  Get the surface current.
            task_wait = self.services.launch_task(self.NPROC,
                                                  self.services.get_working_dir(),
                                                  self.SURFACE_EXE,
                                                  '-woutf={}'.format(self.current_wout_file),
                                                  '-surff={}'.format(self.current_cariddi_input),
                                                  '-para=-1',
                                                  logfile = 'surface.log')

#  Wait for SURFACE to finish.
            self.services.wait_task(task_wait)

            task_wait = self.services.launch_task(self.NPROC,
                                                  self.services.get_working_dir(),
                                                  self.CARIDDI_EXE,
                                                  '--mgrid_file={}'.format(self.mgrid_file),
                                                  '--matrix_path={}'.format(self.cariddi_matrix_path),
                                                  '--vmec_current={}'.format(self.current_cariddi_input),
                                                  logfile = 'cariddi.log')
        else:
#  Clear eddy currents from the mgrid file.
            task_wait = self.services.launch_task(self.NPROC,
                                                  self.services.get_working_dir(),
                                                  self.CARIDDI_EXE,
                                                  '--mgrid_file={}'.format(self.mgrid_file),
                                                  '--matrix_path={}'.format(self.cariddi_matrix_path),
                                                  logfile = 'cariddi.log')

        self.services.wait_task(task_wait)

#-------------------------------------------------------------------------------
#
#  CARIDDI Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'cariddi: finalize')
