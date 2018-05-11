#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for V3FIT init component. This wapper only takes a V3FIT input
#  file and runs V3FIT.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ZipState
import os

#-------------------------------------------------------------------------------
#
#  V3FIT init Component Constructor
#
#-------------------------------------------------------------------------------
class v3fit_init(Component):
    def __init__(self, services, config):
        print('v3fit_init: Construct')
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  V3FIT init Component init method. This method prepairs the namelist input
#  file and creates a dummy out put file. This allows staging the plasma state
#  files. In the v3fit namelist input file configure the v3fit namelist input
#  with the task, internal vmec input name and optional name of the wout file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('v3fit_init: init')

#  Get config filenames.
        current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        current_siesta_namelist = self.services.get_config_param('SIESTA_NAMELIST_INPUT')
        current_siesta_state = self.services.get_config_param('CURRENT_SIESTA_STATE')
        current_v3fit_namelist = self.services.get_config_param('V3FIT_NAMELIST_INPUT')
        current_v3fit_state = self.services.get_config_param('CURRENT_V3FIT_STATE')

#  Stage input files. Remove an old namelist input if it exists.
        if os.path.exists(current_v3fit_namelist):
            os.remove(current_v3fit_namelist)
        if os.path.exists(current_siesta_namelist):
            os.remove(current_siesta_namelist)
        if os.path.exists(current_vmec_namelist):
            os.remove(current_vmec_namelist)
        self.services.stage_input_files(self.INPUT_FILES)

#  Create plasma state from files. Input files can either be a new plasma state,
#  namelist input file or both. If both files were staged, replace the namelist
#  input file. If the namelist file is present flag the plasma state as needing
#  to be updated.
        with ZipState.ZipState(current_v3fit_state, 'a') as zip_ref:
            if os.path.exists(current_v3fit_namelist):
                zip_ref.write(current_v3fit_namelist)
                zip_ref.set_state(state='needs_update')

            if os.path.exists(current_siesta_state) or os.path.exists(current_siesta_namelist):
                with ZipState.ZipState(current_siesta_state, 'a') as zip_siesta_ref:
                    if os.path.exists(current_siesta_namelist):
                        zip_siesta_ref.write(current_siesta_namelist)
                        zip_siesta_ref.set_state(state='needs_update')
                    with ZipState.ZipState(current_vmec_state, 'a') as zip_vmec_ref:
                        if os.path.exists(current_vmec_namelist):
                            zip_vmec_ref.write(current_vmec_namelist)
                            zip_vmec_ref.set_state(state='needs_update')
                    zip_siesta_ref.write(current_vmec_state)
                zip_ref.write(current_siesta_state)
            else:
                with ZipState.ZipState(current_vmec_state, 'a') as zip_vmec_ref:
                    if os.path.exists(current_vmec_namelist):
                        zip_vmec_ref.write(current_vmec_namelist)
                        zip_vmec_ref.set_state(state='needs_update')
                zip_ref.write(current_vmec_state)

        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  V3FIT init Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('v3fit_init: step')

#-------------------------------------------------------------------------------
#
#  V3FIT init Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('v3fit_init: finalize')
