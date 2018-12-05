#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for V3FIT component. This wapper only takes a V3FIT input file
#  and runs V3FIT.
#
#-------------------------------------------------------------------------------

from component import Component
from omfit.classes.omfit_namelist import OMFITnamelist
from omfit.classes.omfit_nc import OMFITnc
from utilities import ZipState
from utilities import ScreenWriter
import json
import os

#-------------------------------------------------------------------------------
#
#  V3FIT Component Constructor
#
#-------------------------------------------------------------------------------
class v3fit(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  V3FIT Component init method. This method prepairs the namelist input file.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'v3fit: init')
        self.services.stage_plasma_state()

#  Get config filenames.
        self.current_v3fit_namelist = self.services.get_config_param('V3FIT_NAMELIST_INPUT')
        self.current_v3fit_state = self.services.get_config_param('CURRENT_V3FIT_STATE')
        self.result_file = 'result.{}.nc'.format(self.current_v3fit_namelist)
        current_siesta_namelist = self.services.get_config_param('SIESTA_NAMELIST_INPUT')
        current_siesta_state = self.services.get_config_param('CURRENT_SIESTA_STATE')
        current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        current_wout_file = 'wout_{}.nc'.format(current_vmec_namelist.replace('input.','',1))

#  Stage plasma state.
        self.services.stage_plasma_state()

#  Unzip files from the plasma state. Use mode a so files can be read and
#  written to.
        self.zip_ref = ZipState.ZipState(self.current_v3fit_state, 'a')
        self.zip_ref.extract(self.current_v3fit_namelist)
        if self.result_file in self.zip_ref:
            self.zip_ref.extract(self.result_file)

        if current_siesta_state in self.zip_ref:
            self.zip_ref.extract(current_siesta_state)

            with ZipState.ZipState(current_siesta_state, 'r') as siesta_zip_ref:
                siesta_zip_ref.extract(current_siesta_namelist)
                namelist = OMFITnamelist(current_siesta_namelist)
                current_restart_file = 'siesta_{}.nc'.format(namelist['siesta_info']['restart_ext'])

                siesta_zip_ref.extract(current_restart_file)
                flags = siesta_zip_ref.get_state()
                if 'state' in flags and flags['state'] == 'updated':
                    self.zip_ref.set_state(state='needs_update')

                siesta_zip_ref.extract(current_vmec_state)

                with ZipState.ZipState(current_vmec_state, 'r') as vmec_zip_ref:
                    vmec_zip_ref.extract(current_wout_file)
                    flags = vmec_zip_ref.get_state()
                    if 'state' in flags and flags['state'] == 'updated':
                        self.zip_ref.set_state(state='needs_update')

                keywords['siesta_nli_filename'] = current_siesta_namelist
                keywords['siesta_restart_filename'] = current_restart_file
                keywords['vmec_nli_filename'] = current_vmec_namelist
                keywords['vmec_wout_input'] = current_wout_file
                keywords['model_eq_type'] = 'siesta'
        else:
            self.zip_ref.extract(current_vmec_state)

            with ZipState.ZipState(current_vmec_state, 'r') as vmec_zip_ref:
                vmec_zip_ref.extract(current_wout_file)
                vmec_zip_ref.extract(current_vmec_namelist)
                flags = vmec_zip_ref.get_state()
                if 'state' in flags and flags['state'] == 'updated':
                    self.zip_ref.set_state(state='needs_update')

            keywords['vmec_nli_filename'] = current_vmec_namelist
            keywords['vmec_wout_input'] = current_wout_file
            keywords['model_eq_type'] = 'vmec'

#  Update parameters in the namelist.
        self.set_namelist(**keywords)

#-------------------------------------------------------------------------------
#
#  V3FIT Component step method. This runs V3FIT.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'v3fit: step')

        flags = self.zip_ref.get_state()

        if 'state' in flags and flags['state'] == 'needs_update':
            self.set_namelist(my_task='v3post')
            self.task_wait = self.services.launch_task(self.NPROC,
                                                       self.services.get_working_dir(),
                                                       self.V3FIT_EXE,
                                                       self.current_v3fit_namelist,
                                                       logfile = 'v3fit.log')

#  Update flags.
            self.zip_ref.set_state(state='updated')

#  Wait for V3FIT to finish.
            if (self.services.wait_task(self.task_wait) and not os.path.exists(self.result_file)):
                self.services.error('v3fit: step failed.')

#  Add the result file to the plasma state.
            self.zip_ref.write([self.current_v3fit_namelist, self.result_file])

        else:
#  Update flags.
            self.zip_ref.set_state(state='unchanged')

        if 'result_file' in keywords:
            result_nc = OMFITnc(self.result_file)
            nsteps = result_nc['nsteps']['data']
            result = {'signal_model': result_nc['signal_model_value']['data'][nsteps,:,0].tolist()}
            with open(keywords['result_file'], 'w') as result_ref:
                json.dump(result, result_ref)
            self.zip_ref.write(keywords['result_file'])

        self.zip_ref.close()
        self.services.update_plasma_state()

#-------------------------------------------------------------------------------
#
#  V3FIT Component finalize method. This cleans up afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'v3fit: finalize')

#-------------------------------------------------------------------------------
#
#  V3FIT Component set_namelist method. This sets the namelist input file from
#  the keywords.
#
#-------------------------------------------------------------------------------
    def set_namelist(self, **keywords):
#  Update parameters in the namelist.
        namelist = OMFITnamelist(self.current_v3fit_namelist)

        if 'vmec_nli_filename' in keywords:
            namelist['v3fit_main_nli']['vmec_nli_filename'] = keywords['vmec_nli_filename']
            self.update = True
            del keywords['vmec_nli_filename']
        if 'vmec_wout_input' in keywords:
            namelist['v3fit_main_nli']['vmec_wout_input'] = keywords['vmec_wout_input']
            self.update = True
            del keywords['vmec_wout_input']
        if 'siesta_nli_filename' in keywords:
            namelist['v3fit_main_nli']['siesta_nli_filename'] = keywords['siesta_nli_filename']
            self.update = True
            del keywords['siesta_nli_filename']
        if 'siesta_restart_filename' in keywords:
            namelist['v3fit_main_nli']['siesta_restart_filename'] = keywords['siesta_restart_filename']
            self.update = True
            del keywords['siesta_restart_filename']
        if 'model_eq_type' in keywords:
            namelist['v3fit_main_nli']['model_eq_type'] = keywords['model_eq_type']
            self.update = True
            del keywords['model_eq_type']
        if 'my_task' in keywords:
            namelist['v3fit_main_nli']['my_task'] = keywords['my_task']
            self.update = True
            del keywords['my_task']

        if len(keywords) > 0:
            self.zip_ref.set_state(state='needs_update')

            for key, value in keywords.iteritems():
                namelist['v3fit_main_nli'][key] = value

        namelist.save()
