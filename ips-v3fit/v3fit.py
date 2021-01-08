#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for V3FIT component. This wapper only takes a V3FIT input file
#  and runs V3FIT.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from omfit.classes.omfit_namelist import OMFITnamelist
from omfit.classes.omfit_nc import OMFITnc
from ips_component_utilities import ZipState
from ips_component_utilities import ScreenWriter
from ips_component_utilities import NamelistItem
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
        self.services.stage_state()

#  Get config filenames.
        self.current_v3fit_namelist = self.services.get_config_param('V3FIT_NAMELIST_INPUT')
        self.current_v3fit_state = self.services.get_config_param('CURRENT_V3FIT_STATE')
        self.result_file = 'result.{}.nc'.format(self.current_v3fit_namelist)
        current_siesta_namelist = self.services.get_config_param('SIESTA_NAMELIST_INPUT')
        current_siesta_state = self.services.get_config_param('CURRENT_SIESTA_STATE')
        current_vmec_namelist = self.services.get_config_param('VMEC_NAMELIST_INPUT')
        self.current_vmec_state = self.services.get_config_param('CURRENT_VMEC_STATE')
        self.current_wout_file = 'wout_{}.nc'.format(current_vmec_namelist.replace('input.','',1))

#  Stage state.
        self.services.stage_state()

#  Unzip files from the state. Use mode a so files can be read and written to.
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

                siesta_zip_ref.extract(self.current_vmec_state)

                with ZipState.ZipState(self.current_vmec_state, 'r') as vmec_zip_ref:
                    vmec_zip_ref.extract(self.current_wout_file)
                    flags = vmec_zip_ref.get_state()
                    if 'state' in flags and flags['state'] == 'updated':
                        self.zip_ref.set_state(state='needs_update')

                keywords['siesta_nli_filename'] = current_siesta_namelist
                keywords['siesta_restart_filename'] = current_restart_file
                keywords['vmec_nli_filename'] = current_vmec_namelist
                keywords['vmec_wout_input'] = self.current_wout_file
                keywords['model_eq_type'] = 'siesta'
        else:
            self.zip_ref.extract(self.current_vmec_state)

            with ZipState.ZipState(self.current_vmec_state, 'r') as vmec_zip_ref:
                vmec_zip_ref.extract(self.current_wout_file)
                vmec_zip_ref.extract(current_vmec_namelist)
                flags = vmec_zip_ref.get_state()
                if 'state' in flags and flags['state'] == 'updated':
                    self.zip_ref.set_state(state='needs_update')

            keywords['vmec_nli_filename'] = current_vmec_namelist
            keywords['vmec_wout_input'] = self.current_wout_file
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

        if 'state' in flags and flags['state'] == 'needs_update' or 'force_update' in keywords:
            task_wait = self.services.launch_task(self.NPROC,
                                                  self.services.get_working_dir(),
                                                  self.V3FIT_EXE,
                                                  self.current_v3fit_namelist,
                                                  logfile = 'v3fit_{}.log'.format(timeStamp))

#  Update flags.
            self.zip_ref.set_state(state='updated')

#  Wait for V3FIT to finish.
            if (self.services.wait_task(task_wait) and not os.path.exists(self.result_file)):
                self.services.error('v3fit: step failed.')

            if 'force_update' in keywords:
                with ZipState.ZipState(self.current_vmec_state, 'a') as vmec_zip_ref:
                    vmec_zip_ref.write(self.current_wout_file)
                self.zip_ref.write(self.current_vmec_state)

#  Add the result file to the state.
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
        self.services.update_state()

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
        namelist = OMFITnamelist(self.current_v3fit_namelist,
                                 collect_arrays={
                                 'pp_ne_b'              : {'default' : 0.0,        'shape' : (22,),       'offset' : (0,),    'sparray' : True},
                                 'pp_ne_as'             : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 'pp_ne_af'             : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 'pp_sxrem_b'           : {'default' : 0.0,        'shape' : (22,),       'offset' : (0,),    'sparray' : True},
                                 'pp_sxrem_as'          : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 'pp_sxrem_af'          : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 #'pp_sxrem_ptype_a'     : {'default' : '',         'shape' : (10,),       'offset' : (1,)},
                                 'pp_sxrem_b_a'         : {'default' : 0.0,        'shape' : (10, 22),    'offset' : (1,0),   'sparray' : True},
                                 'pp_sxrem_as_a'        : {'default' : 0.0,        'shape' : (10, 101),   'offset' : (1,1),   'sparray' : True},
                                 'pp_sxrem_af_a'        : {'default' : 0.0,        'shape' : (10, 101),   'offset' : (1,1),   'sparray' : True},
                                 'pp_te_b'              : {'default' : 0.0,        'shape' : (22,),       'offset' : (0,),    'sparray' : True},
                                 'pp_te_as'             : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 'pp_te_af'             : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 'pp_ti_b'              : {'default' : 0.0,        'shape' : (22,),       'offset' : (0,),    'sparray' : True},
                                 'pp_ti_as'             : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 'pp_ti_af'             : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 'pp_ze_b'              : {'default' : 0.0,        'shape' : (22,),       'offset' : (0,),    'sparray' : True},
                                 'pp_ze_as'             : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 'pp_ze_af'             : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 'sxrem_te_a'           : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 'sxrem_ratio_a'        : {'default' : 0.0,        'shape' : (101,),      'offset' : (1,),    'sparray' : True},
                                 #'model_sxrem_type_a'   : {'default' : '',         'shape' : (10,),       'offset' : (1,)},
                                 'sxrem_min'            : {'default' : 0.0,        'shape' : (10,),       'offset' : (1,),    'sparray' : True},
                                 'coosig_wgts'          : {'default' : 0.0,        'shape' : (100,),      'offset' : (1,),    'sparray' : True},
                                 #'rc_type'              : {'default' : '',         'shape' : (100,),      'offset' : (1,)},
                                 'rc_index'             : {'default' : 0,          'shape' : (100,),      'offset' : (1,),    'sparray' : True},
                                 'rc_value'             : {'default' : 0.0,        'shape' : (100,),      'offset' : (1,),    'sparray' : True},
                                 #'dp_type'              : {'default' : '',         'shape' : (100,),      'offset' : (1,)},
                                 'dp_index'             : {'default' : 0,          'shape' : (100,2),     'offset' : (1,1),   'sparray' : True},
                                 #'rp_type'              : {'default' : '',         'shape' : (100,),      'offset' : (1,)},
                                 'rp_index'             : {'default' : 0,          'shape' : (100,),      'offset' : (1,),    'sparray' : True},
                                 'rp_index2'            : {'default' : 0,          'shape' : (100,),      'offset' : (1,),    'sparray' : True},
                                 'rp_vrnc'              : {'default' : 0.0,        'shape' : (100,),      'offset' : (1,),    'sparray' : True},
                                 #'rp_range_type'        : {'default' : 'infinity', 'shape' : (100,2),     'offset' : (1,1)},
                                 'rp_range_value'       : {'default' : 0.0,        'shape' : (100,2),     'offset' : (1,1),   'sparray' : True},
                                 'rp_range_index'       : {'default' : 0,          'shape' : (100,2,2),   'offset' : (1,1,1), 'sparray' : True},
                                 #'lp_type'              : {'default' : '',         'shape' : (100,),      'offset' : (1,)},
                                 'lp_index'             : {'default' : 0,          'shape' : (100,),      'offset' : (1,),    'sparray' : True},
                                 'lp_index2'            : {'default' : 0,          'shape' : (100,),      'offset' : (1,),    'sparray' : True},
                                 #'lp_sets'              : {'default' : '',         'shape' : (100,100),   'offset' : (1,1)},
                                 'lp_sets_index'        : {'default' : 0,          'shape' : (100,100),   'offset' : (1,1),   'sparray' : True},
                                 'lp_sets_index2'       : {'default' : 0,          'shape' : (100,100),   'offset' : (1,1),   'sparray' : True},
                                 'lp_sets_coeff'        : {'default' : 0.0,        'shape' : (100,100),   'offset' : (1,1),   'sparray' : True},
                                 'sdo_data_a'           : {'default' : 0.0,        'shape' : (4012,),     'offset' : (1,),    'sparray' : True},
                                 'sdo_sigma_a'          : {'default' : 0.0,        'shape' : (4012,),     'offset' : (1,),    'sparray' : True},
                                 'sdo_weight_a'         : {'default' : 0.0,        'shape' : (4012,),     'offset' : (1,),    'sparray' : True},
                                 #'mag_a'                : {'default' : True,       'shape' : (1000,),     'offset' : (1,)},
                                 #'mag_3d_a'             : {'default' : False,      'shape' : (1000,),     'offset' : (1,)},
                                 'sdo_s_spec_imin'      : {'default' : 0,          'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'sdo_s_spec_imax'      : {'default' : 0,          'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'sdo_s_spec_floor'     : {'default' : 0.0,        'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'sdo_s_spec_fraction'  : {'default' : 0.0,        'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'sdo_w_spec_imin'      : {'default' : 0,          'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'sdo_w_spec_imax'      : {'default' : 0,          'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'sdo_w_spec_weight'    : {'default' : 0.0,        'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'mag_spec_imin'        : {'default' : 0,          'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'mag_spec_imax'        : {'default' : 0,          'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 #'mag_spec_use_induced' : {'default' : True,       'shape' : (150,),      'offset' : (1,)},
                                 'sfactor_spec_imin'    : {'default' : 0,          'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'sfactor_spec_imax'    : {'default' : 0,          'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'sfactor_spec_fac'     : {'default' : 0.0,        'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'soffset_spec_imin'    : {'default' : 0,          'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'soffset_spec_imax'    : {'default' : 0,          'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'soffset_spec_fac'     : {'default' : 0.0,        'shape' : (150,),      'offset' : (1,),    'sparray' : True},
                                 'n_phi_lif'            : {'default' : 0,          'shape' : (1000,),     'offset' : (1,),    'sparray' : True},
                                 'lif_arz'              : {'default' : 0.0,        'shape' : (1000,5,5),  'offset' : (1,0,0), 'sparray' : True},
                                 'lif_rc'               : {'default' : 0.0,        'shape' : (1000,),     'offset' : (1,),    'sparray' : True},
                                 'lif_zc'               : {'default' : 0.0,        'shape' : (1000,),     'offset' : (1,),    'sparray' : True},
                                 'lif_sigma'            : {'default' : 0.001,      'shape' : (1000,),     'offset' : (1,),    'sparray' : True},
                                 'lif_phi_degree'       : {'default' : 0.0,        'shape' : (1000,1000), 'offset' : (1,1),   'sparray' : True},
                                 #'lif_on_edge'          : {'default' : False,      'shape' : (1000,),     'offset' : (1,)},
                                 #'prior_name'           : {'default' : '',         'shape' : (1000,),     'offset' : (1,)},
                                 #'prior_param_name'     : {'default' : '',         'shape' : (1000,),     'offset' : (1,)},
                                 'prior_indices'        : {'default' : 0,          'shape' : (1000,2),    'offset' : (1,1),   'sparray' : True},
                                 #'prior_units'          : {'default' : '',         'shape' : (1000,),     'offset' : (1,)},
                                 'n_sig_coosig'         : {'default' : 0,          'shape' : (1000,),     'offset' : (1,),    'sparray' : True},
                                 'coosig_indices'       : {'default' : 0,          'shape' : (1000,100),  'offset' : (1,1),   'sparray' : True},
                                 'coosig_coeff'         : {'default' : 0.0,        'shape' : (1000,100),  'offset' : (1,1),   'sparray' : True},
                                 #'coosig_type'          : {'default' : '',         'shape' : (1000,),     'offset' : (1,)},
                                 #'coosig_name'          : {'default' : '',         'shape' : (1000,),     'offset' : (1,)},
                                 #'coosig_units'         : {'default' : '',         'shape' : (1000,),     'offset' : (1,)},
                                 'coosig_wgts_id'       : {'default' : -1,         'shape' : (1000,),     'offset' : (1,),    'sparray' : True},
                                 'n_gp_signal'          : {'default' : 0,          'shape' : (12,),       'offset' : (1,),    'sparray' : True},
                                 'gp_signal_indices'    : {'default' : 0,          'shape' : (12,4012),   'offset' : (1,1),   'sparray' : True},
                                 #'gp_model_type'        : {'default' : '',         'shape' : (12,),       'offset' : (1,)},
                                 'gp_model_index'       : {'default' : 0,          'shape' : (12,),       'offset' : (1,),    'sparray' : True},
                                 'gp_param_vrnc'        : {'default' : 0.001,      'shape' : (12,100),    'offset' : (1,1),   'sparray' : True},
                                 'gp_tolerance'         : {'default' : 1.0E-4,     'shape' : (12,),       'offset' : (1,),    'sparray' : True},
                                 'gp_cholesky_fact'     : {'default' : 0.0,        'shape' : (12,),       'offset' : (1,),    'sparray' : True}
                                 })

        if 'vmec_nli_filename' in keywords:
            namelist['v3fit_main_nli']['vmec_nli_filename'] = keywords['vmec_nli_filename']
            self.update = True
            del keywords['vmec_nli_filename']
        if 'my_task' in keywords and keywords['my_task'] == 'v3post' and 'vmec_wout_input' in keywords:
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

            for key, value in keywords.items():
                NamelistItem.set(namelist['v3fit_main_nli'], key, value)

        if namelist['v3fit_main_nli']['my_task'] == 'reconstruct' or namelist['v3fit_main_nli']['my_task'] == 'reconstruct_a1':
            del namelist['v3fit_main_nli']['vmec_wout_input']

        namelist.save()
