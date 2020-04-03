#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive Serial Runner component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
import os
import json

#-------------------------------------------------------------------------------
#
#  Dictionary of fastran outputs and their prefix codes.
#
#-------------------------------------------------------------------------------
parameter = {
    'q95'         : {'prefix' : 'aeqdisk',        'name' : 'q95',         'depends' : []                                                     },
    'peak'        : {'prefix' : 'aeqdisk',        'name' : 'peak',        'depends' : []                                                     },
    'li'          : {'prefix' : 'aeqdisk',        'name' : 'li',          'depends' : []                                                     },
    'qmin'        : {'prefix' : 'aeqdisk',        'name' : 'qmin',        'depends' : []                                                     },
    'betap'       : {'prefix' : 'aeqdisk',        'name' : 'betap',       'depends' : []                                                     },
    'rho_qmin'    : {'prefix' : 'aeqdisk',        'name' : 'rho_qmin',    'depends' : []                                                     },
    'betan_efit'  : {'prefix' : 'aeqdisk',        'name' : 'betan',       'depends' : []                                                     },
    'ne'          : {'prefix' : 'fastran',        'name' : 'ne',          'depends' : []                                                     },
    'ni'          : {'prefix' : 'fastran',        'name' : 'ni',          'depends' : []                                                     },
    'te'          : {'prefix' : 'fastran',        'name' : 'te',          'depends' : []                                                     },
    'ti'          : {'prefix' : 'fastran',        'name' : 'ti',          'depends' : []                                                     },
    'q'           : {'prefix' : 'fastran',        'name' : 'q',           'depends' : []                                                     },
    'rmajor'      : {'prefix' : 'fastran',        'name' : 'rmajor',      'depends' : []                                                     },
    'aminor'      : {'prefix' : 'fastran',        'name' : 'aminor',      'depends' : []                                                     },
    'zeff'        : {'prefix' : 'fastran',        'name' : 'zeff',        'depends' : []                                                     },
    'kappa'       : {'prefix' : 'fastran',        'name' : 'delta',       'depends' : []                                                     },
    'delta'       : {'prefix' : 'fastran',        'name' : 'kappa',       'depends' : []                                                     },
    'ip'          : {'prefix' : 'fastran',        'name' : 'ip',          'depends' : []                                                     },
    'bt'          : {'prefix' : 'fastran',        'name' : 'b0',          'depends' : []                                                     },
    'r'           : {'prefix' : 'fastran',        'name' : 'r0',          'depends' : []                                                     },
    'a'           : {'prefix' : 'fastran',        'name' : 'a0',          'depends' : []                                                     },
    'we'          : {'prefix' : 'fastran',        'name' : 'we',          'depends' : []                                                     },
    'wi'          : {'prefix' : 'fastran',        'name' : 'wi',          'depends' : []                                                     },
    'wb'          : {'prefix' : 'fastran',        'name' : 'wb',          'depends' : []                                                     },
    'betan'       : {'prefix' : 'fastran',        'name' : 'betan',       'depends' : []                                                     },
    'nebar'       : {'prefix' : 'fastran',        'name' : 'nebar',       'depends' : []                                                     },
    'tea'         : {'prefix' : 'fastran',        'name' : 'tea',         'depends' : []                                                     },
    'tia'         : {'prefix' : 'fastran',        'name' : 'tia',         'depends' : []                                                     },
    'taue'        : {'prefix' : 'fastran',        'name' : 'taue',        'depends' : []                                                     },
    'taui'        : {'prefix' : 'fastran',        'name' : 'taui',        'depends' : []                                                     },
    'tauth'       : {'prefix' : 'fastran',        'name' : 'tauth',       'depends' : []                                                     },
    'tautot'      : {'prefix' : 'fastran',        'name' : 'tautot',      'depends' : []                                                     },
    'tautot'      : {'prefix' : 'fastran',        'name' : 'tautot',      'depends' : []                                                     },
    'tau89'       : {'prefix' : 'fastran',        'name' : 'tau89',       'depends' : []                                                     },
    'tau98'       : {'prefix' : 'fastran',        'name' : 'tau98',       'depends' : []                                                     },
    'pnbe'        : {'prefix' : 'fastran',        'name' : 'pnbe',        'depends' : []                                                     },
    'pnbi'        : {'prefix' : 'fastran',        'name' : 'pnbi',        'depends' : []                                                     },
    'prfe'        : {'prefix' : 'fastran',        'name' : 'prfe',        'depends' : []                                                     },
    'prfi'        : {'prefix' : 'fastran',        'name' : 'prfi',        'depends' : []                                                     },
    'prad'        : {'prefix' : 'fastran',        'name' : 'prad',        'depends' : []                                                     },
    'pei'         : {'prefix' : 'fastran',        'name' : 'pei',         'depends' : []                                                     },
    'poh'         : {'prefix' : 'fastran',        'name' : 'poh',         'depends' : []                                                     },
    'ibs'         : {'prefix' : 'fastran',        'name' : 'ibs',         'depends' : []                                                     },
    'inb'         : {'prefix' : 'fastran',        'name' : 'inb',         'depends' : []                                                     },
    'irf'         : {'prefix' : 'fastran',        'name' : 'irf',         'depends' : []                                                     },
    'pfuse_equiv' : {'prefix' : 'fastran',        'name' : 'pfuse_equiv', 'depends' : []                                                     },
    'pfusi_equiv' : {'prefix' : 'fastran',        'name' : 'pfusi_equiv', 'depends' : []                                                     },
    'rho'         : {'prefix' : 'fastran',        'name' : 'rho',         'depends' : []                                                     },
    'rhob'        : {'prefix' : 'fastran',        'name' : 'rhob',        'depends' : []                                                     },
    'volp'        : {'prefix' : 'fastran',        'name' : 'volp',        'depends' : []                                                     },
    'pe_rf'       : {'prefix' : 'fastran',        'name' : 'pe_rf',       'depends' : []                                                     },
    'pfus0'       : {'prefix' : 'outone',         'name' : 'pfus0',       'depends' : []                                                     },
    'tauR'        : {'prefix' : 'get_tauR',       'name' : 'tauR',        'depends' : ['sigma', 'rho', 'rhob']                               },
    'sigma'       : {'prefix' : 'sigma_hirshman', 'name' : 'sigma',       'depends' : ['ne', 'te', 'aminor', 'rmajor', 'zeff', 'q']          },
    'neavg'       : {'prefix' : 'get_neavg',      'name' : 'sigma',       'depends' : ['volp', 'rhob', 'ne']                                 },
    'rho_ec'      : {'prefix' : 'get_rho_ec',     'name' : 'rho_ec',      'depends' : ['rho', 'pe_rf']                                       },
    'q0'          : {'prefix' : 'get_q0',         'name' : 'q0',          'depends' : ['q']                                                  },
    'nus_core'    : {'prefix' : 'nues',           'name' : 'nus_core',    'depends' : ['rmajor', 'aminor', 'ne', 'te', 'zeff', 'q', 'i_core']},
    'kappa1'      : {'prefix' : 'get_kappa1',     'name' : 'kappa1',      'depends' : ['kappa']                                              },
    'delta1'      : {'prefix' : 'get_delta1',     'name' : 'delta1',      'depends' : ['delta']                                              },
    'area'        : {'prefix' : 'get_area',       'name' : 'area',        'depends' : []                                                     },
    'area1'       : {'prefix' : 'get_area1',      'name' : 'area1',       'depends' : ['r', 'a', 'kappa', 'delta']                           },
    'pinj'        : {'prefix' : 'get_pinj',       'name' : 'pinj',        'depends' : ['prfe', 'prfi', 'pnbe', 'pnbi']                       },
    'pinj_ei'     : {'prefix' : 'get_pinj_ei',    'name' : 'pink_ei',     'depends' : ['prfe', 'pnbe', 'pinj']                               },
    'pfus'        : {'prefix' : 'get_pfus',       'name' : 'pfus',        'depends' : ['pfuse_equiv', 'pfusi_equiv']                         },
    'qdt'         : {'prefix' : 'get_qdt',        'name' : 'qdt',         'depends' : ['pfus', 'pinj']                                       },
    'wmhd'        : {'prefix' : 'get_wmhd',       'name' : 'wmhd',        'depends' : ['we', 'wi', 'wb']                                     },
    'fni'         : {'prefix' : 'get_fni',        'name' : 'fni',         'depends' : ['inb', 'irf', 'ibs' 'ip']                             },
    'fbs'         : {'prefix' : 'get_fbs',        'name' : 'fbs',         'depends' : ['ibs', 'ip']                                          },
    'fnb'         : {'prefix' : 'get_fnb',        'name' : 'fnb',         'depends' : ['inb', 'ip']                                          },
    'frf'         : {'prefix' : 'get_frf',        'name' : 'frf',         'depends' : ['irf', 'ip']                                          },
    'h98'         : {'prefix' : 'get_h98',        'name' : 'h98',         'depends' : ['tauth', 'tau98']                                     },
    'fgw'         : {'prefix' : 'get_fgw',        'name' : 'fgw',         'depends' : ['nebar', 'ip', 'a']                                   },
    'fgw_ped'     : {'prefix' : 'get_fgw_ped',    'name' : 'fgw_ped',     'depends' : ['neped', 'ip', 'a']                                   },
    'pmhd'        : {'prefix' : 'get_pmhd',       'name' : 'pmhd',        'depends' : ['b', 'ip', 'a', 'betan', 'mu0']                       },
    'nus_edge'    : {'prefix' : 'constant',       'name' : 'nus_edge',    'depends' : []                                                     },
    'i_core'      : {'prefix' : 'constant',       'name' : 'i_core',      'depends' : []                                                     },
    'mu0'         : {'prefix' : 'constant',       'name' : 'mu0',         'depends' : []                                                     }
}

constants = {
    'i_core'   : 50,
    'nus_edge' : 0.0,
    'mu0'      : math.pi*4.0E-7
}

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component Constructor
#
#-------------------------------------------------------------------------------
class massive_serial_runner(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component init method. This method prepairs the state.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: init')

#  Get config filenames.
        if (timeStamp == 0.0):
            self.current_state = self.services.get_config_param('CURRENT_MSR_STATE')
            self.current_batch = self.services.get_config_param('CURRENT_MSR_BATCH')
            self.model_config = self.services.get_config_param('MODEL_CONFIG')
            self.fastran_config = self.services.get_config_param('FASTRAN_CONFIG')

#  Stage state.
        self.services.stage_state()

#  Unzip files from the state. Use mode a so files an be read and written to.
        self.zip_ref = ZipState.ZipState(self.current_state, 'a')
        self.zip_ref.extract('inscan')

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component step method. Not used.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: step')

        flags = self.zip_ref.get_state()

#  Run the massive serial workflow.
        if 'state' in flags and flags['state'] == 'needs_update':
            with open(current_batch, 'r') as json_ref:
                keys = json.load(json_ref).keys()

            task_wait = self.services.launch_task(MAX(self.NPROC, len(keys)),
                                                  self.services.get_working_dir(),
                                                  self.services.MASSIVE_SERIAL_EXE,
                                                  'inscan',
                                                  self.fastran_config,
                                                  MAX(self.NPROC, len(keys)),
                                                  1)

#  Collect results of the workflow.
            if (self.services.wait_task(task_wait)):
                self.services.error('massive_serial_runner: step failed.')



#  Append output data to batch. Need to know what data to read.

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner Component finalize method. This cleans up
#  afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: finalize')

#-------------------------------------------------------------------------------
#
#  Load data parameter. Data parameters are listed in the parameters dictionary.
#  for a given parameter, start by loading the dependencies. Then based on the
#  prefix value load the data from the appropriate source. To avoid file io
#  overhead, a file will be first checked it is opened. Opened files will be
#  added to a list to cache the file reference. This updates the data base with
#  the newly loaded parameter.
#
#  FIXME: Not sure how I should update for the mutilple directories so start with just one.
#
#-------------------------------------------------------------------------------
    def load_parameter(self, directory, parameter, data_base):

#  First check that the parameter is not already loaded.
        if parameter in data_base:
            return data_base

#  Next load any dependencies.
        for key in parameters[parameter][depends]:
            data_base = self.load_parameter(directory, key, parameters)

#  Now load the parameter.

#-------------------------------------------------------------------------------
#
#  Load a parameter from a constant.
#
#-------------------------------------------------------------------------------
    def constant(parameter):
        return constants[parameter]

#-------------------------------------------------------------------------------
#
#
#
#-------------------------------------------------------------------------------
