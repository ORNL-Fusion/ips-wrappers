#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive Serial Runner init component.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ScreenWriter
from utilities import ZipState
import os
import json

#-------------------------------------------------------------------------------
#
#  Dictionary of fastran inputs and their prefix codes.
#
#-------------------------------------------------------------------------------
parameter = {
    'r'           : {'prefix' : 'fastran_init', 'name' : 'r0'              },
    'a'           : {'prefix' : 'fastran_init', 'name' : 'a0'              },
    'kappa'       : {'prefix' : 'fastran_init', 'name' : 'kappa'           },
    'delta'       : {'prefix' : 'fastran_init', 'name' : 'delta'           },
    'bt'          : {'prefix' : 'fastran_init', 'name' : 'b0'              },
    'ip'          : {'prefix' : 'fastran_init', 'name' : 'ip'              },
    'zeffped'     : {'prefix' : 'fastran_init', 'name' : 'zeffped'         },
    'm'           : {'prefix' : 'fastran_init', 'name' : 'm'               },
    'xmid'        : {'prefix' : 'fastran_init', 'name' : 'xmid'            },
    'xwid'        : {'prefix' : 'fastran_init', 'name' : 'xwid'            },
    'ne0'         : {'prefix' : 'fastran_init', 'name' : 'ne0'             },
    'neped'       : {'prefix' : 'fastran_init', 'name' : 'neped'           },
    'nesep'       : {'prefix' : 'fastran_init', 'name' : 'nesep'           },
    'alpha_ne'    : {'prefix' : 'fastran_init', 'name' : 'alpha_ne'        },
    'beta_ne'     : {'prefix' : 'fastran_init', 'name' : 'beta_ne'         },
    'te0'         : {'prefix' : 'fastran_init', 'name' : 'te0'             },
    'teped'       : {'prefix' : 'fastran_init', 'name' : 'teped'           },
    'tesep'       : {'prefix' : 'fastran_init', 'name' : 'tesep'           },
    'alpha_te'    : {'prefix' : 'fastran_init', 'name' : 'alpha_te'        },
    'beta_te'     : {'prefix' : 'fastran_init', 'name' : 'beta_te'         },
    'ti0'         : {'prefix' : 'fastran_init', 'name' : 'ti0'             },
    'tiped'       : {'prefix' : 'fastran_init', 'name' : 'tiped'           },
    'tisep'       : {'prefix' : 'fastran_init', 'name' : 'tisep'           },
    'alpha_ti'    : {'prefix' : 'fastran_init', 'name' : 'alpha_ti'        },
    'beta_ti'     : {'prefix' : 'fastran_init', 'name' : 'beta_ti'         },
    'omega0'      : {'prefix' : 'fastran_init', 'name' : 'omega0'          },
    'alpha_omega' : {'prefix' : 'fastran_init', 'name' : 'alpha_omega'     },
    'beta_omega'  : {'prefix' : 'fastran_init', 'name' : 'beta_omega'      },
    'zeff0'       : {'prefix' : 'fastran_init', 'name' : 'zeff0'           },
    'pe'          : {'prefix' : 'hcd_model'   , 'name' : 'pe'              },
    'pi'          : {'prefix' : 'hcd_model'   , 'name' : 'pi'              },
    'pinja_0'     : {'prefix' : 'nubeam'      , 'name' : 'pinja_0'         },
    'pinja_1'     : {'prefix' : 'nubeam'      , 'name' : 'pinja_1'         },
    'rfpow_0'     : {'prefix' : 'toray'       , 'name' : 'rfpow_0'         },
    'thet_0'      : {'prefix' : 'toray'       , 'name' : 'rfpow_0'         },
    'phai_0'      : {'prefix' : 'toray'       , 'name' : 'phai_0'          },
    'thgrill_0'   : {'prefix' : 'genray'      , 'name' : 'grill_thgrill_0' },
    'frqncy_0'    : {'prefix' : 'genray'      , 'name' : 'WAVE_FRQNCY_0'   },
    'powers_0'    : {'prefix' : 'genray'      , 'name' : 'GRILL_POWERS_0'  },
    'rho_jpeak'   : {'prefix' : 'modeleq_init', 'name' : 'RHO_JPEAK'       },
    'jaxis'       : {'prefix' : 'modeleq_init', 'name' : 'JAXIS'           }
}

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component Constructor
#
#-------------------------------------------------------------------------------
class massive_serial_runner_init(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component init method. This method prepairs the
#  state.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_init: init')

#  Get config filenames.
        current_batch = self.services.get_config_param('CURRENT_MSR_BATCH')
        current_state = self.services.get_config_param('CURRENT_MSR_STATE')
        model_config = self.services.get_config_param('MODEL_CONFIG')

#  Remove old inputs.
        for file in os.listdir('.'):
            os.remove(file)

#  Stage input files and setup inital state.
        self.services.stage_input_files(self.INPUT_FILES)

        if os.path.exists(current_batch):
            if not os.path.exists(current_batch):
                raise Exception('Model config {} not found.'.format(model_config))

            with open(current_batch, 'r') as json_ref:
                create_input(json.load(json_ref))

#  Create plasma state from files. Input files can either be a new plasma state,
#  training data file or both. If both file were staged, replace the training
#  data input file. If the training data file is present flag the plasma state
#  as needing to be updated.
        with ZipState.ZipState(current_state, 'a') as zip_ref:
            if os.path.exists('inscan'):
                zip_ref.write('inscan')
                zip_ref.write(current_batch)
                zip_ref.write(model_config)
                zip_ref.set_state(state='needs_update')

        self.services.update_state()

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component step method. Not used.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_init: step')

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component finalize method. This cleans up
#  afterwards. Not used.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner_init: finalize')

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component create_input method. This takes an data
#  base of inputs and formats the inscan file.
#
#-------------------------------------------------------------------------------
    def create_input(data_base):
        keys = data_base.keys()

        input = '{}'.format(set_header('TIME_ID'))

        for key in keys:
            input = '{} {}'.format(input, set_header(key))

        input = '{}\n'.format(input)

        for i in range(len(data_base[keys[0]])):
            input = '{}{:05d}'.format(input, i)

            for key in keys:
                input = '{} {}'.format(input, data_base[key][i])

            input = '{}\n'.format(input)

        with open('inscan', 'w') as input_ref:
            input_ref.write(input)

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component set_header method. Sets the header for
#  a key. This method assumes the key is a float.
#
#-------------------------------------------------------------------------------
    def set_header(key):
        return '{}:{}:{}'.format(get_prefix(key), key.toupper(), 'float')

#-------------------------------------------------------------------------------
#
#  Massive Serial Runner init Component get_prefix method. This gets the prefix
#  code for the header. Prefixes for the key are defined in the parameter
#  dictionary.
#
#-------------------------------------------------------------------------------
    def get_prefix(key):
        return parameter[key.tolower()]['prefix']
