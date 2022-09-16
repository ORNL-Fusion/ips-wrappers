#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for the Massive VMEC Runner component.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from omfit_classes.omfit_namelist import OMFITnamelist
from ips_component_utilities import ZipState
from ips_component_utilities import ScreenWriter
from ips_component_utilities import NamelistItem
from ml import adaptive
import shutil
import netCDF4
import os
import re

#-------------------------------------------------------------------------------
#
#  Massive VMEC Runner Component Constructor
#
#-------------------------------------------------------------------------------
class massive_vmec_runner(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)

#-------------------------------------------------------------------------------
#
#  Massive VMEC Runner Component init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: init')

#  Stage state.
        self.services.stage_state()

        if timeStamp == 0.0:
            self.current_state = self.services.get_config_param('CURRENT_MVR_STATE')
            self.current_batch = self.services.get_config_param('CURRENT_BATCH')
            self.constraint_path = self.services.get_config_param('MODULE_PATH')
            self.constraint_name = self.services.get_config_param('MODULE_NAME')
            self.batch_size = int(self.services.get_config_param('BATCH_SIZE'))
            self.namelist_template = self.services.get_config_param('VMEC_NAMELIST_INPUT')
            self.cores_per_node = self.services.get_config_param('CORES_PER_NODE')

#  Unzip files from the state. Use mode a so files can be read and written to.
        self.zip_ref = ZipState.ZipState(self.current_state, 'a')
        self.zip_ref.extract(self.current_batch)
        
        if timeStamp == 0.0:
            model_config = self.services.get_config_param('MODEL_CONFIG')
            self.zip_ref.extract(model_config)
            self.model_config = adaptive.load_json(model_config)
            self.zip_ref.extract(self.namelist_template)

        self.batch_data = adaptive.load_json(self.current_batch)

        self.services.create_task_pool('vmec_pool')

        for i in range(self.batch_size):
            namelist_name = 'input.{}.vmec'.format(i)
            
            if timeStamp == 0.0:
                shutil.copy(self.namelist_template, namelist_name)
                self.services.add_task('vmec_pool', namelist_name, 1, os.getcwd(), self.VMEC_EXE, namelist_name, logfile='{}.log'.format(namelist_name))

            namelist = OMFITnamelist(namelist_name,
                                     collect_arrays={
                                        'ns_array'    : {'default' : 0, 'shape' : (100,),    'offset' : (1,),     'sparray' : True},
                                        'niter_array' : {'default' : 0, 'shape' : (100,),    'offset' : (1,),     'sparray' : True},
                                        'rbs'         : {'default' : 0, 'shape' : (203,101), 'offset' : (-101,0), 'sparray' : True},
                                        'rbc'         : {'default' : 0, 'shape' : (203,101), 'offset' : (-101,0), 'sparray' : True},
                                        'zbs'         : {'default' : 0, 'shape' : (203,101), 'offset' : (-101,0), 'sparray' : True},
                                        'zbc'         : {'default' : 0, 'shape' : (203,101), 'offset' : (-101,0), 'sparray' : True},
                                        'am'          : {'default' : 0, 'shape' : (21,),     'offset' : (0,),     'sparray' : True},
                                        'ai'          : {'default' : 0, 'shape' : (21,),     'offset' : (0,),     'sparray' : True},
                                        'ac'          : {'default' : 0, 'shape' : (21,),     'offset' : (0,),     'sparray' : True},
                                        'am_aux_s'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                        'am_aux_f'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                        'ai_aux_s'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                        'ai_aux_f'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                        'ac_aux_s'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                        'ac_aux_f'    : {'default' : 0, 'shape' : (10001,),  'offset' : (1,),     'sparray' : True},
                                        'raxis'       : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                        'zaxis'       : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                        'raxis_cc'    : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                        'raxis_cs'    : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                        'zaxis_cc'    : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                        'zaxis_cs'    : {'default' : 0, 'shape' : (102,),    'offset' : (0,),     'sparray' : True},
                                        'ftol_array'  : {'default' : 0, 'shape' : (100,),    'offset' : (1,),     'sparray' : True},
                                        'extcur'      : {'default' : 0, 'shape' : (300,),    'offset' : (1,),     'sparray' : True}
                                     })
        
            for key, value in self.batch_data['inputs'].items():
                NamelistItem.set(namelist['indata'], key, value[i])
                
            namelist.save()

#-------------------------------------------------------------------------------
#
#  Massive VMEC Runner Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0, **keywords):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: step')

        self.services.submit_tasks('vmec_pool')

        self.batch_data['outputs'] = {}
        
        self.services.get_finished_tasks('vmec_pool')
 
        for i in range(self.batch_size):
            wout_name = 'wout_{}.vmec.nc'.format(i)
            print(wout_name)

            with netCDF4.Dataset(wout_name, 'r') as wout_ref:
                for entry in self.model_config['outputs']:
                    if entry['name'] not in self.batch_data['outputs']:
                        print(entry['name'] not in self.batch_data['outputs'])
                        self.batch_data['outputs'] = {entry['name'] : []}
                    name, index = parse_name(entry['name'])
                    self.batch_data['outputs'][entry['name']].append(wout_ref.variables[name][:][index])

        adaptive.save_json(self.current_batch, self.batch_data)

        self.zip_ref.write(self.current_batch)
        self.zip_ref.close()

#-------------------------------------------------------------------------------
#
#  Massive VMEC Runner Component finalsize method.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'massive_serial_runner: finalize')
        self.services.remove_task_pool('vmec_pool')

#-------------------------------------------------------------------------------
#
#  Parse name to separate the indicies from the name.
#
#-------------------------------------------------------------------------------
def parse_name(name):    
    return re.search('(.*)\(', name).group()[:-1], eval(re.search('\((.*)\)', name).group())
