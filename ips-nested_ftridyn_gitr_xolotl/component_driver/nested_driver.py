
#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for NESTED Driver component.
#
#-------------------------------------------------------------------------------

from component import Component
import os
import shutil

#-------------------------------------------------------------------------------
#
#  NESTED Driver Class
#
#-------------------------------------------------------------------------------
class nested_driver(Component):

#-------------------------------------------------------------------------------
#
#  nested_driver Component Constructor
#
#-------------------------------------------------------------------------------
    def __init__(self, services, config):
        print('nested_driver: Construct')
        Component.__init__(self, services, config)
        self.async_queue = {}
    
#-------------------------------------------------------------------------------
#
#  nested_driver Component init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('nested_driver: init')

        self.nested_components = {}

        #stage input files
        self.services.stage_input_files(self.INPUT_FILES)
    
#  Sub workflows require manual setup. First a sub directory must be created.
#  Then copying of the input files must be performed manually. The first
#  argument of create sub workflow doesn't appear to do anything.

        self.nested_components['component_a'] = {'sim_name': None, 'init': None, 'driver': None, 'sub_working_dir': 'component_a_init'}
        print(self.nested_components['component_a'])
        print(self.nested_components['component_a']['sub_working_dir'])
 
        os.mkdir(self.nested_components['component_a']['sub_working_dir'])
        #shutil.copy2(self.services.get_config_param('COMPONENT_A_NAMELIST_INPUT'), self.nested_components['component_a']['sub_working_dir'])
        (self.nested_components['component_a']['sim_name'],
         self.nested_components['component_a']['init'],
         self.nested_components['component_a']['driver']) = self.services.create_sub_workflow('component_a_sub',
                                                                                              self.services.get_config_param('COMPONENT_A_CONF'),
                                                                                              {'PWD' : self.services.get_config_param('PWD'),
                                                                                               #'COMPONENT_A_NAMELIST_INPUT' : self.services.get_config_param('COMPONENT_A_NAMELIST_INPUT'),
                                                                                               'LOG_FILE' : 'log.component_a.warning'
                                                                                           })
                                                                                              
        self.nested_components['component_ftx'] = {'sim_name': None, 'init': None, 'driver': None, 'sub_working_dir': 'xolotl-ftridyn_driver'}
        os.mkdir(self.nested_components['component_ftx']['sub_working_dir'])
        #shutil.copy2(self.services.get_config_param('GITR_OUTPUT'), self.nested_components['component_ftx']['sub_working_dir'])
        paramFiles=self.services.get_config_param('XOLOTL_INPUT_FILES').split()
        for f in paramFiles:
            shutil.copy2(f, self.nested_components['component_ftx']['sub_working_dir'])

        (self.nested_components['component_ftx']['sim_name'],
         self.nested_components['component_ftx']['init'],
         self.nested_components['component_ftx']['driver']) = self.services.create_sub_workflow('ftx_sub',
                                                                                                self.services.get_config_param('FTX_CONF'),
                                                                                                {'PWD'                   : self.services.get_config_param('PWD'),
                                                                                                 'XOLOTL_INPUT_FILES' : self.services.get_config_param('XOLOTL_INPUT_FILES'),
                                                                                                 'GITR_OUTPUT'  : self.services.get_config_param('GITR_OUTPUT'),
                                                                                                 'LOG_FILE'                   : 'log.ftx.warning'
                                                                                             })
        
    
#-------------------------------------------------------------------------------
#
#  nested_driver Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('nested_driver: step')

        #  Inialized the plasma state by calling the init components of the sub
        #  workflows.

        #  Initalize the sub workflow drivers. We want to pass arguments to the sub
        #  workflow components. Create a dictionary of dictionaries for the stuff to
        #  override. Need some standard naming convection here.

        #  Step the drivers. Here Component a doesn't depend on component b so both can 
        #  be run in parallel.


        print('\n')
        print('FTridyn-GITR:init')
        print('\n')
        self.async_queue['component_a:driver:init'] = self.services.call(self.nested_components['component_a']['driver'], 'init', timeStamp)
        #self.async_queue['component_a:driver:init'] = self.services.call_nonblocking(self.nested_components['component_a']['driver'], 'init', timeStamp)
                                                                                     #, override = {'message': 'called from nested'})

        print('\n')
        print('FTridyn-GITR:step')
        print('\n')
        #self.services.wait_call_list([self.async_queue['component_a:driver:init']], True)
        del self.async_queue['component_a:driver:init']
        self.async_queue['component_a:driver:step'] = self.services.call(self.nested_components['component_a']['driver'], 'step', 0.0)
        #self.async_queue['component_a:driver:step'] = self.services.call_nonblocking(self.nested_components['component_a']['driver'], 'step', 0.0)


        gitrOutFile=self.services.get_config_param('GITR_OUTPUT')
        print ' '
        print 'Copy GITRs output ', gitrOutFile
        print '\t from ', self.SUBMIT_DIR  
        print '\t to ', self.nested_components['component_ftx']['sub_working_dir']
        print ' '
        shutil.copy2(self.SUBMIT_DIR+'/'+gitrOutFile, self.nested_components['component_ftx']['sub_working_dir'])


        print('\n')
        print('FTridyn-Xolotl:init')
        print('\n')
        self.async_queue['component_ftx:driver:init'] = self.services.call(self.nested_components['component_ftx']['driver'], 'init', timeStamp)
        #self.async_queue['component_ftx:driver:init'] = self.services.call_nonblocking(self.nested_components['component_ftx']['driver'],'init',timeStamp)
                                                                                     #, override = {'message': 'called from nested'})

        print('\n')
        print('FTridyn-Xolotl:step')
        print('\n')
        #self.services.wait_call_list([self.async_queue['component_ftx:driver:init']], True)
        del self.async_queue['component_ftx:driver:init']
        self.async_queue['component_ftx:driver:step'] = self.services.call(self.nested_components['component_ftx']['driver'], 'step', 0.0)
        #self.async_queue['component_ftx:driver:step'] = self.services.call_nonblocking(self.nested_components['component_ftx']['driver'], 'step', 0.0)
    
#-------------------------------------------------------------------------------
#
#  nested_driver Component finalize method.
#
#-------------------------------------------------------------------------------

    def finalize(self, timeStamp=0.0):
        print('nested_driver: finalize')

#  Wait until everything is finished before finalizing.
        self.services.wait_call_list(self.async_queue.values(), True)
        self.async_queue = {}

#  Call the finalize methods.
        self.async_queue['component_a:driver:finalize'] = self.services.call_nonblocking(self.nested_components['component_a']['driver'], 'finalize', timeStamp)
        self.async_queue['component_ftx:driver:finalize'] = self.services.call_nonblocking(self.nested_components['component_ftx']['driver'], 'finalize', timeStamp)

#  Wait until everything is finished before finalizing.
        self.services.wait_call_list(self.async_queue.values(), True)
        self.async_queue = {}
