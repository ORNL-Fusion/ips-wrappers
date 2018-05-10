#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for PARENT Driver component.
#
#-------------------------------------------------------------------------------

from component import Component
import os
import shutil
import subprocess

import sys
#from contextlib import redirect_stdout

#-------------------------------------------------------------------------------
#
#  PARENT Driver Class
#
#-------------------------------------------------------------------------------
class parent_driver(Component):

#-------------------------------------------------------------------------------
#
#  parent_driver Component Constructor
#
#-------------------------------------------------------------------------------
    def __init__(self, services, config):
        print('parent_driver: Construct')
        Component.__init__(self, services, config)
        #dictionaries to store the FT-X and FT-GITR workflows and parametrization (config file, etc.) 
        self.running_components = {}
        self.child_components = {}
        self.async_queue = {}
        self.nested_components = {} 

#-------------------------------------------------------------------------------
#
#  parent_driver Component init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print ' '
        print('parent_driver: init')
 
        import sys
        print 'running python ', sys.version_info
   
        num_children = int(self.services.get_config_param('NUM_CHILDREN'))
        child_conf = self.services.get_config_param('CHILD_COMPONENT_CONF')
        
        keys = {
               'PWD' : self.services.get_config_param('PWD')
               }


        #stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        #  Sub workflows require manual setup. First a sub directory must be created.
        #  Then copying of the input files must be performed manually. The first
        #  argument of create sub workflow doesn't appear to do anything.


        #CREATE FT-GITR SUB-WORKFLOW: 
        
        self.nested_components['component_a'] = {'sim_name': None, 'init': None, 'driver': None, 'sub_working_dir': 'component_a_init'}
        print(self.nested_components['component_a'])
        print(self.nested_components['component_a']['sub_working_dir'])

        os.mkdir(self.nested_components['component_a']['sub_working_dir'])
        #shutil.copy2(self.services.get_config_param('COMPONENT_A_NAMELIST_INPUT'), self.child_components['component_a']['sub_working_dir'])
        print 'Creating sub workflo with: '
        print '\t name: component_a'
        print '\t config: ', self.services.get_config_param('COMPONENT_A_CONF')
        print '\t PWD: ', self.services.get_config_param('PWD')
        print '\t COMPONENT_A_NAMELIST_INPUT ', self.services.get_config_param('COMPONENT_A_NAMELIST_INPUT')
        print '\t LOG_FILE: log.component_a.warning'

        (self.nested_components['component_a']['sim_name'],
         self.nested_components['component_a']['init'],
         self.nested_components['component_a']['driver']) = self.services.create_sub_workflow('component_a',
                                                                                             self.services.get_config_param('COMPONENT_A_CONF'),
                                                                                             {'PWD' : self.services.get_config_param('PWD'),
                                                                                              'SIM_NAME' : self.services.get_config_param('INPUT_DIR'),
                                                                                              'COMPONENT_A_NAMELIST_INPUT' : self.services.get_config_param('COMPONENT_A_NAMELIST_INPUT'),
                                                                                              'LOG_FILE' : 'log.component_a.warning'})

        print 'creating first sub-workflow DONE!'
        ## CREATE SUB-WORKFLOW FOR MULTIPLE (num_children) FTRIDYN-XOLOTL RUNS, IN PARALLEL ##

        for i in range(0, num_children):
            child_comp = 'ftx_{}'.format(i)
            
            keys['LOG_FILE'] = 'log.{}'.format(child_comp)
            keys['SIM_NAME'] = child_comp
            keys['INPUT_DIR'] = '{}_dir'.format(child_comp)
            
            self.child_components[child_comp] = {
                                                'sim_name'  : None, #keys['SIM_NAME'],
                                                'init'      : None,
                                                'driver'    : None,
                                                'INPUT_DIR' : keys['INPUT_DIR'], #'{}_dir'.format(child_comp)
                                                'LOG_FILE'  : keys['LOG_FILE']
                                                }
        
#  Input files will be staged from this directory.

            print 'In directory', os.getcwd()
            print '\t create input directory: ', self.child_components[child_comp]['INPUT_DIR']

            if os.path.exists(self.child_components[child_comp]['INPUT_DIR']):
                shutil.rmtree(self.child_components[child_comp]['INPUT_DIR'])
            os.mkdir(self.child_components[child_comp]['INPUT_DIR'])

            #  Copy files to the created directory.
            #by hand for now --> turn this into loop over input files
            print '\t and copy input files:'
            print '\t \t parameter file 1D: ', self.SUBMIT_DIR+'/paramXolotl_1D.txt'
            print '\t \t parameter file 2D: ', self.SUBMIT_DIR+'/paramXolotl_2D.txt'
            shutil.copyfile(self.SUBMIT_DIR+'/paramXolotl_1D.txt',self.child_components[child_comp]['INPUT_DIR']+'/paramXolotl_1D.txt')
            shutil.copyfile(self.SUBMIT_DIR+'/paramXolotl_2D.txt',self.child_components[child_comp]['INPUT_DIR']+'/paramXolotl_2D.txt')

            print ' '
            print 'Create_sub_workflow with parameters:'
            print '\t child_comp = ', child_comp
            print '\t child_conf = ', child_conf
            print '\t keys = ', keys
            print '\t self.child_components[child_comp][INPUT_DIR] = ', self.child_components[child_comp]['INPUT_DIR']
            print '\n'
            (self.child_components[child_comp]['sim_name'],
             self.child_components[child_comp]['init'],
             self.child_components[child_comp]['driver']) = self.services.create_sub_workflow(child_comp, 
                                                                                              child_conf, 
                                                                                              keys, 
                                                                                              self.child_components[child_comp]['INPUT_DIR'])

        print('\n')
        print('DONE SETTING UP WORKFLOWS')
        print('\n')

#-------------------------------------------------------------------------------
#
#  parent_driver Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('parent_driver: step')

        #RUN init AND step OF COMPONENT_A (FT-GITR) SUB-WORKFLOW
        print('\n')
        print('FTRIDYN-GITR:init')
        print('\n')

        self.async_queue['component_a:driver:init'] = self.services.call(self.nested_components['component_a']['driver'], 'init', timeStamp)
        print('\n')
        print('FTRIDYN-GITR:step')
        print('\n')

        del self.async_queue['component_a:driver:init']
        self.async_queue['component_a:driver:step'] = self.services.call(self.nested_components['component_a']['driver'], 'step', 0.0)        

        ## THIS IS JUST A TEST AT READING/MODIFYING/USING GITR's OUTPUT

        #for i in range(0, num_children):
        for child_comp, child in self.child_components.items():
            i=int(filter(str.isdigit, child_comp))
            print 'index of ', child_comp, ' is ', i
            sedInFile=self.SUBMIT_DIR+'/gitrOut.txt'
            sedOutFile=self.SUBMIT_DIR+'/gitrOut_'+str(i)+'.txt'
            flux=10000+i*1000
            print 'modifying ', sedInFile, ' to write flux=',flux, ' in ', sedOutFile , '\n'
            gitrOutSedString="sed    -e 's/flux=[^ ]*/flux=%e/' < %s > %s"   % (flux, sedInFile, sedOutFile)

            subprocess.call([gitrOutSedString], shell=True)
            #shutil.copyfile(self.SUBMIT_DIR+'/gitrOut_'+str(i)+'.txt',self.child_components[child_comp]['INPUT_DIR']+'/gitrOut.txt')
            shutil.copyfile(sedOutFile,self.child_components[child_comp]['INPUT_DIR']+'/gitrOut.txt')

        #RUN init AND step OF CHILD (FT-X) SUB-WORKFLOW

        #No INIT in the FT-X workflow

        print('\n')
        print('FTRIDYN-Xolotl:init')
        print('\n')

        #  Loop over the children and all the initize component.
        for child_comp, child in self.child_components.items():
            print ' '
            print 'Call driver:init for child ', child_comp
            print '\t with dictionary ', child
            self.running_components['{}:driver:init'.format(child['sim_name'])] = self.services.call_nonblocking(child['driver'],
                                                                                                                 'init',
                                                                                                                 timeStamp,
                                                                                                                 **child) #**keys
        print ' '
        #  Loop over the children and all the initize driver component.
        for child_comp, child in self.child_components.items(): #.values():

            print ' '
            print 'for child ', child_comp, ' with sim_name ', child['sim_name'] 
            print '\t Wait for driver:init to be done'
            self.services.wait_call(self.running_components['{}:driver:init'.format(child['sim_name'])], True)
            print '\t Done waiting for driver:init'
            print '\t Call driver:step now'

            #child['LOG_FILE']='log.step.{}'.format(child_comp)
            self.running_components['{}:driver:step'.format(child['sim_name'])] = self.services.call_nonblocking(child['driver'],
                                                                                                                 'step', 
                                                                                                                 timeStamp,
                                                                                                                 **child) #**keys)
            del self.running_components['{}:driver:init'.format(child['sim_name'])]

        print ' '
#-------------------------------------------------------------------------------
#
#  parent_driver Component finalize method.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('parent_driver: finalize')

#  Wait until all the dependent components are finished.
        self.services.wait_call_list(self.running_components.values(), True)
        self.running_components = {}

#  Call finalize on all components. Need to manually call finalize on init components.

        print 'finalize FTRIDYN-GITR'
        self.async_queue['component_a:driver:finalize'] = self.services.call_nonblocking(self.nested_components['component_a']['driver'], 'finalize', timeStamp)

        print 'finalize FTRIDYN-Xolotl'
        for child_comp, child in self.child_components.items():
            self.running_components['{}:driver:finalize'.format(child['sim_name'])] = self.services.call_nonblocking(child['driver'],
                                                                                                                     'finalize', timeStamp)

            print ' '
            print '\t Child ', child_comp, ' FINALIZED!'
            print ' '


#  Wait until all the dependent components are finished.
        self.services.wait_call_list(self.running_components.values(), True)
