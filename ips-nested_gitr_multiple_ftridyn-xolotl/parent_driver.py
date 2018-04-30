#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for PARENT Driver component.
#
#-------------------------------------------------------------------------------

from component import Component
import os
import shutil

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
        self.running_components = {}
        self.child_components = {}
    
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
        
        for i in range(0, num_children):
            child_comp = 'child_{}'.format(i)
            
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
            print '\t \t gitrs output: ', self.SUBMIT_DIR+'/gitrOut_'+str(i)+'.txt'
            print '\t \t parameter file 1D: ', self.SUBMIT_DIR+'/paramXolotl_1D.txt'
            print '\t \t parameter file 2D: ', self.SUBMIT_DIR+'/paramXolotl_2D.txt'
            #print '\t \t to ', self.child_components[child_comp]['INPUT_DIR'], '\n'
            shutil.copyfile(self.SUBMIT_DIR+'/gitrOut_'+str(i)+'.txt',self.child_components[child_comp]['INPUT_DIR']+'/gitrOut.txt')
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
             self.child_components[child_comp]['driver']) = self.services.create_sub_workflow(child_comp, child_conf, keys, self.child_components[child_comp]['INPUT_DIR'])

#No INIT in the FT-X workflow
#  Loop over the children and all the initize component.
        #for child in self.child_components.values():
        #    self.running_components['{}:init:init'.format(child['sim_name'])] = self.services.call_nonblocking(child['init'],
        #                                                                                                       'init', timeStamp)
            
#  Loop over the children and all the initize component.
        for child in self.child_components.values():
#            keys = {'message' : 'Hello from {}'.format(child['sim_name'])}
        #    keys = {'message' : '10'}
        #    self.services.wait_call(self.running_components['{}:init:init'.format(child['sim_name'])], True)
            print ' '
            print 'Call driver:init for child ', child_comp
            self.running_components['{}:driver:init'.format(child['sim_name'])] = self.services.call_nonblocking(child['driver'],
                                                                                                                 'init', 
                                                                                                                 timeStamp,
                                                                                                                 **child) #**keys
        #    del self.running_components['{}:init:init'.format(child['sim_name'])]
            
#-------------------------------------------------------------------------------
#
#  parent_driver Component step method.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('parent_driver: step')
        
#  Loop over the children and all the initize driver component.
        for child_comp, child in self.child_components.items(): #.values():

            #keys = {
            #    'PWD' : self.services.get_config_param('PWD')
            #}
            
            #keys['LOG_FILE'] = 'log.{}'.format(child_comp)
            #keys['SIM_NAME'] = child_comp #self.child_components[child]['sim_name']) #child_comp
            #keys['INPUT_DIR'] = child['INPUT_DIR'] #'{}_dir'.format(child_comp)

            #print 'For child ', child_comp 
            #print '\t characterized by ', child
            #print '\t calling driver:step with keys ', keys

            print 'Wait for driver:init of child ', child_comp,' to be done'
            self.services.wait_call(self.running_components['{}:driver:init'.format(child['sim_name'])], True)

            print 'Done waiting for driver:init of child ',child_comp
            print ' '
            print 'Call driver:step now'
            print '\t for child ', child['sim_name']

            #child['LOG_FILE']='log.step.{}'.format(child_comp)
            self.running_components['{}:driver:step'.format(child['sim_name'])] = self.services.call_nonblocking(child['driver'],
                                                                                                                 'step', 
                                                                                                                 timeStamp,
                                                                                                                 **child) #**keys)

            #print ' '
            #print '\t Driver:step of child ', child_comp, ' done!'
            #print ' '
            del self.running_components['{}:driver:init'.format(child['sim_name'])]

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
        for child_comp, child in self.child_components.items():
            #self.running_components['{}:init:finalize'.format(child['sim_name'])] = self.services.call_nonblocking(child['init'],
            #                                                                                                       'finalize', timeStamp)
            self.running_components['{}:driver:finalize'.format(child['sim_name'])] = self.services.call_nonblocking(child['driver'],
                                                                                                                     'finalize', timeStamp)

            print ' '
            print '\t Child ', child_comp, ' FINALIZED!'
            print ' '


#  Wait until all the dependent components are finished.
        self.services.wait_call_list(self.running_components.values(), True)
