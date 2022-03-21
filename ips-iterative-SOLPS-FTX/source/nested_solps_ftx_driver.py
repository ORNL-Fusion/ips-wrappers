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
        #dictionaries to store the SOLPS and FTX workflows and parametrization (config file, etc.) 
        self.running_components = {}
        self.child_components = {}
        self.async_queue = {}
        self.nested_components = {} 

#-------------------------------------------------------------------------------
#
#  parent_driver Component init method.
#
#-------------------------------------------------------------------------------


    ###########  NOTES CHANGING TO SOLPS-FTX   ###########
    ##  no big changes should be needed in init method  ## 
    ## only creates the components & stage input files  ##

    def init(self, timeStamp=0.0):
        print(' ')
        print('parent_driver: init')
 
        import sys
        print('running python ', sys.version_info)
   
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


        #CREATE THE SOLPS SUB-WORKFLOW: 
        
        self.nested_components['component_SOLPS'] = {'sim_name': None, 'init': None, 'driver': None, 'sub_working_dir': 'component_SOLPS_init'}
        print((self.nested_components['component_SOLPS']))
        print((self.nested_components['component_SOLPS']['sub_working_dir']))

        os.mkdir(self.nested_components['component_SOLPS']['sub_working_dir'])
        #shutil.copy2(self.services.get_config_param('COMPONENT_SOLPS_NAMELIST_INPUT'), self.child_components['component_SOLPS']['sub_working_dir'])
        print('Creating sub workflo with: ')
        print('\t name: component_SOLPS')
        print('\t config: ', self.services.get_config_param('COMPONENT_SOLPS_CONF'))
        print('\t PWD: ', self.services.get_config_param('PWD'))
        print('\t COMPONENT_SOLPS_NAMELIST_INPUT ', self.services.get_config_param('COMPONENT_SOLPS_NAMELIST_INPUT'))
        print('\t LOG_FILE: log.component_SOLPS.warning')

        (self.nested_components['component_SOLPS']['sim_name'],
         self.nested_components['component_SOLPS']['init'],
         self.nested_components['component_SOLPS']['driver']) = self.services.create_sub_workflow('component_SOLPS',
                                                                                             self.services.get_config_param('COMPONENT_SOLPS_CONF'),
                                                                                             {'PWD' : self.services.get_config_param('PWD'),
                                                                                              'SIM_NAME' : self.services.get_config_param('INPUT_DIR'),
                                                                                              'COMPONENT_SOLPS_NAMELIST_INPUT' : self.services.get_config_param('COMPONENT_SOLPS_NAMELIST_INPUT'),
                                                                                              'LOG_FILE' : 'log.component_SOLPS.warning'})

        print('creating first sub-workflow DONE!')


        ## CREATE SUB-WORKFLOW FOR MULTIPLE (num_children) FTRIDYN-XOLOTL RUNS, IN PARALLEL ##

        for i in range(0, num_children):
            child_comp = 'ftx_{}'.format(i)
            
            keys['LOG_FILE'] = 'log.{}'.format(child_comp)
            keys['SIM_NAME'] = child_comp
            keys['INPUT_DIR'] = '{}_dir'.format(child_comp)
            keys['INPUT_FILES'] = self.services.get_config_param('INPUT_FTX')
            
            self.child_components[child_comp] = {
                                                'sim_name'  : None, #keys['SIM_NAME'],
                                                'init'      : None,
                                                'driver'    : None,
                                                'INPUT_DIR' : keys['INPUT_DIR'], #'{}_dir'.format(child_comp)
                                                'LOG_FILE'  : keys['LOG_FILE']
                                                }
        
#  Input files will be staged from this directory.

            print('In directory', os.getcwd())
            print('\t create input directory: ', self.child_components[child_comp]['INPUT_DIR'])

            if os.path.exists(self.child_components[child_comp]['INPUT_DIR']):
                shutil.rmtree(self.child_components[child_comp]['INPUT_DIR'])
            os.mkdir(self.child_components[child_comp]['INPUT_DIR'])

            #  Copy files to the created directory.
            input_file_list=keys['INPUT_FILES'].split()
            print('\t and copy input files:' , input_file_list) #self.child_components[child_comp]['INPUT_FILES']
            for input_file in input_file_list: #self.child_components[child_comp]['INPUT_FILES'])):
                print('\t \t', input_file , ' from ', self.SUBMIT_DIR, ' to ', self.child_components[child_comp]['INPUT_DIR'])
                shutil.copyfile(self.SUBMIT_DIR+'/'+input_file,self.child_components[child_comp]['INPUT_DIR']+'/'+input_file)
            print('\t ...done copying input files')

            print(' ')
            print('Create_sub_workflow with parameters:')
            print('\t child_comp = ', child_comp)
            print('\t child_conf = ', child_conf)
            print('\t keys = ', keys)
            print('\t self.child_components[child_comp][INPUT_DIR] = ', self.child_components[child_comp]['INPUT_DIR'])
            print('\n')
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

   #############  NOTES CHANGING TO SOLPS-FTX  #############
   ##  step-method is now a loop over time for SOLPS+FTX  ##
   ##  SOLPS-init within loop for now. Move out as needed ##
   ##  e.g.: if solps:init models inter-ELM steady-state  ## 
   
    def step(self, timeStamp=0.0):
        print('parent_driver: step')

        ## TO-DO CHANGING SOLPS-FTX:
        ## INSERT LOOP OVER TIME HERE
        ## INIT TIME, END TIME, STEP TIME FROM CONFIG        
        
        #RUN init AND step OF COMPONENT_SOLPS SUB-WORKFLOW
        print('\n')
        print('SOLPS:init')
        print('\n')

        self.async_queue['component_SOLPS:driver:init'] = self.services.call(self.nested_components['component_SOLPS']['driver'], 'init', timeStamp)
        print('\n')
        print('SOLPS:step')
        print('\n')

        del self.async_queue['component_SOLPS:driver:init']

        ## TO-DO CHANGING TO SOLPS-FTX
        ## change 0.0 in call to timeStamp, as/if needed 
        self.async_queue['component_SOLPS:driver:step'] = self.services.call(self.nested_components['component_SOLPS']['driver'], 'step', 0.0)        

        print('\n')
        print(' COMPLETED SOLPS:step')
        print('\n')

        ## TO-DO CHANGING TO SOLPS-FTX
        ## change copying gitrOut.txt for:
        ## 1 - read output of SOLPS from solpsOut.txt (one per child) ; if possible as dictionary
        ## 2 - do ops to calculate/format inputs for FTX:
        ##        e.g., [Ti, Te, Z] --> E_in ; B-field --> A_in
        ## 3 - print ftxIn.txt
        ##        include reformating entries
        ##        e.g., plasmaSpecies = D C --> plasmaSpecies = He W D T C
        ## 4 - figure out how to pass time of loop to FTX (in ftxIn.txt?)
        
        for child_comp, child in list(self.child_components.items()):
            i=int(list(filter(str.isdigit, child_comp)))
            print('index of ', child_comp, ' is ', i)
            sedOutFile=self.SUBMIT_DIR+'/gitrOut'+str(i)+'.txt'
            print('copying gitr file from ', sedOutFile, 'to ', self.child_components[child_comp]['INPUT_DIR']+'/gitrOut.txt')
            shutil.copyfile(sedOutFile,self.child_components[child_comp]['INPUT_DIR']+'/gitrOut.txt')

        ## END OF TO-DO CHANGING TO SOLPS-FTX

            
        #RUN init AND step OF CHILD (FT-X) SUB-WORKFLOW

        print('\n')
        print('FTRIDYN-Xolotl:init')
        print('\n')

        ## TO-DO CHANGING TO SOLPS-FTX
        ## make sure current (simulated) time is passed to FTX
        ## make sure end-time of this loop/run is passed (see change 4 above)
        ##      preferably to the FTX driver
        ##      at the very least to Xolotl
        
        #  Loop over the children and all the initize component.
        for child_comp, child in list(self.child_components.items()):
            print(' ')
            print('Call driver:init for child ', child_comp)
            print('\t with dictionary ', child)
            self.running_components['{}:driver:init'.format(child['sim_name'])] = self.services.call_nonblocking(child['driver'],
                                                                                                                 'init',
                                                                                                                 timeStamp,
                                                                                                                 **child) #**keys
        print(' ')
        #  Loop over the children and all the initize driver component.
        for child_comp, child in list(self.child_components.items()): #.values():

            print(' ')
            print('for child ', child_comp, ' with sim_name ', child['sim_name']) 
            print('\t Wait for driver:init to be done')
            self.services.wait_call(self.running_components['{}:driver:init'.format(child['sim_name'])], True)
            print('\t Done waiting for driver:init')
            print('\t Call driver:step now')

            #child['LOG_FILE']='log.step.{}'.format(child_comp)
            self.running_components['{}:driver:step'.format(child['sim_name'])] = self.services.call_nonblocking(child['driver'],
                                                                                                                 'step', 
                                                                                                                 timeStamp,
                                                                                                                 **child) #**keys)
            del self.running_components['{}:driver:init'.format(child['sim_name'])]


             ## TO-DO CHANGING TO SOLPS-FTX ##
             ## END LOOP OVER TIME HERE


        print(' ')
#-------------------------------------------------------------------------------
#
#  parent_driver Component finalize method.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('parent_driver: finalize')

#  Wait until all the dependent components are finished.
        self.services.wait_call_list(list(self.running_components.values()), True)
        self.running_components = {}

#  Call finalize on all components. Need to manually call finalize on init components.

        print('finalize SOLPS')
        self.async_queue['component_SOLPS:driver:finalize'] = self.services.call_nonblocking(self.nested_components['component_SOLPS']['driver'], 'finalize', timeStamp)

        print('finalize FTRIDYN-Xolotl')
        for child_comp, child in list(self.child_components.items()):
            self.running_components['{}:driver:finalize'.format(child['sim_name'])] = self.services.call_nonblocking(child['driver'],
                                                                                                                     'finalize', timeStamp)

            print(' ')
            print('\t Child ', child_comp, ' FINALIZED!')
            print(' ')


#  Wait until all the dependent components are finished.
        self.services.wait_call_list(list(self.running_components.values()), True)
