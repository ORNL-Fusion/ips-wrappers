#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS wrapper for PARENT Driver component.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
import os
import shutil
import subprocess
import sys
from .python_scripts_for_coupling import plasmaOut2ftxIn

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
        cwd = self.services.get_working_dir()

        print('running python ', sys.version_info)
   
        num_children = int(self.services.get_config_param('NUM_CHILDREN'))
        child_conf = self.services.get_config_param('CHILD_COMPONENT_CONF')

        self.solps_output_file=list(self.services.get_config_param('SOLPS_OUTPUT_FORMAT'))
        print('will be reading solps output from:', self.solps_output_file)
        
        self.init_time = float(self.services.get_config_param('INIT_TIME'))
        self.end_time = float(self.services.get_config_param('END_TIME'))
        self.time_step = float(self.services.get_config_param('TIME_STEP'))
        self.init_loop_n = int(self.services.get_config_param('INIT_LOOP_N'))
        
        keys = {
               'PWD' : self.services.get_config_param('PWD')
               }


        #stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        #  Sub workflows require manual setup. First a sub directory must be created.
        #  Then copying of the input files must be performed manually. The first
        #  argument of create sub workflow doesn't appear to do anything.


        #CREATE THE SOLPS SUB-WORKFLOW:

        print('CREATE THE SOLPS SUB-WORKFLOW:')
        print('\n')
        self.nested_components['component_SOLPS'] = {'sim_name': None, 'init': None, 'driver': None, 'sub_working_dir': 'component_SOLPS_init'}
        print('\t component: ', (self.nested_components['component_SOLPS']))
        print('\t working dir: ', (self.nested_components['component_SOLPS']['sub_working_dir']))

        os.mkdir(self.nested_components['component_SOLPS']['sub_working_dir'])
        ####shutil.copy2(self.services.get_config_param('COMPONENT_SOLPS_NAMELIST_INPUT'), self.child_components['component_SOLPS']['sub_working_dir'])
        print('Creating sub workflow with: ')
        print('\t name: component_SOLPS')
        print('\t config: ', self.services.get_config_param('COMPONENT_SOLPS_CONF'))
        print('\t PWD: ', self.services.get_config_param('PWD'))
        # TEST COMMENT: COMMENT OUT SOLPS SECTION FOR NOW
        #print('\t COMPONENT_SOLPS_NAMELIST_INPUT ', self.services.get_config_param('COMPONENT_SOLPS_NAMELIST_INPUT'))
        #print('\t LOG_FILE: log.component_SOLPS.warning')

        # TEST COMMENT: COMMENT OUT SOLPS SECTION FOR NOW
        #(self.nested_components['component_SOLPS']['sim_name'],
        # self.nested_components['component_SOLPS']['init'],
        # self.nested_components['component_SOLPS']['driver']) = self.services.create_sub_workflow('component_SOLPS',
        #                                                                                     self.services.get_config_param('COMPONENT_SOLPS_CONF'),
        #                                                                                     {'PWD' : self.services.get_config_param('PWD'),
        #                                                                                      'SIM_NAME' : self.services.get_config_param('INPUT_DIR'),
        #                                                                                      'COMPONENT_SOLPS_NAMELIST_INPUT' : self.services.get_config_param('COMPONENT_SOLPS_NAMELIST_INPUT'),
        #                                                                                      'LOG_FILE' : 'log.component_SOLPS.warning'})

        print('creating first sub-workflow DONE!')
        print('\n')
        print('\n')

        ## CREATE SUB-WORKFLOW FOR MULTIPLE (num_children) FTRIDYN-XOLOTL RUNS, IN PARALLEL ##

        print('CREATE SUB-WORKFLOW FOR MULTIPLE FTRIDYN-XOLOTL RUNS ')
        print('\n')
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

        cwd = self.services.get_working_dir()

        #Insert time-loop here:
        timeStamp = self.init_time                
        t_count = self.init_loop_n
        while timeStamp < self.end_time:
        
            #RUN init AND step OF COMPONENT_SOLPS SUB-WORKFLOW
            print('\n')
            print('SOLPS:init')
            print('\n')

            # TEST COMMENT: COMMENT OUT SOLPS SECTION FOR NOW
            #self.async_queue['component_SOLPS:driver:init'] = self.services.call(self.nested_components['component_SOLPS']['driver'], 'init', timeStamp)
            print('\n')
            print('SOLPS:step')
            print('\n')
            
            # TEST COMMENT: COMMENT OUT SOLPS SECTION FOR NOW
            #del self.async_queue['component_SOLPS:driver:init']

            # TEST COMMENT: COMMENT OUT SOLPS SECTION FOR NOW
            #self.async_queue['component_SOLPS:driver:step'] = self.services.call(self.nested_components['component_SOLPS']['driver'], 'step', timeStamp) 
            
            print('\n')
            print(' COMPLETED SOLPS:step')
            print('\n')

            ## change copying gitrOut.txt for:
            ## 1 - read output of SOLPS from solpsOut.txt (one per child) ; if possible as dictionary
            ## 2 - do ops to calculate/format inputs for FTX:
            ##     a) check inputs expected from SOLPS / that need special care
            ##     i.  inputs that are single values:
            ##         particle flux, heat flux
            ##     ii. inputs that are given for each species:
            ##         for now: always include information for all 4 species: simpler to handle
            ##                  and doesn't add comp time as the ftx workflow skips anything with fluxFraction=0
            ##         fluxFraction, [Ti, Te, Z] --> E_in ; B-field --> A_in
            ##     b) print all others 'as is'
            ## 3 - print ftxIn.txt
            ## 4 - figure out how to pass time of loop to FTX (in ftxIn.txt?)
            ## 5 - copy ftxIn to child's input directory

            print('\n')
            print('READ AND FORMAT INPUT FTX WORKFLOWS')
            print('\n')
            
            for child_comp, child in list(self.child_components.items()):

                #to get i, these two ways should work:
                #i=int(list(filter(str.isdigit, child_comp))[0]) #not tested for cases with more than 1 child
                i=int(''.join(list(filter(str.isdigit, child_comp))))
                print('\t index of ', child_comp, ' is ', i)
                sys.stdout.flush()
                
                ## 1 - read output of SOLPS from solpsOut.txt (one per child) - DONE?
                print('\t re-format SOLPS output into FTX input: ')
                solps_outFile = self.solps_output_file[0]+'_'+str(i)+'.'+self.solps_output_file[1] #0=filename (solpsOut) ; 1=fomat (txt) 
                ftxInFileFormat=list(self.services.get_config_param('FTX_INPUT_FORMAT'))
                ftxInFileName=ftxInFileFormat[0]+str(i)+'.'+ftxInFileFormat[1]                
                print('\t \t call plasmaOut2ftxIn with: ')
                print('\t \t \t TEST: input file: ', solps_outFile)
                print('\t \t \t TEST: output file: ', ftxInFileName)
                plasmaOut2ftxIn.plasmaOut2ftxIn(plasmaOutFile=self.INPUT_DIR+'/'+solps_outFile, ftxInFile=cwd+'/'+ftxInFileName)


                # 4 - add timeStap to its own file
                #  For now keep this section in driver 
                #  Move to its own function is time parameter operations become more complex
                print('\t write time input file:')
                timeFileName=self.services.get_config_param('TIME_FILE_NAME')
                print('\t \t file name: ', timeFileName)
                timeFile=open(timeFileName, "w")
                print('\t \t INIT_TIME : ', timeStamp)
                print('\t \t END_TIME : ', timeStamp + self.time_step)
                print('\t \t LOOP_TIME_STEP : ', self.time_step)
                print('\t \t LOOP_N : ', t_count)
                timeFile.write('INIT_TIME={0}\n'.format(timeStamp))
                timeFile.write('END_TIME={0}\n'.format(timeStamp + self.time_step))
                timeFile.write('LOOP_TIME_STEP={0}\n'.format(self.time_step))
                timeFile.write('LOOP_N={0}\n'.format(t_count))
                print(' ')

                #set start mode based on loop number ; change to driver's start mode is we implement restarting capabilities
                if (t_count == 0):
                    print('\t \t START_MODE = INIT')
                    timeFile.write('START_MODE=INIT\n')
                else:
                    print('\t \t START_MODE = RESTART')
                    timeFile.write('START_MODE=RESTART\n')
                timeFile.close()
                print(' ')

                
                # 5 - copy ftxIn and timeFile to ftx input directory
                print('\t copying input to FTX file')
                print('\t \t from: ', cwd+'/'+ftxInFileName)
                print('\t \t to: ', self.child_components[child_comp]['INPUT_DIR']+'/ftxInput.txt')
                shutil.copyfile(cwd+'/'+ftxInFileName,self.child_components[child_comp]['INPUT_DIR']+'/ftxInput.txt')
                print('\t copying time parameter ')
                print('\t \t from: ', timeFileName)
                print('\t \t to: ', self.child_components[child_comp]['INPUT_DIR']+'/'+timeFileName)
                shutil.copyfile(timeFileName,self.child_components[child_comp]['INPUT_DIR']+'/'+timeFileName)
                print('\t copying output file of SOLPS')
                print('\t \t from: ',self.INPUT_DIR+'/'+solps_outFile)
                print('\t \t to: ', self.child_components[child_comp]['INPUT_DIR']+'/solpsOut.txt')
                shutil.copyfile(self.INPUT_DIR+'/'+solps_outFile,self.child_components[child_comp]['INPUT_DIR']+'/solpsOut.txt')
                
                print('\n')
                sys.stdout.flush()

                #save ftxIn and solpsOut for each loop (with time-stamp)
                shutil.copyfile(self.INPUT_DIR+'/'+solps_outFile, solps_outFile+'_t'+str(timeStamp))
                shutil.copyfile(cwd+'/'+ftxInFileName,cwd+'/'+ftxInFileName+'_t'+str(timeStamp))
                shutil.copyfile(timeFileName,timeFileName+'_t'+str(timeStamp))
                print('\t to save input files for each time-loop, copied:')
                print('\t \t ', self.INPUT_DIR+'/'+solps_outFile, ' as ', solps_outFile+'_t'+str(timeStamp))
                print('\t \t ', cwd+'/'+ftxInFileName, ' as ', ftxInFileName+'_t'+str(timeStamp))
                print('\t \t ', timeFileName , ' as ', timeFileName+'_t'+str(timeStamp))

                
                #update plasma state file:
                self.services.update_state()
                
            ## END OF solpsOut --> ftxIn and time-parameters
            
            #RUN init AND step OF CHILD (FT-X) SUB-WORKFLOW
            
            print('\n')
            print('FTRIDYN-Xolotl:init')
            print('\n')
            
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
            sys.stdout.flush()
            
            #  Loop over the children and all the initize driver component.
            for child_comp, child in list(self.child_components.items()): #.values():
                print(' ')
                print('for child ', child_comp, ' with sim_name ', child['sim_name']) 
                print('\t Wait for driver:init to be done')
                self.services.wait_call(self.running_components['{}:driver:init'.format(child['sim_name'])], True)
                print('\t Done waiting for driver:init')
                print('\t Call driver:step now')
                sys.stdout.flush()
                
                #child['LOG_FILE']='log.step.{}'.format(child_comp)
                self.running_components['{}:driver:step'.format(child['sim_name'])] = self.services.call_nonblocking(child['driver'],
                                                                                                                     'step', 
                                                                                                                     timeStamp,
                                                                                                                     **child) #**keys)
                sys.stdout.flush()
                
            for child_comp, child in list(self.child_components.items()):
                print(' ')
                print('for child ', child_comp, ' with sim_name ', child['sim_name'])
                print('\t Wait for driver:step to be done')
                self.services.wait_call(self.running_components['{}:driver:step'.format(child['sim_name'])], True)
                print('\t Done waiting for driver:step')
                sys.stdout.flush()

            for child_comp, child in list(self.child_components.items()):
                del self.running_components['{}:driver:init'.format(child['sim_name'])]
                del self.running_components['{}:driver:step'.format(child['sim_name'])] 
                sys.stdout.flush()

            timeStamp+=self.time_step
            t_count+=1
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
        # TEST COMMENT: COMMENT OUT SOLPS SECTION FOR NOW
        #self.async_queue['component_SOLPS:driver:finalize'] = self.services.call_nonblocking(self.nested_components['component_SOLPS']['driver'], 'finalize', timeStamp)

        print('finalize FTRIDYN-Xolotl')
        for child_comp, child in list(self.child_components.items()):
            self.running_components['{}:driver:finalize'.format(child['sim_name'])] = self.services.call_nonblocking(child['driver'],
                                                                                                                     'finalize', timeStamp)

            print(' ')
            print('\t Child ', child_comp, ' FINALIZED!')
            print(' ')


#  Wait until all the dependent components are finished.
        self.services.wait_call_list(list(self.running_components.values()), True)
