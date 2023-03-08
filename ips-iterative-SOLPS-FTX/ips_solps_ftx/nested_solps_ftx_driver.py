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
from ips_solps_ftx.python_scripts_for_coupling import plasmaOut2ftxIn
import pickle

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
        self.ftx_components = {}
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
        self.print_test = self.services.get_config_param('PRINT_TEST')

        #test giving explicit wrapper path in modernized FTX workflow
        if self.print_test:
            try:
                self.SCRIPT
                if self.SCRIPT == "":
                    print('no explicit script path provided. use module loaded in environment')
                else:
                    print('using explicit path to wrapper')
                    print(self.SCRIPT)
            except Exception as e:
                print(e)
                print('no script variable defined. use module loaded in environment')
            print(' ')

        
        if (self.print_test):
            print('running python ', sys.version_info)
   
        num_ftx = int(self.services.get_config_param('NUM_SUBWF'))
        ftx_conf = self.services.get_config_param('SUBWF_COMPONENT_CONF')
        solps_conf = self.services.get_config_param('SOLPS_COMPONENT_CONF')
        
        self.solps_output_file=list(self.services.get_config_param('SOLPS_OUTPUT_FORMAT'))

        if (self.print_test):
            print('will be reading solps output from:', self.solps_output_file)
        
        self.init_time = float(self.services.get_config_param('INIT_TIME'))
        self.end_time = float(self.services.get_config_param('END_TIME'))
        self.time_step = float(self.services.get_config_param('TIME_STEP'))
        self.init_loop_n = int(self.services.get_config_param('INIT_LOOP_N'))

        ftx_keys = {} # test without -- not sure it's used'PWD' : self.services.get_config_param('PWD')


        #stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        #  Sub workflows require manual setup. First a sub directory must be created.
        #  Then copying of the input files must be performed manually. The first
        #  argument of create sub workflow doesn't appear to do anything.


        #CREATE THE SOLPS SUB-WORKFLOW:
        print('\n')
        print('CREATE THE SOLPS SUB-WORKFLOW:')
        print('\n')
        self.nested_components['component_SOLPS'] = {'sim_name': None,
                                                     'init': None,
                                                     'driver': None,
                                                     'INPUT_DIR' : 'solps_dir',
                                                     'LOG_FILE' : 'log.component_SOLPS.warning'
                                                     }
        if (self.print_test):
            print('Defined dictionary : ')
            print('\t', (self.nested_components['component_SOLPS']))

        print('Creating sub workflow with: ')
        print('\t name: component_SOLPS')
        print('\t config: ', solps_conf)
        print('\t keys: SIM_NAME : solps_iter', 'LOG_FILE : log.component_SOLPS.warning')
        print('\t INPUT DIR ', self.services.get_config_param('SOLPS_INPUT_DIR'))

        (self.nested_components['component_SOLPS']['sim_name'],
         self.nested_components['component_SOLPS']['init'],
         self.nested_components['component_SOLPS']['driver']) = self.services.create_sub_workflow('component_SOLPS',
                                                                                                  solps_conf, 
                                                                                                  {'SIM_NAME' : 'solps_iter',
                                                                                                   'LOG_FILE' : 'log.component_SOLPS.warning'},
                                                                                                  self.services.get_config_param('SOLPS_INPUT_DIR'))

        print('creating SOLPS sub-workflow DONE!')
        if (self.print_test):
            print('\t SOLPS component is:') #TEST PRINT
            print('\t', self.nested_components['component_SOLPS'])
        print('\n')

        ## CREATE SUB-WORKFLOW FOR MULTIPLE (num_ftx) FTRIDYN-XOLOTL RUNS, IN PARALLEL ##
        print('CREATE SUB-WORKFLOW FOR MULTIPLE FTRIDYN-XOLOTL RUNS ')
        print('\n')
        for i in range(0, num_ftx):
            ftx_comp = 'ftx_{}'.format(i)
            
            ftx_keys['LOG_FILE'] = 'log.{}'.format(ftx_comp)
            ftx_keys['SIM_NAME'] = ftx_comp
            ftx_keys['INPUT_DIR'] = 'input_{}'.format(ftx_comp)
            ftx_keys['INPUT_FILES'] = self.services.get_config_param('INPUT_FTX')
            
            self.ftx_components[ftx_comp] = {
                                                'sim_name'  : None, #keys['SIM_NAME'],
                                                'init'      : None,
                                                'driver'    : None,
                                                'INPUT_DIR' : ftx_keys['INPUT_DIR'], #'{}_dir'.format(child_comp)
                                                'LOG_FILE'  : ftx_keys['LOG_FILE']
                                                }
            if (self.print_test):
                print('Defined dictionary : ')
                print('\t', (self.ftx_components[ftx_comp]))
                
            #  Input files will be staged from this directory.

            print('In directory', os.getcwd())
            print('\t create input directory: ', self.ftx_components[ftx_comp]['INPUT_DIR'])

            if os.path.exists(self.ftx_components[ftx_comp]['INPUT_DIR']):
                shutil.rmtree(self.ftx_components[ftx_comp]['INPUT_DIR'])
            os.mkdir(self.ftx_components[ftx_comp]['INPUT_DIR'])

            #  Copy files to the created directory.
            input_file_list=ftx_keys['INPUT_FILES'].split()
            print('\t and copy input files:' , input_file_list) #self.child_components[child_comp]['INPUT_FILES']
            for input_file in input_file_list: #self.child_components[child_comp]['INPUT_FILES'])):
                print('\t \t', input_file , ' from ', self.SUBMIT_DIR, ' to ', self.ftx_components[ftx_comp]['INPUT_DIR'])
                shutil.copyfile(self.SUBMIT_DIR+'/'+input_file,self.ftx_components[ftx_comp]['INPUT_DIR']+'/'+input_file)
            print('\t ...done copying input files')

            print(' ')
            print('Create_sub_workflow with parameters:')
            print('\t ftx_comp = ', ftx_comp)
            print('\t ftx_conf = ', ftx_conf)
            print('\t keys = ', ftx_keys)
            print('\t self.ftx_components[ftx_comp][INPUT_DIR] = ', self.ftx_components[ftx_comp]['INPUT_DIR'])
            print('\n')
            (self.ftx_components[ftx_comp]['sim_name'],
             self.ftx_components[ftx_comp]['init'],
             self.ftx_components[ftx_comp]['driver']) = self.services.create_sub_workflow(ftx_comp, 
                                                                                              ftx_conf, 
                                                                                              ftx_keys, 
                                                                                              self.ftx_components[ftx_comp]['INPUT_DIR'])
            print('creating FTX sub-workflow DONE!')
            if (self.print_test):
                print('FTX component is:') #TEST PRINT
            print('\t sim_name : ',self.ftx_components[ftx_comp]['sim_name'],'init : ',self.ftx_components[ftx_comp]['init'],'driver : ',self.ftx_components[ftx_comp]['driver']) #TEST PRINT
            print('\t keys : ', ftx_keys, 'INPUT_DIR : ', self.ftx_components[ftx_comp]['INPUT_DIR'])
            print('\n')
            
        #print('\n')
        print('DONE SETTING UP WORKFLOWS')
        print('\n')

        self.services.stage_state()
        self.services.update_state()
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

        #if we ever define specific outFile / logFile, do so here:
        
        #Insert time-loop here:
        timeStamp = self.init_time                
        t_count = self.init_loop_n
        while timeStamp < self.end_time:

            print('\n')
            print('-------------------')
            print('-------------------')
            print('run SOLPS')
            print('-------------------')
            print('-------------------')
            print('\n')


            #test activating the SOLPS environment
            #subprocess.call("conda deactivate", shell=True)
            #subprocess.call("conda activate /global/common/software/atom/cori/cesol_conda/v0.1", shell=True)

            
            #RUN init AND driver OF COMPONENT_SOLPS SUB-WORKFLOW
            print('\t \n')
            print('\t SOLPS:init')
            print('\t \n')

            #solps init:
            print('run SOLPS init:init')
            self.async_queue['component_SOLPS:init:init'] = self.services.call(self.nested_components['component_SOLPS']['init'], 'init', timeStamp)                    
            print('init:init done \n')
            #del self.async_queue['component_SOLPS:init:init']
            self.services.stage_state()

            print('run SOLPS init:step')
            self.async_queue['component_SOLPS:init:step'] = self.services.call(self.nested_components['component_SOLPS']['init'], 'step', timeStamp)
            print('init:step done \n')
            #del self.async_queue['component_SOLPS:init:step']
            self.services.stage_state()

            print('run SOLPS init:finalize')
            self.async_queue['component_SOLPS:init:finalize'] = self.services.call(self.nested_components['component_SOLPS']['init'], 'step', timeStamp)
            print('init:step done \n')
            #del self.async_queue['component_SOLPS:init:finalize']
            self.services.stage_state()

            
            print('\t \n')
            print('\t SOLPS:driver') #-solver-b2plot')
            print('\t \n')
            
            print('run SOLPS driver:init')
            self.async_queue['component_SOLPS:driver:init'] = self.services.call(self.nested_components['component_SOLPS']['driver'], 'init', timeStamp)
            print('driver:init done \n')
            #del self.async_queue['component_SOLPS:driver:init']  
            self.services.stage_state()

            #print('\t run SOLPS solver:init')
            #self.async_queue['component_SOLPS:solver:init'] = self.services.call(self.nested_components['component_SOLPS']['solver'], 'init', timeStamp)
            #print('\t solver:init done')
            #del self.async_queue['component_SOLPS:solver:init']
            #self.services.stage_state()

            #print('\t run SOLPS b2plot:init')
            #self.async_queue['component_SOLPS:b2plot:init'] = self.services.call(self.nested_components['component_SOLPS']['b2plot'], 'init', timeStamp)
            #print('\t b2plot:init done')
            #del self.async_queue['component_SOLPS:b2plot:init']
            #self.services.stage_state()

            
            
            print('run SOLPS driver:step')
            self.async_queue['component_SOLPS:driver:step'] = self.services.call(self.nested_components['component_SOLPS']['driver'], 'step', timeStamp) 
            #del self.async_queue['component_SOLPS:driver:step']
            print('driver:step done \n')
            self.services.stage_state()

            #print('\t run SOLPS solver:step')
            #self.async_queue['component_SOLPS:solver:step'] = self.services.call(self.nested_components['component_SOLPS']['solver'], 'step', timeStamp)
            #print('\t solver:step done')
            #del self.async_queue['component_SOLPS:solver:step']
            #self.services.stage_state()

            #print('\t run SOLPS b2plot:step')
            #self.async_queue['component_SOLPS:b2plot:step'] = self.services.call(self.nested_components['component_SOLPS']['b2plot'], 'step', timeStamp)
            #print('\t b2plot:step done')
            #del self.async_queue['component_SOLPS:b2plot:step']
            #self.services.stage_state()
            
            print('run SOLPS driver:finalize')
            self.async_queue['component_SOLPS:driver:finalize'] = self.services.call(self.nested_components['component_SOLPS']['driver'], 'finalize', timeStamp)
            #del self.async_queue['component_SOLPS:driver:finalize']
            print('driver:finalize done \n')
            self.services.stage_state()

            #print('\t run SOLPS solver:finalize')
            #self.async_queue['component_SOLPS:solver:finalize'] = self.services.call(self.nested_components['component_SOLPS']['solver'], 'finalize', timeStamp)
            #print('\t solver:finalize done')
            #del self.async_queue['component_SOLPS:solver:finalize']
            #self.services.stage_state()

            #print('\t run SOLPS b2plot:finalize')
            #self.async_queue['component_SOLPS:b2plot:finalize'] = self.services.call(self.nested_components['component_SOLPS']['b2plot'], 'finalize', timeStamp)
            #print('\t b2plot:finalize done')
            #del self.async_queue['component_SOLPS:b2plot:finalize']
            #self.services.stage_state()

            
            print('\n')
            print('COMPLETED SOLPS')
            #print('\n')

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
            print('--------------------')
            print('--------------------')
            print('SOLPS out --> FTX in')
            print('--------------------')
            print('--------------------')
            print('\n')
            
            #print('\n')
            print('READ AND FORMAT INPUT FTX WORKFLOWS')
            print('\n')
            
            for ftx_comp, ftx in list(self.ftx_components.items()):

                #to get i, these two ways should work:
                #i=int(list(filter(str.isdigit, child_comp))[0]) #not tested for cases with more than 1 child
                i=int(''.join(list(filter(str.isdigit, ftx_comp))))
                print('index of ', ftx_comp, ' is ', i)
                sys.stdout.flush()
                
                ## 1 - read output of SOLPS from solpsOut.txt (one per child) - DONE?
                print('re-format SOLPS output into FTX input: ')
                #these three are used later as well, so define here
                ftxInFileFormat=list(self.services.get_config_param('FTX_INPUT_FORMAT'))
                ftxInFileName=ftxInFileFormat[0]+str(i)+'.'+ftxInFileFormat[1]
                solps_outFile=self.solps_output_file[0]+'_'+str(i)+'.'+self.solps_output_file[1] #0=filename (solpsOut) ; 1=fomat (txt)
                p2ftx_log=cwd+'/log.plasmaOut2ftxIn'+'_t'+str(timeStamp)

                self.plasmaOut2ftxIn={}
                self.plasmaOut2ftxIn['plasmaOutFile'] = self.INPUT_DIR+'/'+solps_outFile
                self.plasmaOut2ftxIn['ftxInFile'] = cwd+'/'+ftxInFileName
                self.plasmaOut2ftxIn['logFile']= p2ftx_log #or logFile (=None) 
                self.plasmaOut2ftxIn['print_test'] = self.print_test

                #launch plasmaOut2ftxIn instead of calling function:            
                #plasmaOut2ftxIn.plasmaOut2ftxIn(plasmaOutFile=self.INPUT_DIR+'/'+solps_outFile, ftxInFile=cwd+'/'+ftxInFileName, print_test=self.print_test, logFile=plasmaOut2ftxIn_log)
                pkl_p2ftx_file=cwd+'/plasmaOut2ftxIn.pkl'  ## TO-DO: define name in config file, here give abs path
                pickle.dump(self.plasmaOut2ftxIn, open(pkl_p2ftx_file, "wb" ) )

                if (self.print_test):
                    print('call plasmaOut2ftxIn with: ')
                    print('\t input file: ', self.plasmaOut2ftxIn['plasmaOutFile'])
                    print('\t output file: ', self.plasmaOut2ftxIn['ftxInFile'])
                    print('\t print_test: ', self.plasmaOut2ftxIn['print_test'])
                    print('\t output log: ', self.plasmaOut2ftxIn['logFile'])
                
                sys.stdout.flush()
                try:
                    self.TRANSLATE_P2FTX
                    plasma2ftx_script=self.TRANSLATE_P2FTX
                    print('\t Launch user-defined python script :  ')
                    print('\t', plasma2ftx_script)
                    
                except: #DEFAULT PATH:                                                                                                                                                                        
                    plasma2ftx_script = 'plasmaOut2ftxIn.py'
                    print('\t Launch default python script: ')
                    print('\t ',plasma2ftx_script)
                    
                sys.stdout.flush()
                if self.print_test:
                    print('\t launch task : ', plasma2ftx_script)
                    print('\t with pkl file: ', pkl_p2ftx_file)
                    print('\t logFile : ', p2ftx_log)
                    print('\t self.services.get_working_dir() = ', self.services.get_working_dir())
                
                task_id_p2ftx = self.services.launch_task(1,self.services.get_working_dir(),
                                                          plasma2ftx_script, logFile=p2ftx_log)
                ret_val_p2ftx = self.services.wait_task(task_id_p2ftx)
                if self.print_test:
                    print('after running ', plasma2ftx_script)
                    print('ret_val_p2ftx = ', ret_val_p2ftx)
                sys.stdout.flush()
                #self.plasmaOut2ftxIn.clear()
                
                
                # 4 - add timeStap to its own file
                #  For now keep this section in driver 
                #  Move to its own function is time parameter operations become more complex
                print('write time input file:')
                timeFileName=self.services.get_config_param('TIME_FILE_NAME')
                print('\t file name: ', timeFileName)
                timeFile=open(timeFileName, "w")
                print('\t INIT_TIME : ', timeStamp)
                print('\t END_TIME : ', timeStamp + self.time_step)
                print('\t LOOP_TIME_STEP : ', self.time_step)
                print('\t LOOP_N : ', t_count)
                timeFile.write('INIT_TIME={0}\n'.format(timeStamp))
                timeFile.write('END_TIME={0}\n'.format(timeStamp + self.time_step))
                timeFile.write('LOOP_TIME_STEP={0}\n'.format(self.time_step))
                timeFile.write('LOOP_N={0}\n'.format(t_count))
                #print(' ')
                sys.stdout.flush()
                
                #set start mode based on loop number ; change to driver's start mode is we implement restarting capabilities
                if (t_count == 0):
                    print('\t START_MODE = INIT')
                    timeFile.write('START_MODE=INIT\n')
                else:
                    print('\t START_MODE = RESTART')
                    timeFile.write('START_MODE=RESTART\n')
                timeFile.close()
                print(' ')
                sys.stdout.flush()
                
                # 5 - copy ftxIn and timeFile to ftx input directory
                print('copying input to FTX file')
                print('\t from: ', cwd+'/'+ftxInFileName)
                print('\t to: ', self.ftx_components[ftx_comp]['INPUT_DIR']+'/ftxInput.txt')
                shutil.copyfile(cwd+'/'+ftxInFileName,self.ftx_components[ftx_comp]['INPUT_DIR']+'/ftxInput.txt')
                print('copying time parameter ')
                print('\t from: ', timeFileName)
                print('\t to: ', self.ftx_components[ftx_comp]['INPUT_DIR']+'/'+timeFileName)
                shutil.copyfile(timeFileName,self.ftx_components[ftx_comp]['INPUT_DIR']+'/'+timeFileName)
                print('copying output file of SOLPS')
                print('\t from: ',self.INPUT_DIR+'/'+solps_outFile)
                print('\t to: ', self.ftx_components[ftx_comp]['INPUT_DIR']+'/solpsOut.txt')
                shutil.copyfile(self.INPUT_DIR+'/'+solps_outFile,self.ftx_components[ftx_comp]['INPUT_DIR']+'/solpsOut.txt')
                
                print('\n')
                sys.stdout.flush()

                #save ftxIn and solpsOut for each loop (with time-stamp)
                shutil.copyfile(self.INPUT_DIR+'/'+solps_outFile, solps_outFile+'_t'+str(timeStamp))
                shutil.copyfile(cwd+'/'+ftxInFileName,cwd+'/'+ftxInFileName+'_t'+str(timeStamp))
                shutil.copyfile(timeFileName,timeFileName+'_t'+str(timeStamp))
                print('to save input files for each time-loop, copied:')
                print('\t ', self.INPUT_DIR+'/'+solps_outFile, ' as ', solps_outFile+'_t'+str(timeStamp))
                print('\t ', cwd+'/'+ftxInFileName, ' as ', ftxInFileName+'_t'+str(timeStamp))
                print('\t ', timeFileName , ' as ', timeFileName+'_t'+str(timeStamp))

                sys.stdout.flush()
                #update plasma state file:
                self.services.update_state()
                
            ## END OF solpsOut --> ftxIn and time-parameters

            print('\n')
            print('-------------------')
            print('-------------------')
            print('run FTRIDYN-Xolotl')
            print('-------------------')
            print('-------------------')
            
            #RUN init AND step OF CHILD (FT-X) SUB-WORKFLOW            
            print('\n')
            print('FTRIDYN-Xolotl:init')
            print('\n')
            
            #  Loop over the children and all the initize component.
            for ftx_comp, ftx in list(self.ftx_components.items()):
                print(' ')
                print('Call driver:init for subworkflow ', ftx_comp)
                print('with dictionary ', ftx)
                self.running_components['{}:driver:init'.format(ftx['sim_name'])] = self.services.call_nonblocking(ftx['driver'],
                                                                                                                     'init',
                                                                                                                     timeStamp,
                                                                                                                     **ftx) #**keys
            print(' ')
            sys.stdout.flush()
            self.services.stage_state()
            
            #  Loop over the children and all the initialize driver component.
            print('\n')
            print('FTRIDYN-Xolotl:step')
            print('\n')
            for ftx_comp, ftx in list(self.ftx_components.items()): #.values():
                print(' ')
                print('\t for subworkflow ', ftx_comp, ' with sim_name ', ftx['sim_name']) 
                print('\t Wait for driver:init to be done')
                self.services.wait_call(self.running_components['{}:driver:init'.format(ftx['sim_name'])], True)
                self.services.update_state()
                sys.stdout.flush()
                print('\t Done waiting for driver:init')
                print('\t Free to call driver:step')
                sys.stdout.flush()
                
                #child['LOG_FILE']='log.step.{}'.format(child_comp)
                self.running_components['{}:driver:step'.format(ftx['sim_name'])] = self.services.call_nonblocking(ftx['driver'],
                                                                                                                     'step', 
                                                                                                                     timeStamp,
                                                                                                                     **ftx) #**keys)
                sys.stdout.flush()
                
            for ftx_comp, ftx in list(self.ftx_components.items()):
                print(' ')
                print('\t for subworkflow ', ftx_comp, ' with sim_name ', ftx['sim_name'])
                print('\t Wait for driver:step to be done')
                self.services.wait_call(self.running_components['{}:driver:step'.format(ftx['sim_name'])], True)
                self.services.update_state()
                sys.stdout.flush()
                print('\t Done waiTing for driver:step')
                sys.stdout.flush()

            for ftx_comp, ftx in list(self.ftx_components.items()):
                del self.running_components['{}:driver:init'.format(ftx['sim_name'])]
                del self.running_components['{}:driver:step'.format(ftx['sim_name'])] 
                sys.stdout.flush()


            print('\n')
            print('--------------------')
            print('FTX out --> SOLPS in')
            print('--------------------')
            print('\n')

            #print('\n')
            print('NOT IMPLEMENTED YET')
            print('\n')

            self.services.update_state()
            
            timeStamp+=self.time_step
            t_count+=1
            ## END LOOP OVER TIME HERE
            print('\n \t ------------')
            print('\t END OF LOOP ', t_count)
            print('\t TIME NOW IS t = ',timeStamp)
            print(' \t ------------\n')
            self.services.update_state()

            sys.stdout.flush()
            
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
        for ftx_comp, ftx in list(self.ftx_components.items()):
            self.running_components['{}:driver:finalize'.format(ftx['sim_name'])] = self.services.call_nonblocking(ftx['driver'],
                                                                                                                     'finalize', timeStamp)

            print(' ')
            print('\t Subworkflow ', ftx_comp, ' FINALIZED!')
            print(' ')


#  Wait until all the dependent components are finished.
        self.services.wait_call_list(list(self.running_components.values()), True)
