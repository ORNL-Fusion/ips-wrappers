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
import pickle
import difflib
import numpy as np
from ips_solps_ftx.python_scripts_for_coupling import plasmaOut2ftxIn
from ips_solps_ftx.python_scripts_for_coupling import write_ftxOut
#from ips_solps_ftx.python_scripts_for_coupling import SOLPS_heatflux_for_FTX #maybe generalize name
from ips_solps_ftx.python_scripts_for_coupling import SOLPS_outputs_for_FTX
from ips_solps_ftx.python_scripts_for_coupling import average_SOLPS_input
from ips_solps_ftx.python_scripts_for_coupling import updateSOLPSinput
from ips_solps_ftx.python_scripts_for_coupling import read_impact_angle
from ips_solps_ftx.python_scripts_for_coupling import read_impact_energy

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
        print('Running in cwd:')
        print(cwd)
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
        ftx_locs_string = self.services.get_config_param('SUBWF_LOC')
        #ftx_conf = self.services.get_config_param('SUBWF_COMPONENT_CONF')
        solps_conf = self.services.get_config_param('SOLPS_COMPONENT_CONF')
        ftx_input_dir=self.services.get_config_param('FTX_INPUT_DIR')
        self.solps_input_dir=self.services.get_config_param('SOLPS_INPUT_DIR')

        print('read various parent config file parameters:')
        self.ftx_locs=[0]*num_ftx
        if len(ftx_locs_string.split())<num_ftx:
            print('\t WARNING: less grid points than number of ftx runs')
            print('\t len(ftx_locs_string.split()) = ', len(ftx_locs_string.split()), 'and num_ftx = ', num_ftx )
            print('\t other grid points set to zero ; continue at your own risk')
        elif len(ftx_locs_string.split())>num_ftx:
            print('\t WARNING: more grid points than number of ftx ; ignore the rest of gridpoints!')
            print('\t len(ftx_locs_string.split()) = ', len(ftx_locs_string.split()), 'and num_ftx = ', num_ftx )            
        for i in range(0, num_ftx):        
            self.ftx_locs[i]=int(ftx_locs_string.split()[i])
        print('\t Running ', num_ftx, ' FTX simulation(s), for SOLPS grid point(s): ', self.ftx_locs)
        sys.stdout.flush()
        
        self.solps_output_file=list(self.services.get_config_param('SOLPS_OUTPUT_FORMAT'))

        if (self.print_test):
            print(' ')
            print('will be reading solps output from:', self.solps_output_file)
        
        self.init_time = float(self.services.get_config_param('INIT_TIME'))
        self.end_time = float(self.services.get_config_param('END_TIME'))
        self.time_step = float(self.services.get_config_param('TIME_STEP'))
        self.init_loop_n = int(self.services.get_config_param('INIT_LOOP_N'))
        self.time_decimal = int(self.services.get_config_param('TIME_DECIMAL'))

        try:
            self.H_plasma = str(self.services.get_config_param('H_PLASMA'))
            print('TEST: defined as H plasma = ', self.H_plasma, 'in config file')
            print('TEST: will pass to FTX and cross-check there against H being in plasmaSpecies')
        except Exception as e:
            print(e)
            print('TEST: No H plasma value defined in config file. ')
            print('TEST: Will not assign a value or pass it to FTX (default in FTX is no)')            
            
        ftx_keys = {} 


        #stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        #  Sub workflows require manual setup. First a sub directory must be created.
        #  Then copying of the input files must be performed manually. The first
        #  argument of create sub workflow doesn't appear to do anything.

        #create empty files for those defined in plasma state to avoid errors:
        #to only stage/update what the driver needs, use PLASMA_STATE_FILES ; otherwise use self.STATE_FILES 
        plasma_state_file = self.services.get_config_param('PLASMA_STATE_FILES')  
        plasma_state_list = plasma_state_file.split()
        print(' ')
        print('create empty plasma state files if they donnot exist...')
        for index in range(len(plasma_state_list)):
            open(plasma_state_list[index], 'a').close()
            if self.print_test:
                print('\t created: ', plasma_state_list[index])
        print('... plasma state files created')
        print(' ')
        
        
        #CREATE THE SOLPS SUB-WORKFLOW:
        print('\n')
        print('------------------------------')
        print('CREATE THE SOLPS SUB-WORKFLOW:')
        print('------------------------------')
        print('\n')
        #To-Do: remove this is it's not being used
        self.nested_components['component_SOLPS'] = {'sim_name': None,
                                                     'init': None,
                                                     'driver': None #,
                                                     }
                                                     #'INPUT_DIR' : 'solps_dir',
                                                     #'LOG_FILE' : 'log.component_SOLPS.warning'
                                                     #}
        if (self.print_test):
            print('Defined dictionary : ')
            print('\t', (self.nested_components['component_SOLPS']))

        print('Creating sub workflow with: ')
        print('\t name: component_SOLPS')
        print('\t config: ', solps_conf)
        print('\t keys: SIM_NAME : solps_iter', 'LOG_FILE : log.component_SOLPS.warning')
        print('\t INPUT DIR ', self.solps_input_dir)
        print(' ')

        #make copy of solps_input_dir as init_solps_input_dir
        print('Some SOLPS input files will be updated during the simulation')
        print('Keep a copy of initial SOLPS input dir...')
        shutil.copytree(self.solps_input_dir, self.solps_input_dir+'_initBackup')
        print('... done copying input dir as: ', self.solps_input_dir+'_initBackup')
        print(' ')
        
        (self.nested_components['component_SOLPS']['sim_name'],
         self.nested_components['component_SOLPS']['init'],
         self.nested_components['component_SOLPS']['driver']) = self.services.create_sub_workflow('component_SOLPS',
                                                                                                  solps_conf, 
                                                                                                  {'SIM_NAME' : 'solps_iter',
                                                                                                   'LOG_FILE' : 'log.component_SOLPS.warning'},
                                                                                                  self.solps_input_dir)

        print('creating SOLPS sub-workflow DONE!')
        if (self.print_test):
            print('\t SOLPS component is:') 
            print('\t', self.nested_components['component_SOLPS'])
        print('\n')

        ## CREATE SUB-WORKFLOW FOR MULTIPLE (num_ftx) FTRIDYN-XOLOTL RUNS, IN PARALLEL ##
        print('----------------------------------------------------')
        print('CREATE SUB-WORKFLOW FOR MULTIPLE FTRIDYN-XOLOTL RUNS ')
        print('----------------------------------------------------')
        print('\n')
        for i in range(0, num_ftx):
            ftx_comp = 'ftx_{}'.format(i)

            print('-------- ', ftx_comp, ' --------\n')
            
            ftx_keys['LOG_FILE'] = 'log.{}'.format(ftx_comp)
            ftx_keys['SIM_NAME'] = ftx_comp
            ftx_keys['INPUT_DIR'] = 'input_{}'.format(ftx_comp)
            ftx_keys['INPUT_FILES'] = self.services.get_config_param('INPUT_FTX')
            ftx_keys['time_decimal'] = int(self.time_decimal)
            try:
                ftx_keys['H_PLASMA'] = str(self.H_plasma)
                print('H_PLASMA=',self.H_plasma,'passed to ftx_keys')
            except Exception as e:
                print(e)
                print('H_PLASMA will not be defined in ftx_keys.')
            print(' ')
                
            self.ftx_components[ftx_comp] = {
                                                'sim_name'  : None, 
                                                'init'      : None,
                                                'driver'    : None,
                                                'INPUT_DIR' : ftx_keys['INPUT_DIR'],
                                                'LOG_FILE'  : ftx_keys['LOG_FILE']
                                                }
            if (self.print_test):
                print('Defined dictionary : ')
                print('\t', (self.ftx_components[ftx_comp]))

            #  Input files will be staged from this directory.
            print('In cwd:')
            print('...create input directory: ', self.ftx_components[ftx_comp]['INPUT_DIR'])

            if os.path.exists(self.ftx_components[ftx_comp]['INPUT_DIR']):
                shutil.rmtree(self.ftx_components[ftx_comp]['INPUT_DIR'])
            os.mkdir(self.ftx_components[ftx_comp]['INPUT_DIR'])

            #  Copy files to the created directory.
            input_file_list=ftx_keys['INPUT_FILES'].split()
            print('... and copy input files:') # , input_file_list)
            print('    from ftx input directory :', ftx_input_dir)
            print('    to sub-wf input file directory : ' , self.ftx_components[ftx_comp]['INPUT_DIR'])
            for input_file in input_file_list:
                print('\t', input_file)
                #print('\t from submission directory to ', input_file , ' from ', self.SUBMIT_DIR, ' to ', self.ftx_components[ftx_comp]['INPUT_DIR']) #<-- without long paths 
                shutil.copyfile(ftx_input_dir+'/'+input_file,self.ftx_components[ftx_comp]['INPUT_DIR']+'/'+input_file)
            print('...done copying input files')

            ftx_conf = self.services.get_config_param('SUBWF_COMPONENT_CONF')+'_{}'.format(i)
            if self.print_test:
                print('config file fot this subworkflow is: ', ftx_conf)
            
            print(' ')
            print('Create_sub_workflow with parameters:')
            print('\t ftx_comp = ', ftx_comp)
            print('\t ftx_conf = ', ftx_input_dir+'/'+ftx_conf)
            print('\t keys = ', ftx_keys)
            print('\t input_dir = ', self.ftx_components[ftx_comp]['INPUT_DIR'])
            print('\n')
            sys.stdout.flush()
            (self.ftx_components[ftx_comp]['sim_name'],
             self.ftx_components[ftx_comp]['init'],
             self.ftx_components[ftx_comp]['driver']) = self.services.create_sub_workflow(ftx_comp, 
                                                                                              ftx_input_dir+'/'+ftx_conf, 
                                                                                              ftx_keys, 
                                                                                              self.ftx_components[ftx_comp]['INPUT_DIR'])
            print('creating FTX sub-workflow DONE!')
            if (self.print_test):
                print('FTX component is:') 
                print('\t sim_name : ',self.ftx_components[ftx_comp]['sim_name'],'init : ',self.ftx_components[ftx_comp]['init'],'driver : ',self.ftx_components[ftx_comp]['driver'], ', INPUT_DIR : ', self.ftx_components[ftx_comp]['INPUT_DIR']) 
                print('\t keys : ', ftx_keys)
            print('\n')

        print('-------------------------')
        print('DONE SETTING UP WORKFLOWS')
        print('-------------------------')
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
        print('Running in cwd = ', cwd, ':')

        ftx_input_dir=self.services.get_config_param('FTX_INPUT_DIR')
        #if we ever define specific outFile / logFile, do so here:
        
        #Insert time-loop here:
        time = round(self.init_time, self.time_decimal)         
        t_count = self.init_loop_n


        #solps init component before while loop:
        print('\n')
        print('\t SOLPS:init')
        print('\n')

        print('This is run once, in the beginning of the simulation')

        print('run SOLPS init:init')
        self.async_queue['component_SOLPS:init:init'] = self.services.call(self.nested_components['component_SOLPS']['init'], 'init', time)
        print('init:init done \n')
        self.services.stage_state()
        
        print('run SOLPS init:step')
        self.async_queue['component_SOLPS:init:step'] = self.services.call(self.nested_components['component_SOLPS']['init'], 'step', time)
        print('init:step done \n')
        self.services.stage_state()
        
        print('run SOLPS init:finalize')
        self.async_queue['component_SOLPS:init:finalize'] = self.services.call(self.nested_components['component_SOLPS']['init'], 'finalize', time)
        print('init:step done \n')
        self.services.stage_state()
        
        while time < self.end_time:

            print('\n Starting iteration for time : ', time)
            
            print('\n')
            print('-------------------')
            print('-------------------')
            print('run SOLPS driver')
            print('-------------------')
            print('-------------------')
            print('\n')
            
            print('run SOLPS driver:init')
            self.async_queue['component_SOLPS:driver:init'] = self.services.call(self.nested_components['component_SOLPS']['driver'], 'init', time)
            print('driver:init done \n')
            self.services.stage_state()

            print('run SOLPS driver:step')
            self.async_queue['component_SOLPS:driver:step'] = self.services.call(self.nested_components['component_SOLPS']['driver'], 'step', time) 
            print('driver:step done \n')
            self.services.stage_state()

            print('run SOLPS driver:finalize')
            self.async_queue['component_SOLPS:driver:finalize'] = self.services.call(self.nested_components['component_SOLPS']['driver'], 'finalize', time)
            print('driver:finalize done \n')
            self.services.stage_state()

            print('stage subworkflow outputs of component_SOLPS... ')
            sys.stdout.flush()
            try:
                self.services.stage_subflow_output_files(subflow_name='component_SOLPS')
                print('... succesfully staged component_SOLPS outputs')
            except:
                print('... could not stage component_SOLPS outputs')
            print('... done staging subworflow outputs after driver:step\n')
            sys.stdout.flush()            

            #save solps log file so that it's not overwritten in every iteration:
            print(' ')
            print('Copy SOLPS log file to avoid overwriting in next iteration')
            print('solps.log as solps.log_t'+str(time))
            shutil.copyfile('solps.log', 'solps.log_t'+str(time))

            print('\n')
            print('COMPLETED SOLPS')
            
            
            ## pass plasma information by:
            ## 1 - read output of SOLPS from solpsOut.txt (one per child) using SOLPS_outputs_for_FTX script 
            ## 2 - do ops to calculate/format inputs for FTX in plasmaOut2ftxIn:
            ##     a) check inputs expected from SOLPS / that need special care
            ##     i.  inputs that are single values:
            ##         particle flux, heat flux
            ##     ii. inputs that are given for each species:
            ##         for now: always include information for all 4 species: simpler to handle
            ##                  and doesn't add comp time as the ftx workflow skips anything with fluxFraction=0
            ##         fluxFraction, [Ti, Te, Z] --> E_in ; B-field --> A_in
            ##     b) print all others 'as is'
            ## 3 - print ftxIn.txt
            ## 4 - pass time of loop to FTX in its own file
            ## 5 - copy ftxIn and time-file to child's input directory

            print('\n')
            print('--------------------')
            print('--------------------')
            print('SOLPS out --> FTX in')
            print('--------------------')
            print('--------------------')
            print('\n')


            #here implement parsing SOLPS output

            #check that the files needed for analysis were staged correctly:
            #for now hard coded ; define files in config file if we want
            if (self.print_test):
                print('Check that staged files exist and are not empty')
                print(' ')
                for f in ['b2fgmtry', 'b2fstate', 'fort.44']: #'b2fpardf', 'b2frates',  also exist but not needed here
                    if os.path.exists(f):
                        open_file = open(f)
                        lines_file = open_file.readlines()
                        if len(lines_file) > 0:
                            print('\t Succesfully opened file (not empty) :', open_file.name)
                        #print('\t len(',open_file.name,') = ', len(lines_file))
                        open_file.close()
                    else:
                        print('\t WARNING: no ', f ,' found!')
                        print("\t SOLPS out --> FTX in won't work")
                print(' ')
                sys.stdout.flush()

                    
            # 2 - Call SOLPS_heatflux_for_FTX to write SOLPS output into FTX readable form

            #these parameters don't change with the FTX location
            self.SOLPSoutputs={}
            self.SOLPSoutputs['b2fstate_file'] = 'b2fstate'
            self.SOLPSoutputs['b2fgmtry_file'] = 'b2fgmtry'
            self.SOLPSoutputs['fort44_file'] = 'fort.44'
            SOLPSoutputs_log='log.SOLPSoutputs_t'+str(time)

            
            for ftx_comp, ftx in list(self.ftx_components.items()):
                i=int(''.join(list(filter(str.isdigit, ftx_comp))))
                solps_outFile=self.solps_output_file[0]+'_'+str(i)+'.'+self.solps_output_file[1] #0=filename (solpsOut) ; 1=fomat (txt)
                #solps_outFile='solpsOutput_'+str(i)+'.test' #test file to avoid overwritting the one from input

                #these parameters are location-dependent
                self.SOLPSoutputs['rad_grid'] = self.ftx_locs[i]
                self.SOLPSoutputs['pickleOutputFile']=solps_outFile
                self.SOLPSoutputs['print_test']=self.print_test
                self.SOLPSoutputs['logFile']=SOLPSoutputs_log

                #launch SOLPSoutputs instead of calling function:
                pkl_SOLPSoutputs_file='SOLPSoutputs.pkl' #rm cwd+'/'
                pickle.dump(self.SOLPSoutputs, open(pkl_SOLPSoutputs_file, "wb" ) )
                shutil.copyfile(pkl_SOLPSoutputs_file, 'SOLPSoutputs_'+str(i)+'.pkl')

                print('Launch SOLPS_outputs_for_FTX ... ')
                if (self.print_test):
                    print('\t with inputs: ')
                    print('\t rad_grid : ', self.SOLPSoutputs['rad_grid'])
                    print('\t b2fstate_file :', self.SOLPSoutputs['b2fstate_file'])
                    print('\t b2fgmtry_file :', self.SOLPSoutputs['b2fgmtry_file'])
                    print('\t fort44_file :', self.SOLPSoutputs['fort44_file'])
                    print('\t SOLPSoutput file : ', self.SOLPSoutputs['pickleOutputFile'])
                    print('\t print_test : ', self.SOLPSoutputs['print_test'])
                    print('\t logFile : ', self.SOLPSoutputs['logFile'])
                    print(' ')
                sys.stdout.flush()

                SOLPSoutputs_script = 'SOLPS_outputs_for_FTX.py'
                print('\t Launch python script: ')
                print('\t', SOLPSoutputs_script)                
                sys.stdout.flush()
                
                if self.print_test:
                    print('\t launch task : ', SOLPSoutputs_script)
                    print('\t with pkl file: ', 'SOLPSoutputs.pkl')
                    print('\t logFile : ', SOLPSoutputs_log)
                    print(' ')
                task_id_SOLPSoutputs = self.services.launch_task(1,self.services.get_working_dir(),
                                                                 SOLPSoutputs_script, logFile=SOLPSoutputs_log)
                ret_val_SOLPSoutputs = self.services.wait_task(task_id_SOLPSoutputs)

                if self.print_test:
                    print('\t after running ', SOLPSoutputs_script)
                    print('\t ret_val_SOLPSoutputs = ', ret_val_SOLPSoutputs)
                    print(' ')
                sys.stdout.flush()

                print('... done with SOLPSoutputs for ftx '+str(i)+' \n')
                sys.stdout.flush()

                print('READ AND FORMAT INPUTS TO FTX WORKFLOW(S)')
                print(' ')
            
            #for ftx_comp, ftx in list(self.ftx_components.items()): ## I think we can just continue loop, no need to break and start again
            #    print(' ')
            #    i=int(''.join(list(filter(str.isdigit, ftx_comp))))
            #    if (self.print_test):
            #        print('index of ', ftx_comp, ' is ', i)
            #        sys.stdout.flush()

                # 3 - use inputs from hPIC:
                #     implemented for both angle scalar & distribution; energy scalar & distribution
                #     different impact angle/energy models: solps, readFileScalar, readFileDistrib
                #     implement different functions within read_impact_angle : scalar & distrib (idem for energy)
                #       for a single angle value, read_impact_angle.scalar returns value of angle
                #       for distributions [not implemented yet!] read_impact_angle.distrib reads file name of distrib, instead of the value

                #3.a angle from hPIC
                impactAngleModel=self.services.get_config_param('IMPACT_ANGLE_MODEL')
                if self.print_test:
                    print('impact angle model is ', impactAngleModel)
                if impactAngleModel=='solps':
                    print('will use angle provided in input from SOLPS (= B-field)')
                elif impactAngleModel=='readFileScalar':
                    #read a single value for angle
                    # ADD (i) to these two file names if we need a different distribution for each FTX
                    hpicAngleFile=self.services.get_config_param('HPIC_ANGLE_FILE_NAME') #e.g. 'hpicAngles.txt'
                    readAngle_log='log.readImpactAngle_t'+str(time) 
                    print('read single angle value from file: ', hpicAngleFile)
                    print('\t using script read_impact_angle.scalar with inputs:')
                    print('\t t = ', time, ' ; print_test =  ; ', self.print_test, ' ; logFile =', readAngle_log)
                    [ftxImpactAngle_D,ftxImpactAngle_C] = read_impact_angle.scalar(time=time, inputAngleFile=hpicAngleFile, print_test=self.print_test, logFile=readAngle_log)
                    sys.stdout.flush()
                    
                    print('\t for time ', time, ' impact angles of D and C are ', ftxImpactAngle_D, ' and ', ftxImpactAngle_C,' deg')
                    #replace angle in solpsOut pkl file
                    orig_solps_outFile=self.solps_output_file[0]+'_'+str(i)+'_orig.'+self.solps_output_file[1]+'_t'+str(time)
                    temp_solps_outFile=self.solps_output_file[0]+'_'+str(i)+'_temp.'+self.solps_output_file[1]
                    shutil.copy(solps_outFile,orig_solps_outFile) #save copy of original
                    shutil.move(solps_outFile,temp_solps_outFile) #mv solps_outFile to temp
                    solpsOut_dic=pickle.load( open(temp_solps_outFile, "rb" ))
                    #if self.print_test:
                    #    print('\t before changing the angle, solpsOut dictionary is:')
                    #    print('\t', solpsOut_dic)
                    solpsOut_dic['inputAngle']=[ftxImpactAngle_D, ftxImpactAngle_C]
                    pickle.dump(solpsOut_dic, open(solps_outFile, "wb" ) ) #same file name, updated values
                    print('changed solps dictionarys inputAngle to ', [ftxImpactAngle_D, ftxImpactAngle_C])
                    #if self.print_test:
                    #    print('\t after changing the angle, solpsOut dictionary is:')
                    #    print('\t', solpsOut_dic)
                    print(' ')
                    sys.stdout.flush()
                    
                elif impactAngleModel=='readFileDistrib':
                    print('angle distribution model not implemented yet')
                    print('do nothing')
                else:
                    print('angle model not recognized')
                    print('do nothing')
                sys.stdout.flush()

                #3.b energy from hPIC
                impactEnergyModel=self.services.get_config_param('IMPACT_ENERGY_MODEL')
                if self.print_test:
                    print('impact energy model is ', impactEnergyModel)
                if impactEnergyModel=='solps':
                    print('will use impact energy provided in input from SOLPS (based on Ti and Te)')
                elif impactEnergyModel=='readFileScalar':
                    #read a single value for energy (per species)
                    # ADD (i) to these two file names if we need a different distribution for each FTX
                    hpicEnergyFile=self.services.get_config_param('HPIC_ENERGY_FILE_NAME') #e.g. 'hpicEnergies.txt'
                    readEnergy_log='log.readImpactEnergy_t'+str(time)
                    print('read single energy value from file: ', hpicEnergyFile)
                    print('\t using script read_impact_energy.scalar with inputs:')
                    print('\t t = ', time, ' ; print_test =  ; ', self.print_test, ' ; logFile =', readEnergy_log)
                    [ftxImpactEnergy_D,ftxImpactEnergy_C] = read_impact_energy.scalar(time=time, inputEnergyFile=hpicEnergyFile, print_test=self.print_test, logFile=readEnergy_log)
                    sys.stdout.flush()                    

                    print('\t for time ', time, ' impact energies of D and C are : ', ftxImpactEnergy_D, ' and ', ftxImpactEnergy_C, ' eV' )
                    #replace energy in solpsOut pkl file
                    orig_solps_outFile=self.solps_output_file[0]+'_'+str(i)+'_orig.'+self.solps_output_file[1]+'_t'+str(time)
                    temp_solps_outFile=self.solps_output_file[0]+'_'+str(i)+'_temp.'+self.solps_output_file[1]
                    #make copy of original (if file doesn't exist ; otherwise it was already copied in replacing angle)
                    if os.path.exists(orig_solps_outFile):
                        if self.print_test:
                            print('\t ', orig_solps_outFile, 'already exists. no need to copy')
                    else:
                        if self.print_test:
                            print('\t ', orig_solps_outFile, 'does not exists. make copy')
                        shutil.copy(solps_outFile,orig_solps_outFile)
                    shutil.move(solps_outFile,temp_solps_outFile)
                    solpsOut_dic=pickle.load( open(temp_solps_outFile, "rb" ))
                    #if self.print_test:
                    #    print('\t before changing the energy, solpsOut dictionary is:')
                    #    print('\t', solpsOut_dic)
                    solpsOut_dic['inputEnergy']=[ftxImpactEnergy_D,ftxImpactEnergy_C]
                    pickle.dump(solpsOut_dic, open(solps_outFile, "wb" ) ) #same file name, updated values
                    print('changed solps dictionarys inputEnergy to ', [ftxImpactEnergy_D,ftxImpactEnergy_C])
                    #if self.print_test:
                    #    print('\t after changing the energy, solpsOut dictionary is:')
                    #    print('\t', solpsOut_dic)
                    print(' ')
                    sys.stdout.flush()

                elif impactEnergyModel=='readFileDistrib':
                    print('energy distribution model not implemented yet')
                    print('do nothing')

                else:
                    print('energy model not recognized')
                    print('do nothing')
                sys.stdout.flush()

                #END of 3 - use inputs from hPIC

                #save solps_outFile, regardless of modified by hPIC or not
                if self.print_test:
                    print('save SOLPS output file ', solps_outFile, 'as ', solps_outFile+'_t'+str(time))
                shutil.copyfile(solps_outFile, solps_outFile+'_t'+str(time))

                
                # 4 - REFORMAT SOLPS OUTPUT FILE                
                
                ## read output of SOLPS from solpsOut.txt (one per child) - DONE?
                print('re-format SOLPS output into FTX input: ')
                #these three are used later as well, so define here
                ftxInFileFormat=list(self.services.get_config_param('FTX_INPUT_FORMAT'))
                ftxInFileName=ftxInFileFormat[0]+'_'+str(i)+'.'+ftxInFileFormat[1]
                #solps_outFile=self.solps_output_file[0]+'_'+str(i)+'.'+self.solps_output_file[1] #not needed anymore, bc we're in the same loop 0=filename (solpsOut) ; 1=fomat (pkl)
                p2ftx_log='log.plasmaOut2ftxIn'+'_t'+str(time) #rm cwd+'/'

                self.plasmaOut2ftxIn={}
                self.plasmaOut2ftxIn['plasmaOutFile'] = solps_outFile #ftx_input_dir+'/'+solps_outFile #for now, located in ftx input dir ; when it's solps output, it should be in driver dir
                self.plasmaOut2ftxIn['ftxInFile'] = ftxInFileName #rm cwd+'/'
                self.plasmaOut2ftxIn['logFile']= p2ftx_log
                self.plasmaOut2ftxIn['print_test'] = self.print_test

                #launch plasmaOut2ftxIn instead of calling function:            
                pkl_p2ftx_file='plasmaOut2ftxIn.pkl' #rm cwd+'/'
                pickle.dump(self.plasmaOut2ftxIn, open(pkl_p2ftx_file, "wb" ) )

                print('Launch plasmaOut2ftxIn... ')
                if (self.print_test):
                    print('\t with inputs: ')
                    print('\t input file: ', solps_outFile )#<-- without long paths: self.plasmaOut2ftxIn['plasmaOutFile'])
                    print('\t output file: ', ftxInFileName) #<-- without cwd: self.plasmaOut2ftxIn['ftxInFile'])
                    print('\t print_test: ', self.plasmaOut2ftxIn['print_test'])
                    print('\t output log: ', 'log.plasmaOut2ftxIn'+'_t'+str(time)) #<-- without cwd: self.plasmaOut2ftxIn['logFile'])
                    print(' ')
                sys.stdout.flush()

                plasma2ftx_script = 'plasmaOut2ftxIn.py'
                print('\t Launch python script: ')
                print('\t',plasma2ftx_script)
                    
                sys.stdout.flush()
                if self.print_test:
                    print('\t launch task : ', plasma2ftx_script)
                    print('\t with pkl file: ', 'plasmaOut2ftxIn.pkl') # <-- without cwd: pkl_p2ftx_file)
                    print('\t logFile : ', 'log.plasmaOut2ftxIn'+'_t'+str(time)) #<-- without cwd: p2ftx_log)
                    #<-- withiout long paths: print('\t self.services.get_working_dir() = ', self.services.get_working_dir())
                
                task_id_p2ftx = self.services.launch_task(1,self.services.get_working_dir(),
                                                          plasma2ftx_script, logFile=p2ftx_log)
                ret_val_p2ftx = self.services.wait_task(task_id_p2ftx)
                
                if self.print_test:
                    print('\t after running ', plasma2ftx_script)
                    print('\t ret_val_p2ftx = ', ret_val_p2ftx)
                    print(' ')
                sys.stdout.flush()
                #self.plasmaOut2ftxIn.clear()
                print('... done with plasmaOut2ftxIn \n')                
                ## END OF 4 - solpsOut --> ftxIn

                #5- copy ftxIn to ftx input directory
                print('copy input file to FTX input directory')
                print('\t from: ', ftxInFileName) ##<-- without cwd: cwd+ftxInFileName
                print('\t to: ', self.ftx_components[ftx_comp]['INPUT_DIR']+'/ftxInput.txt')
                shutil.copyfile(ftxInFileName,self.ftx_components[ftx_comp]['INPUT_DIR']+'/ftxInput.txt')
                print('copy SOLPS output file to FTX input directory')
                print('\t from: ', solps_outFile)
                print('\t to: ', self.ftx_components[ftx_comp]['INPUT_DIR']+'/solpsOut.pkl')
                shutil.copyfile(solps_outFile,self.ftx_components[ftx_comp]['INPUT_DIR']+'/solpsOut.pkl') #modify ftx_input_dir to be driver dir when output comes from solps
                
                print('\n')
                sys.stdout.flush()

                #save ftxIn and solpsOut for each loop (with time-stamp)
                shutil.copyfile(solps_outFile, solps_outFile+'_t'+str(time)) #modify ftx_input_dir to be driver dir when output comes from solps
                shutil.copyfile(ftxInFileName,ftxInFileName+'_t'+str(time)) 
                print('to save input files for each time-loop, copied:')
                print('\t ', solps_outFile, ' as ', solps_outFile+'_t'+str(time))
                print('\t ', ftxInFileName, ' as ', ftxInFileName+'_t'+str(time))
                sys.stdout.flush()
                #end of 5 - save ftxIn
                
                #update plasma state file:
                self.services.update_state()

                    
            # 6 - add time to its own file
            #  For now keep this section in driver
            #  Better move to its own function/script
            print('write time input file:')
            timeFileName=self.services.get_config_param('TIME_FILE_NAME')
            print('\t file name: ', timeFileName)
            timeFile=open(timeFileName, "w")
            print('\t INIT_TIME : ', time)
            print('\t END_TIME : ', time + self.time_step)
            print('\t LOOP_TIME_STEP : ', self.time_step)
            print('\t LOOP_N : ', t_count)
            timeFile.write('INIT_TIME={0}\n'.format(time))
            timeFile.write('END_TIME={0}\n'.format(time + self.time_step))
            timeFile.write('LOOP_TIME_STEP={0}\n'.format(self.time_step))
            timeFile.write('LOOP_N={0}\n'.format(t_count))            
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

            #copy time parameter to the input folder of each FTX subworkflow
            print('copying time parameter: ')
            for ftx_comp, ftx in list(self.ftx_components.items()):
                print(' ')
                #to get i, these two ways should work:
                i=int(''.join(list(filter(str.isdigit, ftx_comp))))
                print('\t from: ', timeFileName)
                print('\t to: ', self.ftx_components[ftx_comp]['INPUT_DIR']+'/'+timeFileName)
                shutil.copyfile(timeFileName,self.ftx_components[ftx_comp]['INPUT_DIR']+'/'+timeFileName)
            shutil.copyfile(timeFileName,timeFileName+'_t'+str(time))
            print('\t ', timeFileName , ' as ', timeFileName+'_t'+str(time))
            sys.stdout.flush()
            
            
            print('\n')
            print('-------------------')
            print('-------------------')
            print('run FTRIDYN-Xolotl')
            print('-------------------')
            print('-------------------')
            
            #RUN init AND step OF CHILD (FT-X) SUB-WORKFLOW            
            print(' ')
            print('FTRIDYN-Xolotl:init')
            print(' ')
            
            #  Loop over the children and all the initize component.
            for ftx_comp, ftx in list(self.ftx_components.items()):
                print(' ')
                print('\t for subworkflow ', ftx_comp)
                print('\t Call driver:init with dictionary ', ftx)
                print(' ')
                self.running_components['{}:driver:init'.format(ftx['sim_name'])] = self.services.call_nonblocking(ftx['driver'],
                                                                                                                     'init',
                                                                                                                     time,
                                                                                                                     **ftx) #**keys
                print(' ')
                sys.stdout.flush()
            self.services.stage_state()
            
            #  Loop over the children and all the initialize driver component.
            print(' ')
            print('FTRIDYN-Xolotl:step')
            print(' ')
            for ftx_comp, ftx in list(self.ftx_components.items()): #.values():
                print('\t for subworkflow ', ftx_comp, ' with sim_name ', ftx['sim_name']) 
                print('\t Wait for driver:init to be done')
                self.services.wait_call(self.running_components['{}:driver:init'.format(ftx['sim_name'])], True)
                self.services.update_state()
                sys.stdout.flush()
                print('\t Done waiting for driver:init')
                print('\t Free to call driver:step')
                print(' ')
                sys.stdout.flush()

                self.running_components['{}:driver:step'.format(ftx['sim_name'])] = self.services.call_nonblocking(ftx['driver'],
                                                                                                                     'step', 
                                                                                                                     time,
                                                                                                                     **ftx) #**keys)
                sys.stdout.flush()
                
            ftxOutFileFormat=list(self.services.get_config_param('FTX_OUTPUT_FORMAT'))
            ftxOutFileName=ftxOutFileFormat[0]+'.'+ftxOutFileFormat[1] # without i, so all ftx output goes into the same file
                
            for ftx_comp, ftx in list(self.ftx_components.items()):
                print('\t for subworkflow ', ftx_comp, ' with sim_name ', ftx['sim_name'])
                print('\t Wait for driver:step to be done')
                self.services.wait_call(self.running_components['{}:driver:step'.format(ftx['sim_name'])], True)
                self.services.update_state()
                sys.stdout.flush()
                print('\t Done waiting for driver:step of ', ftx_comp, '\n')
                print(' ')
                sys.stdout.flush()

                print('--------------------')
                print('FTX out --> SOLPS in')
                print('  for ',ftx_comp)
                print('--------------------\n')

                print('stage ', ftx_comp , ' subworkflow outputs after driver:step') 
                self.services.stage_subflow_output_files(subflow_name=ftx_comp) #tried arguments: ftx_comp, ftx['sim_name'] 
                print('... done staging subworflow outputs after driver:step\n')

                #many file names hard coded for now ; could consider linking to config params
                i=int(''.join(list(filter(str.isdigit, ftx_comp))))

                if (self.print_test):
                    print('Check that staged files exist and are not empty')
                    print(' ')
                    for f in ['last_TRIDYN.dat', 'tridyn.dat', 'retentionOut.txt', 'log.ftx_{}'.format(i)]:
                        if os.path.exists(f):
                            open_file = open(f)                            
                            lines_file = open_file.readlines()
                            if len(lines_file) > 0:
                                print('\t Succesfully opened file (not empty) :', open_file.name)                                
                            #print('\t len(',open_file.name,') = ', len(lines_file))
                            open_file.close()
                        else:
                            print('\t WARNING: no ', f ,' found!')
                            print("\t SOLPS out --> FTX in won't work")
                    print(' ')
                    sys.stdout.flush()

                write_ftxOut_log='log.write_ftxOut_{}'.format(i)+'_t'+str(time) #rm cwd+'/'
                self.write_ftxOut={}
                self.write_ftxOut['grid'] = self.ftx_locs[i]
                self.write_ftxOut['last_tridyn'] = 'last_TRIDYN_{}.dat'.format(i)
                self.write_ftxOut['log_ftx']= 'log.ftx_{}'.format(i)
                self.write_ftxOut['tridyn']='tridyn_{}.dat'.format(i) 
                self.write_ftxOut['retentionFile']='retentionOut_{}.txt'.format(i)
                self.write_ftxOut['print_test'] = self.print_test
                self.write_ftxOut['logFile'] = write_ftxOut_log
                self.write_ftxOut['outFile'] = ftxOutFileName #ftxOutFileStd #rm cwd+'/'

                #to avoid overwritting, copy files with the ftx-index
                if (self.print_test):
                    print('add ftx-index to files staged from subworkflow to avoid overwritting when multiple FTX run in parallel. copy:')
                    print('\t last_TRIDYN.dat as ', self.write_ftxOut['last_tridyn'])
                    print('\t tridyn.dat as ', self.write_ftxOut['tridyn'])
                    print('\t retentionOut.txt as', self.write_ftxOut['retentionFile'])
                    print(' ')
                shutil.copyfile('last_TRIDYN.dat',self.write_ftxOut['last_tridyn'])
                shutil.copyfile('tridyn.dat',self.write_ftxOut['tridyn'])
                shutil.copyfile('retentionOut.txt',self.write_ftxOut['retentionFile'])
                
                #launch write_ftxOut instead of calling function:
                pkl_ftxOut_file='write_ftxOut.pkl' #rm cwd+'/'
                pickle.dump(self.write_ftxOut, open(pkl_ftxOut_file, "wb" ) )

                print('Launch write_ftxOut...')
                if (self.print_test):
                    print('\t with inputs:')
                    print('\t grid/location: ', self.write_ftxOut['grid'] )
                    print('\t last_tridyn: ',self.write_ftxOut['last_tridyn'])
                    print('\t log.ftx_{}'': '.format(i) , self.write_ftxOut['log_ftx'])
                    print('\t tridyn.dat: ',self.write_ftxOut['tridyn'])
                    print('\t retentionFile: ', self.write_ftxOut['retentionFile'])
                    print('\t print_test: ', self.write_ftxOut['print_test'])
                    print('\t output logFile: ', 'log.write_ftxOut'+'_t'+str(time)) #<-- without cwd:
                    print('\t output file: ', ftxOutFileName) #ftxOutFileStd) #<-- without cwd:
                    print(' ')
                sys.stdout.flush()

                write_ftxOut_script = 'write_ftxOut.py'
                print('\t Launch python script: ')
                print('\t',write_ftxOut_script)

                if self.print_test:
                    print('\t launch task : ', write_ftxOut_script)
                    print('\t with pkl file: ', 'write_ftxOut.pkl') # <-- without cwd: 
                    print('\t logFile : ', 'log.write_ftxOut'+'_t'+str(time)) #<-- without cwd                    
                sys.stdout.flush()                    
                    
                task_id_ftxOut = self.services.launch_task(1,self.services.get_working_dir(),
                                                        write_ftxOut_script, logFile=write_ftxOut_log)
                ret_val_ftxOut = self.services.wait_task(task_id_ftxOut)
                if self.print_test:
                    print('\t after running ', write_ftxOut_script)
                    print('\t ret_val_ftxOut = ', ret_val_ftxOut)
                    print(' ')
                sys.stdout.flush()
                print('... done with write_ftxOut \n')
                

            # SOLPS can only take average the FTX recycling factors
            [ave_grid, ave_twall, ave_Rft, ave_Rxol, ave_Rtot] = average_SOLPS_input.average_SOLPS_input(ftxOut_file=self.write_ftxOut['outFile'], print_test=self.print_test, logFile=None)

            #save ftxOut for each loop (with time-stamp) ; move file so that it's not appended next iteration
            shutil.move(self.write_ftxOut['outFile'], ftxOutFileName+'_t'+str(time))
            print('...and save', ftxOutFileName ,' for each time-loop:')
            print('\t move ', ftxOutFileName, ' as ', ftxOutFileName+'_t'+str(time))
            print(' ')

            
            #do whatever needed by SOLPS with these values
            #for now, just print them
            if self.print_test:
                print('from driver, average FTX outputs are:')
                print('\t average grid = ', ave_grid)
                print('\t average T wall = ', ave_twall)
                print('\t average R_FT = ', ave_Rft)
                print('\t average R_Xol = ', ave_Rxol)
                print('\t average R_tot = ', ave_Rtot)
                print(' ')

            #UPDATE SOLPS input.dat
            inputDat='input.dat'
            orig_inputDat=inputDat+'_orig_t'+str(time)
            new_inputDat=inputDat+'_updated_t'+str(time)
            shutil.copyfile(inputDat, orig_inputDat) #save a copy

            #for now, all FTX in DIMES; eventually implement check for grid point falling in rad_grid_name = Left, DIMES or right
            ave_Rft_array = np.array([ave_Rft]) #as array -->, dtype = np.float64)
            #print('TEST TEST: ave_Rft = ', ave_Rft, ' ave_Rft_array =', ave_Rft_array, ' and ave_Rft_array[0] =', ave_Rft_array[0])
            RECYCF = updateSOLPSinput.calc_RECYCF(inputDat, 'fort.44', ['DIMES'], ave_Rft_array) #ave_Rft) 
            ave_Rtot_array = np.array([ave_Rtot])
            ave_twall_array = np.array([ave_twall])
            #print('TEST TEST: ave_Rtot = ', ave_Rtot, ' ave_Rtot_array =', ave_Rtot_array, ' and ave_Rtot_array[0] =', ave_Rtot_array[0])
            #print('TEST TEST: ave_twall = ', ave_twall, ' ave_twall_array =', ave_twall_array, ' and ave_twall_array[0] =', ave_twall_array[0])
            updateSOLPSinput.input_dat_update(orig_inputDat, inputDat, RECYCF, ave_Rtot_array, ave_twall_array, ['DIMES']) 
            shutil.copyfile(inputDat, new_inputDat)

            print('Update values in ', inputDat)
            if self.print_test:
                print('\t copy ', inputDat, 'as ', inputDat+'_updated_t'+str(time))
            
            # JUST FOR TESTING!!
            print('\t compare original ', orig_inputDat, ' and new ', inputDat, ' files')
            with open(orig_inputDat) as file_1:
                file_1_text = file_1.readlines()
 
            with open(new_inputDat) as file_2:
                file_2_text = file_2.readlines()
 
            # Find and print the diff:
            for line in difflib.unified_diff(
                    file_1_text, file_2_text, fromfile='file1.txt',
                    tofile='file2.txt', lineterm=''):
                print('\t found difference between input.dat files: ')
                print('\t', line)
            file_1.close
            file_2.close
            print(' ')
            
            for ftx_comp, ftx in list(self.ftx_components.items()):
                del self.running_components['{}:driver:init'.format(ftx['sim_name'])]
                del self.running_components['{}:driver:step'.format(ftx['sim_name'])] 
                sys.stdout.flush()


            self.services.update_state()


            #copy updated files to input.dat, b2fstate... to self.solps_input_dir
            print('Update SOLPS input file before running next iteration...')
            solps_update_string=self.services.get_config_param('UPDATE_SOLPS_INPUTS')
            update_list=solps_update_string.split()
            for f in update_list:
                print('\t update ', f, ' in solps input dir')
                shutil.copyfile(f, self.solps_input_dir+'/'+f)
            print('... done updating inputs to SOLPS')
            print(' ')
            
            new_time=time+self.time_step
            time = round(new_time, self.time_decimal)
            t_count+=1
            ## END LOOP OVER TIME HERE
            print('\n \t ------------')
            print('\t END OF LOOP ', t_count)
            print('\t TIME NOW IS t = ',time)
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

        time=round(timeStamp, self.time_decimal)
        print('finalize FTRIDYN-Xolotl')
        for ftx_comp, ftx in list(self.ftx_components.items()):
            self.running_components['{}:driver:finalize'.format(ftx['sim_name'])] = self.services.call_nonblocking(ftx['driver'],
                                                                                                                     'finalize', time)

            print(' ')
            print('\t Subworkflow ', ftx_comp, ' FINALIZED!')
            print(' ')


#  Wait until all the dependent components are finished.
        self.services.wait_call_list(list(self.running_components.values()), True)
