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
import param_handler

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
                solps_outFile = self.solps_output_file[0]+'_'+str(i)+'.'+self.solps_output_file[1] #0=filename (solpsOut) ; 1=fomat (txt) 
                #solps_outFile=self.SOLPS_OUTPUT_FILE[0]+'_'+str(i)+'.'+self.SOLPS_OUTPUT_FILE[1] #0=filename (solpsOut) ; 1=fomat (txt)
                self.solpsParams=param_handler.read(self.INPUT_DIR+'/'+solps_outFile)

                print('\t reading output of SOLPS:')
                for k,v, in self.solpsParams.items():
                    print(('\t {0} : {1}'.format(k, v)))
                print(' ')
                sys.stdout.flush()

                #format float to list if needed
                #currently implemented for inputEnergy, Ti, Te, inputAngle, bfieldAngle
                #might need to add Z (all same values or add 1.0?)

                #idea is: if solpsParams[k] is a single value (float or int),
                #turn into list with one Ein value per species (all same value)
                #otherwise (if not float/int, no plasmaSpecies given, etc), leave as is
                print('\t Format Ein, Ti, Te, Ain and/or Bin as lists when posible:')
                if ('inputEnergy' in self.solpsParams):
                    Ein=[]
                    if (isinstance(self.solpsParams['inputEnergy'],float) or isinstance(self.solpsParams['inputEnergy'],int)):
                        if ('plasmaSpecies' in self.solpsParams):                            
                            for n in range(len(self.solpsParams['plasmaSpecies'])):
                                Ein.append(self.solpsParams['inputEnergy'])
                        else: #without plasmaSpecies, just leave as is
                            Ein=self.solpsParams['inputEnergy']
                    else: #if it's not float or int, leave as is (whether it's a list or not)
                        Ein=self.solpsParams['inputEnergy']
                    print('\t Ein =', Ein)
                
                if ('Ti' in self.solpsParams):
                    Ti=[]
                    if (isinstance(self.solpsParams['Ti'],float) or isinstance(self.solpsParams['Ti'],int)):
                        if ('plasmaSpecies' in self.solpsParams):
                            for n in range(len(self.solpsParams['plasmaSpecies'])):
                                Ti.append(self.solpsParams['Ti'])
                        else: #without plasmaSpecies, just leave as is
                            Ti=self.solpsParams['Ti']
                    else: #if it's not float or int, leave as is (whether it's a list or not)
                        Ti=self.solpsParams['Ti']
                    print('\t Ti =', Ti)
                        
                if ('Te' in self.solpsParams):
                    Te=[]
                    if (isinstance(self.solpsParams['Te'],float) or isinstance(self.solpsParams['Te'],int)):
                        if ('plasmaSpecies' in self.solpsParams):
                            for n in range(len(self.solpsParams['plasmaSpecies'])):
                                Te.append(self.solpsParams['Te'])
                        else: #without plasmaSpecies, just leave as is
                            Te=self.solpsParams['Te']
                    else: #if it's not float or int, leave as is (whether it's a list or not)
                        Te=self.solpsParams['Te']
                    print('\t Te =', Te)
                    
                if ('inputAngle' in self.solpsParams):
                    Ain=[]
                    if (isinstance(self.solpsParams['inputAngle'],float) or isinstance(self.solpsParams['inputAngle'],int)):
                        if ('plasmaSpecies' in self.solpsParams):
                            for n in range(len(self.solpsParams['plasmaSpecies'])):
                                Ain.append(self.solpsParams['inputAngle'])
                        else: #without plasmaSpecies, just leave as is
                            Ain=self.solpsParams['inputAngle']
                    else: #if it's not float or int, leave as is (whether it's a list or not)
                        Ain=self.solpsParams['inputAngle']
                    print('\t Ain =', Ain)
                    
                if ('bfieldAngle' in self.solpsParams):
                    Bin=[]
                    if (isinstance(self.solpsParams['bfieldAngle'],float) or isinstance(self.solpsParams['bfieldAngle'],int)):
                        if ('plasmaSpecies' in self.solpsParams):
                            for n in range(len(self.solpsParams['plasmaSpecies'])):
                                Bin.append(self.solpsParams['bfieldAngle'])
                        else: #without plasmaSpecies, just leave as is
                            Bin=self.solpsParams['bfieldAngle']
                    else: #if it's not float or int, leave as is (whether it's a list or not)
                        Bin=self.solpsParams['bfieldAngle']
                    print('\t Bin =', Bin)
                print('\n')
                    
                ## 2 - do ops to calculate/format inputs for FTX:
                print ('\t REFORMAT SOLPS output --> FTX input')
                self.ftxInputs={}            

                print('\t TEST: solps output dict = ', self.solpsParams)
                sys.stdout.flush()

                print('\n')
                
                ## 2.a : check for parameters expected from SOLPS & 'translate' accordingly

                ### 2.a.i: parameters with a single value
                
                ### particle flux: /m2s --> /nm2s done in ftx driver
                if ('flux' in self.solpsParams):
                    self.ftxInputs['flux'] = self.solpsParams['flux'] #/1.0e18
                    print('\t TEST: ftx[flux] = ', self.ftxInputs['flux'])
                sys.stdout.flush()
                    
                ### heat flux: W/m2 --> W/nm2
                ### in some (older) version of Xolotl, heat is one line; in others is two lines
                ### for now implement single version suited for xolotl-tempGrid-build (print warning)
                ### might need to make changes to Xolotl param template too (tbd)
                if 'heatFlux' in self.solpsParams:
                    print ('\t \t WARNING: heat flux in Xolotl handled in the newer way, using 2 lines: temperature model & values')
                    self.ftxInputs['tempHandler'] = 'heat'
                    #if solps heatflux contains 2 fields, assume it's bulkT given, then assume room T
                    if (len(self.solpsParams['heatFlux']) == 2):  
                        self.ftxInputs['tempParam'] = (self.solpsParams['heatFlux'][0], self.solpsParams['heatFlux'][1]) #/m2s --> /nm2s done in ftx driver 
                    #if a single field given, assume bulkT = room T
                    elif (len(self.solpsParams['heatFlux']) == 1):
                        self.ftxInputs['tempParam'] = (self.solpsParams['heatFlux'], 300.0)  #/m2s --> /nm2s done in ftx driver 
                        print ('\t \t WARNING: no bulkT given; assume room T = 300 K')
                    print('\t \t TEST: heat flux tempHandler =', self.ftxInputs['tempHandler'], 'and temparam =', self.ftxInputs['tempParam'])
                    print('\n')
                sys.stdout.flush()


                ### 2.a.ii: inputs that are given for each species:
                    # flux fraction
                    # input energy
                    # input angle
                        
                ### for now, include all
                if ('plasmaSpecies' in  self.solpsParams): 
                    #self.ftxInputs['plasmaSpecies'] = 'He W D T C'
                    #C_f : different options for C in F-TRIDYN: _a, _b, _d... ; _f noted as "C for fusionn"
                    self.ftxInputs['plasmaSpecies'] = ['He','W', 'D', 'T', 'C_f'] 
                    print('\t hard coded plasma species in ftx = ', self.ftxInputs['plasmaSpecies'] )
                    print('\n')
                    
                    #create a variable that's 1 if species exists ; 0 if it doesn't:
                    #probably won't need it, but leave it commented out for now:
                    #for n in range)len(self.ftxInputs['plasmaSpecies'])):
                    #    s = self.ftxInputs['plasmaSpecies'][n]
                    #    if (s in self.solpsParams['plasmaSpecies']):
                    #        speciesExists[n]=1.0
                    #    else:
                    #        speciesExists[n]=0.0
                    
                
                ### reformat to have a vaolue for all species, even if 0.0
                ### might not be the most efficient method, but it should work for all scenarios
                ### example of what I'm trying to do in loop below
                   #if ('He' in self.solpsParams['plasmaSpecies']):
                   #    ftxInputs['fluxFraction'][0] = self.solpsParams['fluxFraction'][nCount-1]
		   #    nCount+=1
                   #else:
                   #    ftxInputs['fluxFraction'][0] = 0.0
                   #
                   #if ('W' in self.solpsParams['plasmaSpecies']):
                   #    ftxInputs['fluxFraction'][1] = self.solpsParams['fluxFraction'][nCount-1]
                   #else:
                   #    ftxInputs['fluxFraction'][1] = 0.0

                print('\t reformat to have Ein, Ain for each species:')
                nCount=0
                for n in range(len(self.ftxInputs['plasmaSpecies'])):
                    s = self.ftxInputs['plasmaSpecies'][n]                    
                    print('\t for n  = ', n , ' ; s = ', s, ' :')

                    #initialize dictionary entries:
                    #create dictionary entries for the first species
                    if (n==0):
                        print('\t \t TEST: initialize Ein, Ti, Te, Ain, Bin lists')
                        #always write flux fraction, even if zero:
                        self.ftxInputs['fluxFraction']=[] #[0.0] DEL
                        print('\t \t TEST: as of loop for ', s, 'ftxInputs[fluxFraction] = ', self.ftxInputs['fluxFraction'])
                        #write inputEnergy if input provided by SOLPS:
                        if  (('inputEnergy' in self.solpsParams) or ('Ti' in self.solpsParams) or ('Te' in self.solpsParams)):
                            self.ftxInputs['inputEnergy']=[] #[0.0] DEL
                            print('\t \t TEST: as of loop for ', s, 'ftxInputs[inputEnergy] = ', self.ftxInputs['inputEnergy'])
                        #write inputAngle if input provided by SOLPS:
                        if (('inputAngle' in self.solpsParams) or ('bfieldAngle' in self.solpsParams)):
                            self.ftxInputs['inputAngle']=[] #[0.0] DEL
                            print('\t \t TEST: as of loop for ', s, 'ftxInputs[inputAngle] = ', self.ftxInputs['inputAngle'])
                        sys.stdout.flush()
                        print('\n')
                        
                    #append values to dictionary entry in subsequent species
                    #DELETE THIS ENTIRE SECTION?
                    #else:                    
                    #    #always write flux fraction, even if zero:
                    #    self.ftxInputs['fluxFraction'].append(0.0) 
                    #    print('TEST: as of loop for ', s, 'ftxInputs[fluxFraction] = ', self.ftxInputs['fluxFraction']) 
                    #    #write inputEnergy if input provided by SOLPS:
                    #    if  (('inputEnergy' in self.solpsParams) or ('Ti' in self.solpsParams) or ('Te' in self.solpsParams)):
                    #        self.ftxInputs['inputEnergy'].append(0.0)
                    #        print('TEST: as of loop for ', s, 'ftxInputs[inputEnergy] = ', self.ftxInputs['inputEnergy']) 
                    #    #write inputAngle if input provided by SOLPS:
                    #    if (('inputAngle' in self.solpsParams) or ('bfieldAngle' in self.solpsParams)):
                    #        self.ftxInputs['inputAngle'].append(0.0)
                    #        print('TEST: as of loop for ', s, 'ftxInputs[inputAngle] = ', self.ftxInputs['inputAngle'])
                    #    sys.stdout.flush()
                
                    if (s in self.solpsParams['plasmaSpecies']):
                        nCount+=1
                        m=nCount-1
                        #flux fraction
                        print('\t check for fluxFraction:')
                        if ('fluxFraction' in self.solpsParams):
                            self.ftxInputs['fluxFraction'].append(self.solpsParams['fluxFraction'][m]) #[n] = self.solpsParams['fluxFraction'][m] DEL                            
                        else:
                            print('\t \t WARNING: no fluxFraction provided by input file. set all to zero')
                            self.ftxInputs['fluxFraction'].append(0.0) #[n] = 0.0 DEL
                        print('\t \t ftx[fluxFraction] = ', self.ftxInputs['fluxFraction'])
                        print('\n')
                        sys.stdout.flush()
                        
                        #inputEnergy
                        #check whether passed as energy, or need to estimate based on Te & Ti
                        print('\t check for inputEnergy:')
                        if  ('inputEnergy' in self.solpsParams):
                            print('\t \t TEST: solpsParams[inputEnergy][m] = ', Ein[m]) # self.solpsParams['inputEnergy'][m]) DEL
                            sys.stdout.flush() #TEST
                            self.ftxInputs['inputEnergy'].append(Ein[m]) #[n] = self.solpsParams['inputEnergy'][m] DEL
                            print('\t \t ftx[inputEnergy] = ', self.ftxInputs['inputEnergy'])
                            print('\n')
                        elif ('Ti' in self.solpsParams):
                            if ('Te' in self.solpsParams):
                                if ('Z' in self.solpsParams):
                                    print('\t \t WARNING: no inputEnegy given ; estimate as 2Ti + 3ZTe')                                
                                    print('\t \t TEST: solpsParams[Ti] = ', Ti[m]) #self.solpsParams['Ti'][m]) DEL
                                    print('\t \t TEST: solpsParams[Z] = ', self.solpsParams['Z'][m]) #may need to change to Z[m]
                                    print('\t \t TEST: solpsParams[Te] = ', Te[m]) #self.solpsParams['Te'][m]) DEL
                                    sys.stdout.flush() #TEST DEL
                                    self.ftxInputs['inputEnergy'].append(2*Ti[m]+3*self.solpsParams['Z'][m]*Te[m]) 
                                    #[n] = 2*self.solpsParams['Ti'][m]+3*self.solpsParams['Z'][m]*self.solpsParams['Te'][m]
                                else:
                                    print('\t \t WARNING: no inputEnegy given or Z given ; estimate as 2Ti + 3Te (Z=1)')
                                    print('\t \t TEST: solpsParams[Ti] = ', Ti[m]) #self.solpsParams['Ti'][m]) DEL
                                    print('\t \t TEST: solpsParams[Te] = ', Te[m]) #self.solpsParams['Te'][m]) DEL
                                    sys.stdout.flush() #TEST                                                                          
                                    self.ftxInputs['inputEnergy'].append(2*Ti[m]+3*Te[m])
                                    #[n] = 2*self.solpsParams['Ti'][m]+3*sself.solpsParams['Te'][m]
                            elif ('Te' not in self.solpsParams):
                                if ('Z' in self.solpsParams):
                                    print('\t \t WARNING: no inputEnegy or Te given ; estimate as Ti * (2+3Z)')
                                    print('\t \t TEST: solpsParams[Ti] = ', Ti[m]) #self.solpsParams['Ti'][m]) DEL
                                    print('\t \t TEST: solpsParams[Z] = ', self.solpsParams['Z'][m]) #may need to change to Z[m]  
                                    sys.stdout.flush() #TEST                                                                          
                                    self.ftxInputs['inputEnergy'].append(Ti*(2+3*self.solpsParams['Z'][m]))
                                    #[n] = self.solpsParams['Ti'][m]*(2+3*self.solpsParams['Z'][m])
                                else:
                                    print('\t \t WARNING: no inputEnegy, Te or Z given ; estimate as 5*Ti (Z=1)')
                                    print('\t \t TEST: solpsParams[Ti] = ', Ti[m]) #self.solpsParams['Ti'][m]) DEL
                                    sys.stdout.flush() #TEST                                                                          
                                    self.ftxInputs['inputEnergy'].append(5*Ti[m])
                                    #[n] = 5*self.solpsParams['Ti'][m]
                            print('\t \t using Ti-Te, ftx[inputEnergy] = ', self.ftxInputs['inputEnergy'])
                            print('\n')
                        elif ('Te' in self.solpsParams):
                            if ('Z' in self.solpsParams):
                                print('\t \t WARNING: no inputEnegy or Ti given ; estimate as Te * (2+3Z)')
                                print('\t \t TEST: solpsParams[Z] = ', self.solpsParams['Z'][m]) #may need to change to Z[m]  
                                print('\t \t TEST: solpsParams[Te] = ', Te[m]) #self.solpsParams['Te'][m]) DEL
                                sys.stdout.flush() #TEST                                                                          
                                self.ftxInputs['inputEnergy'].append(Te*(2+3*self.solpsParams['Z'][m]))
                                #[n] = self.solpsParams['Te'][m]*(2+3*self.solpsParams['Z'][m])
                            else:
                                print('\t \t WARNING: no inputEnegy, Ti or Z given ; estimate as 5*Te (Z=1)')
                                print('\t \t TEST: solpsParams[Te] = ', self.solpsParams['Te'][m])
                                sys.stdout.flush() #TEST                                                                          
                                self.ftxInputs['inputEnergy'].append(5*Te[m])
                                #[n] = 5*self.solpsParams['Te'][m]
                            print('\t \t using Te, ftx[inputEnergy] = ', self.ftxInputs['inputEnergy'])
                            print('\n')
                        else:
                            print ('\t \t WARNING: no inputEnergy, Te or Ti given ; continue with default values in config file')
                            print('\n')
                        sys.stdout.flush()
                        
                        #inputAngle
                        print('\t check for inputAngle:')
                        if ('inputAngle' in self.solpsParams):
                            self.ftxInputs['inputAngle'].append(Ain[m])
                            #[n] = self.solpsParams['inputAngle'][m]
                            print('\t \t ftx[inputAngle] = ', self.ftxInputs['inputAngle'])
                            print('\n')
                        elif ('bfieldAngle' in self.solpsParams):
                            print ('\t \t WARNING: no inputAngle given; use B field angle ')
                            self.ftxInputs['inputAngle'].append(Bin[m])
                            print('\t \t ftx[inputAngle] = ', self.ftxInputs['inputAngle'])
                            print('\n')
                            #[n] = self.solpsParams['bfieldAngle'][m]                            
                        else:
                            print('\t \t WARNING: no inputAngle or bfieldAngle given; continue with default values in config file (likely normal incidence)')
                            print('\n')
                        sys.stdout.flush()    

                    else: # if s not in solpsParams['plasmaSpecies']
                        print('\t ', s, ' is not in solps[plasmaSpecies] ; set fluxFraction, Ein, Ain to zero')
                        #fluxFraction
                        if ('fluxFraction' in self.solpsParams):
                            self.ftxInputs['fluxFraction'].append(0.0) #[n] = 0.0
                            print('\t \t fluxFraction = 0.0')
                        else:
                            print('\t \t WARNING: no fluxFraction provided by input file. set all to zero')
                            self.ftxInputs['fluxFraction'].append(0.0) #[n] = 0.0
                        #inputEnergy
                        if  (('inputEnergy' in self.solpsParams) or ('Ti' in self.solpsParams) or ('Te' in self.solpsParams)):
                            self.ftxInputs['inputEnergy'].append(0.0) #[n] = 0.0
                            print('\t \t inputEnergy = 0.0')
                        else:
                            print ('\t \t WARNING: no inputEnergy, Te or Ti given ; continue with default values in config file')
                        #inputAngle
                        if (('inputAngle' in self.solpsParams) or ('bfieldAngle' in self.solpsParams)):
                            self.ftxInputs['inputAngle'].append(0.0) #[n] = 0.0
                            print('\t \t inputAngle = 0.0')
                        else:
                            print('\t \t WARNING: no inputAngle or bfieldAngle given; continue with default values in config file (likely normal incidence)')
                        print('\n')
                        sys.stdout.flush()
                        
                ## 2.b : write down all other parameters 'as is':

                print('\t write down all other parameters "as is":')
                for k,v, in self.solpsParams.items():

                    if k in self.ftxInputs:
                        print('\t \t TEST: k = ', k)
                        print('\t \t TEST: v = ', v)
                        sys.stdout.flush()
                        print('\t \t ftx[',k,'] = ', self.ftxInputs[k])
                    else:
                        print('\t \t param ', k , 'not in ftx dictionary yet ; add as is')
                        print('\t \t TEST: v = ', v)
                        self.ftxInputs[k] = v
                print('\n')
                sys.stdout.flush()
                
                # 3 - print ftxIn.txt - DONE
                print('\t write FTX input file:')
                ftxInFileFormat=list(self.services.get_config_param('FTX_INPUT_FORMAT'))
                ftxInFileName=ftxInFileFormat[0]+str(i)+'.'+ftxInFileFormat[1]
                print('\t \t file name: ', ftxInFileName)
                ftxInFile=open(ftxInFileName, "w")
                for k,v, in self.ftxInputs.items():
                    if (isinstance(v, int) or isinstance(v, float) or isinstance(v, dict) or isinstance(v, str)):
                        print(('\t \t {0} : {1}'.format(k, v)))
                        ftxInFile.write(('{0}={1}\n'.format(k, v)))
                    else: #not float, int or string; assume it's list
                        print('type of v is ', type(v))
                        sys.stdout.flush()
                        v_string=''
                        for i in range(len(v)):
                            v_string+=str(v[i])+' '
                        v_string = v_string[:-1]
                        print(('\t \t {0} : {1}'.format(k, v_string)))
                        ftxInFile.write(('{0}={1}\n'.format(k, v_string)))
                        del v_string
                ftxInFile.close()
                print(' ')
                
                    
                # 4 - add timeStap to its own file
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
                #ftxInFilePath=self.SUBMIT_DIR+'/'+ftxInFileName
                print('\t copying input to FTX file from ', ftxInFileName, 'to ', self.child_components[child_comp]['INPUT_DIR']+'/ftxInput.txt')
                shutil.copyfile(ftxInFileName,self.child_components[child_comp]['INPUT_DIR']+'/ftxInput.txt')
                print('\t copying time parameter file from ', timeFileName, 'to ', self.child_components[child_comp]['INPUT_DIR']+'/'+timeFileName)
                shutil.copyfile(timeFileName,self.child_components[child_comp]['INPUT_DIR']+'/'+timeFileName)
                print('\t copying output file of SOLPS from ', self.INPUT_DIR+'/'+solps_outFile, 'to ', self.child_components[child_comp]['INPUT_DIR']+'/solpsOut.txt')
                shutil.copyfile(self.INPUT_DIR+'/'+solps_outFile,self.child_components[child_comp]['INPUT_DIR']+'/solpsOut.txt')
                
                print('\n')
                sys.stdout.flush()

                #save ftxIn and solpsOut for each loop (with time-stamp)
                shutil.copyfile(self.INPUT_DIR+'/'+solps_outFile, solps_outFile+'_t'+str(timeStamp))
                shutil.copyfile(ftxInFileName,ftxInFileName+'_t'+str(timeStamp))
                shutil.copyfile(timeFileName,timeFileName+'_t'+str(timeStamp))
                print('\t to save input files for each time-loop, copied:')
                print('\t \t ', self.INPUT_DIR+'/'+solps_outFile, ' as ', solps_outFile+'_t'+str(timeStamp))
                print('\t \t ', ftxInFileName, ' as ', ftxInFileName+'_t'+str(timeStamp))
                print('\t \t ', timeFileName , ' as ', timeFileName+'_t'+str(timeStamp))

                
                #update plasma state file:
                self.services.update_plasma_state()
                
            ## END OF solpsOut --> ftxIn and time-parameters
            
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
