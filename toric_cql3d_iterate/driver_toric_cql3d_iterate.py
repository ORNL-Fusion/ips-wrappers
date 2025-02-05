#! /usr/bin/env python
""" Version 1 (Batchelor 7/23/2017)
Driver for toric/cql3d iteration with logic exposed here rather than buried in the 
 components.  Adapted from generic_driver.py
"""
from __future__ import print_function

import sys
import os
import subprocess
import getopt
import shutil
import math
import toric_reader
import h5py
from ipsframework import Component
import numpy as np
from netCDF4 import *

class toric_driver(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timestamp=0):
      # Driver initialization ? nothing to be done
        return

# ------------------------------------------------------------------------------
#
# step function
#
# ------------------------------------------------------------------------------

    def step(self, timestamp=0):

        services = self.services
        services.stage_state()
        services.stage_input_files(self.INPUT_FILES)

      # get list of ports
        ports = services.get_config_param('PORTS')
        port_names = ports['NAMES'].split()
        print('PORTS =', port_names)
        port_dict = {}
        port_id_list = []

      # Instantiate components in port_names list, except DRIVER itself
        print (' ')

        # INIT is already instantiated by the framework. This adds it to port_dict
        if 'INIT' in port_names:
            initComp = services.get_port('INIT') 
            if(initComp == None):
                print('Error accessing INIT component')
                raise
            port_dict['INIT'] = initComp
            port_id_list.append(initComp)
            print (' ')

        if 'EPA' in port_names:
            epaComp = services.get_port('EPA')
            if(epaComp == None):
                print('Error accessing EPA component')
                raise
            port_dict['EPA'] = epaComp
            port_id_list.append(epaComp)
            print (' ')
       
        if 'RF_EC' in port_names:
            rf_ecComp = services.get_port('RF_EC')
            if(rf_ecComp == None):
                print('Error accessing RF_EC component')
                raise
            port_dict['RF_EC'] = rf_ecComp
            port_id_list.append(rf_ecComp)
            print (' ')
        
        if 'RF_IC' in port_names:
            rf_icComp = services.get_port('RF_IC')
            if(rf_icComp == None):
                print('Error accessing RF_IC component')
                raise
            port_dict['RF_IC'] = rf_icComp
            port_id_list.append(rf_icComp)
            print (' ')
        
        if 'RF_LH' in port_names:
            rf_lhComp = services.get_port('RF_LH')
            if(rf_lhComp == None):
                print('Error accessing RF_LH component')
                raise
            port_dict['RF_LH'] = rf_lhComp
            port_id_list.append(rf_lhComp)
            print (' ')

        if 'NB' in port_names:
            nbComp = services.get_port('NB')
            if(nbComp == None):
                print('Error accessing NB component')
                raise
            port_dict['NB'] = nbComp
            port_id_list.append(nbComp)
            print (' ')

        if 'FUS' in port_names:
            fusComp = services.get_port('FUS')
            if(fusComp == None):
                print('Error accessing FUS component')
                raise
            port_dict['FUS'] = fusComp
            port_id_list.append(fusComp)
            print (' ')
 
        if 'FP' in port_names:
            fpComp = services.get_port('FP')
            if(fpComp == None):
                print('Error accessing FP component')
                raise
            port_dict['FP'] = fpComp
            port_id_list.append(fpComp)
            print (' ')
        
        if 'MONITOR' in port_names:
            monitorComp = services.get_port('MONITOR')
            if(monitorComp == None):
                print('Error accessing MONITOR component')
                raise       
            port_dict['MONITOR'] = monitorComp
            port_id_list.append(monitorComp)
            print (' ')

      # Is this a simulation startup or restart
        sim_mode = services.get_config_param('SIMULATION_MODE')
        print('SIMULATION_MODE =', sim_mode)

      # Get timeloop for simulation
        timeloop = services.get_time_loop()
        tlist_str = ['%.3f'%t for t in timeloop]
        t = tlist_str[0]

      # Initialize components in PORTS list for startup or restart
        print (' ')
        
        init_mode = 'init'
        if sim_mode == 'RESTART' : init_mode = 'restart'

        if 'EPA' in port_names:
            self.component_call(services, 'EPA', epaComp, init_mode, t)
        
        if 'RF_EC' in port_names:
            self.component_call(services, 'RF_EC', rf_ecComp, init_mode, t)
        
        if 'RF_IC' in port_names:
            self.component_call(services, 'RF_IC', rf_icComp, init_mode, t)
        
        if 'RF_LH' in port_names:
            self.component_call(services, 'RF_LH', rf_lhComp, init_mode, t)

        if 'NB' in port_names:
            self.component_call(services, 'NB', nbComp, init_mode, t)

        if 'FUS' in port_names:
            self.component_call(services, 'FUS', fusComp, init_mode, t)

        if 'FP' in port_names:
            self.component_call(services, 'FP', fpComp, init_mode, t)

        if 'MONITOR' in port_names:
            self.component_call(services, 'MONITOR', monitorComp, init_mode, t)

        # Get plasma state files into driver work directory and copy to psn if there is one
        services.stage_state()
        cur_state_file = services.get_config_param('CURRENT_STATE')
        
        # Get Portal RUNID and save to a file
        run_id = services.get_config_param('PORTAL_RUNID')
        sym_root = services.get_config_param('SIM_ROOT')
        path = os.path.join(sym_root, 'PORTAL_RUNID')
        specs = services.get_config_param('SPECS')
        if specs in ['CUSTOM','Custom','custom']:
            custom_specs_list = self.get_config_param(services,'CUSTOM_SPECS').split(' ')
            custom_specs = ", ".join(custom_specs_list)
            N_gen = len(custom_specs_list)

        rfpwr_arg = services.get_config_param('RFPWR_IC')

        # If pwrscale is in config parameters use power rescaling
        arg_pwrscale = self.get_config_param(services,'PWRSCALE',optional= True)
        if arg_pwrscale!=None:
            if arg_pwrscale.strip() in [True, 'true', 'True', 'TRUE']: 
                pwrscale_on = True
                #create powerscale hdf5 file now uses ngena value to set maximum
                #number of values in pwrscale hdf5
                ngena = 4
                pwrscale = np.full(ngena,1.0)
                pfrac_toric = np.full(ngena,1.0)
                pfrac_cql3d = np.full(ngena,1.0)
                pwrtarget = float(rfpwr_arg)
                pwrscale_f =h5py.File('pwrscale.hdf5','w')
                pwrtarget_dset = pwrscale_f.create_dataset("pwrtarget",data=pwrtarget)
                pscale_dset = pwrscale_f.create_dataset("pwrscale",data=pwrscale)
                pfrac_toric_dset = pwrscale_f.create_dataset("pfrac_toric",data=pfrac_toric)
                pfrac_cql3d_dset = pwrscale_f.create_dataset("pfrac_cql3d",data=pfrac_cql3d)
                pwrscale_f.close()
                # Update plasma state files in plasma_state work directory
                try:
                    services.update_state()
                except Exception:
                    logMsg = 'Error in call to update_plasma_state()'
                    self.services.exception(logMsg)
                    raise
            else:
                pwrscale_on = False
        else:
            pwrscale_on=False
        
        # Post init processing: stage plasma state, stage output
        services.stage_output_files(t, self.OUTPUT_FILES)

        print(' init sequence complete')

        INIT_ONLY = self.get_component_param(services, 'INIT_ONLY', optional = True)
        if INIT_ONLY in [True, 'true', 'True', 'TRUE']:   
            message = 'INIT_ONLY: Intentional stop after INIT phase'
            print(message)
            return

# ------------------------------------------------------------------------------
#
# Start Physics Layer
#
# ------------------------------------------------------------------------------


        print('Zeroth step - Run toric only in toric mode and qldci mode for Maxwellian')
        if sim_mode == 'NORMAL' :   # i.e. not RESTART do Maxwellian runs
            t = tlist_str[0]
            print (' ')
            print('\nDriver: starting iteration ', t)
            services.update_time_stamp(t)

            try:
                services.call(rf_icComp, 'step', t, toric_Mode = 'toric', \
                              inumin_Mode = 'Maxwell' , isol_Mode = '1', \
                              save_output = 'False', pwrscale_on=pwrscale_on)
            except Exception:
                message = 'RF_IC toric mode step failed'
                print(message)
                services.exception(message)
                raise 
   
            if(specs=='MIN'):
                try:
                    services.call(rf_icComp, 'step', t, toric_Mode = 'qldci', \
                                  inumin_Mode = 'Maxwell' , isol_Mode = '1', \
                                  save_output = 'True')
                except Exception:
                    message = 'RF_IC qldci mode step failed'
                    print(message)
                    services.exception(message)
                    raise

            if(specs=='MIN+'):
                try:
                    services.call(rf_icComp, 'step', t, toric_Mode = 'qldci1', \
                                  inumin_Mode = 'Maxwell' , isol_Mode = '1', \
                                  save_output = 'False')
                except Exception:
                    message = 'RF_IC qldci mode step failed'
                    print(message)
                    services.exception(message)
                    raise
                try:
                    services.call(rf_icComp, 'step', t, toric_Mode = 'qldci2', \
                                  inumin_Mode = 'Maxwell' , isol_Mode = '1', \
                                  save_output = 'True')
                except Exception:
                    message = 'RF_LH qldci mode step failed'
                    print(message)
                    services.exception(message)
                    raise
            if(specs=='CUSTOM'):
                for i in range(N_gen):
                    toricMode = 'qldci'+str(i+1)
                    if (i+1) == N_gen:
                        save = True
                    else:
                        save = False
                    try:
                        services.call(rf_icComp, 'step', t, toric_Mode = toricMode, \
                                      inumin_Mode = 'Maxwell', isol_Mode='1', \
                                      save_output=save)
                    except Exception:
                        message = 'RF_LH qldci mode step failed'
                        print(message)
                        services.exception(message)
                        raise
                        

        # Iterate through the timeloop, or in this case iteration loop
        for t in tlist_str[1:len(timeloop)]:
            print('\nDriver: starting iteration ', t)
            services.update_time_stamp(t)

            # call pre_step_logic
            services.stage_state()
            self.pre_step_logic(float(t))
            services.update_state()

            # Call step for each component
            if 'EPA' in port_names:
                self.component_call(services, 'EPA', epaComp, 'step', t)

            # Call powerscale iteration
            icount = 0
            running = True
            if pwrscale_on:
                while running:
                    try:
                        services.call(fpComp, 'step', t,rf_code='toric',pwrscale_on=pwrscale_on,icount=icount)
                    except Exception:
                        message = 'FP step failed'
                        print(message)
                        services.exception(message)
                        raise
                    services.stage_state() #get katest state w/ cql powers
                    pwrscale_f = h5py.File('pwrscale.hdf5','r+')
                    pwrtarget_dset = pwrscale_f['pwrtarget']
                    pwrscale_dset = pwrscale_f['pwrscale']
                    pwrscale_toric_dset = pwrscale_f['pfrac_toric']
                    pwrscale_cql3d_dset = pwrscale_f['pfrac_cql3d']
                    pwrscale_min_ratio = pwrscale_toric_dset[0]/pwrscale_cql3d_dset[0]
                    pwrscale_bulk_ratio = pwrscale_toric_dset[1]/pwrscale_cql3d_dset[1]
                    print('pwrscale_min_ratio', pwrscale_min_ratio)
                    print('pwrscale_bulk_ratio', pwrscale_bulk_ratio)
                    
                    if ((pwrscale_min_ratio > 1.05)or(pwrscale_min_ratio<0.95)) or \
                       ((pwrscale_bulk_ratio > 1.05)or(pwrscale_bulk_ratio<0.95)):
                        pwrscale_dset[0] = pwrscale_min_ratio*pwrscale_dset[0]
                        pwrscale_dset[1] = pwrscale_bulk_ratio*pwrscale_dset[1]
                    else:
                        running = False
                    print('pwrscale iter: ', icount, ' pwrscale ', np.array(pwrscale_dset))
                    print('toric_pfrac', np.array(pwrscale_toric_dset))
                    print('cql3d_pfrac', np.array(pwrscale_cql3d_dset))
                    pwrscale_f.close()
                    services.update_state() #need to do a state update here
                    icount = icount+1
            else:
                try:
                    services.call(fpComp, 'step', t,rf_code='toric',pwrscale_on=pwrscale_on)
                except Exception:
                    message = 'FP step failed'
                    print(message)
                    services.exception(message)
                    raise 
            try:
                services.call(rf_icComp, 'step', t, toric_Mode = 'toric', \
                inumin_Mode = 'nonMaxwell' , isol_Mode = '1', \
                                  save_output = 'False',pwrscale_on=pwrscale_on)
            except Exception:
                message = 'RF_IC toric mode step failed'
                print(message)
                services.exception(message)
                raise 

            if(specs=='MIN'):
                try:
                    services.call(rf_icComp, 'step', t, toric_Mode = 'qldci', \
                                  inumin_Mode = 'nonMaxwell' , isol_Mode = '1', \
                                  save_output = 'True')
                except Exception:
                    message = 'RF_IC qldci mode step failed'
                    print(message)
                    services.exception(message)
                    raise

            if(specs=='MIN+'):
                try:
                    services.call(rf_icComp, 'step', t, toric_Mode = 'qldci1', \
                                  inumin_Mode = 'nonMaxwell' , isol_Mode = '1', \
                                  save_output = 'False')
                except Exception:
                    message = 'RF_IC qldci mode step failed'
                    print(message)
                    services.exception(message)
                    raise
                try:
                    services.call(rf_icComp, 'step', t, toric_Mode = 'qldci2', \
                                  inumin_Mode = 'nonMaxwell' , isol_Mode = '1', \
                                  save_output = 'True')
                except Exception:
                    message = 'RF_IC qldci mode step failed'
                    print(message)
                    services.exception(message)
                    raise

            if(specs=='CUSTOM'):
                for i in range(N_gen):
                    toricMode = 'qldci'+str(i+1)
                    if (i+1) == N_gen:
                        save = True
                    else:
                        save = False
                    try:
                        services.call(rf_icComp, 'step', t, toric_Mode = toricMode, \
                                      inumin_Mode = 'Maxwell', isol_Mode='1', \
                                      save_output=save)
                    except Exception:
                        message = 'RF_LH qldci mode step failed'
                        print(message)
                        services.exception(message)
                        raise
            
            if 'MONITOR' in port_names:
                self.component_call(services, 'MONITOR', monitorComp, 'step', t)

            services.stage_state()

         # Post step processing: stage plasma state, checkpoint components and self
            services.stage_output_files(t, self.OUTPUT_FILES)
            services.checkpoint_components(port_id_list, t)
            self.checkpoint(t)

# ------------------------------------------------------------------------------
#
# End of Physics Layer
#
# ------------------------------------------------------------------------------

      # Post simulation: call checkpoint on each component

        services.checkpoint_components(port_id_list, t, Force = True)
        self.checkpoint(t)
      
      # Post simulation: call finalize on each component
        print (' ')
        if 'FP' in port_names:
            self.component_call(services, 'FP', fpComp, 'finalize', t)

        if 'RF_EC' in port_names:
            self.component_call(services, 'RF_EC', rf_ecComp, 'finalize', t)

        if 'RF_IC' in port_names:
            self.component_call(services, 'RF_IC', rf_icComp, 'finalize', t)

        if 'RF_LH' in port_names:
            self.component_call(services, 'RF_LH', rf_lhComp, 'finalize', t)

        if 'NB' in port_names:
            self.component_call(services, 'NB', nbComp, 'finalize', t)

        if 'FUS' in port_names:
            self.component_call(services, 'FUS', fusComp, 'finalize', t)

        if 'EPA' in port_names:
            self.component_call(services, 'EPA', epaComp, 'finalize', t)

        if 'MONITOR' in port_names:
            self.component_call(services, 'MONITOR', monitorComp, 'finalize', t)

# ------------------------------------------------------------------------------
#
# checkpoint function
#
# Does nothing
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print('toric_driver.checkpoint() called')
        

# ------------------------------------------------------------------------------
#
# finalize function 
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp = 0):
      # Driver finalize - nothing to be done
        pass

# "Private" driver methods

    # Component call - wraps the exception handling for all component calls
    def component_call(self, services, port_name, comp, mode, time):
            comp_mode_string = port_name + ' ' + mode
            print (comp_mode_string)

            try:
                services.call(comp, mode, time)
            except Exception:
                message = comp_mode_string + ' failed'
                print(message)
                services.exception(message)
                raise 
            
            return
    
    
    # Pre Step Logic
    def pre_step_logic(self, timeStamp):

        cur_state_file = self.services.get_config_param('CURRENT_STATE')

        ps = Dataset(cur_state_file, 'r+', format = 'NETCDF3_CLASSIC')
        t1 = ps.variables['t1'].getValue()
        ps.variables['t0'].assignValue(t1)
        ps.variables['t1'].assignValue(timeStamp)
        self.services.log('ps%t1 = ', t1)

        power_ic = 0.0
        if ('power_ic' in list(ps.variables.keys())):
            power_ic = ps.variables['power_ic'][:]
            print('toric_driver pre_step_logic: power_ic = ', power_ic)

        ps.close()
        
        print('toric_driver pre_step_logic: timeStamp = ', timeStamp)
        
        return

# ------------------------------------------------------------------------------
#
# "Private"  methods
#
# ------------------------------------------------------------------------------

    # Try to get config parameter - wraps the exception handling for get_config_parameter()
    def get_config_param(self, services, param_name, optional=False):

        try:
            value = services.get_config_param(param_name)
            print(param_name, ' = ', value)
        except Exception :
            if optional: 
                print('config parameter ', param_name, ' not found')
                value = None
            else:
                message = 'required config parameter ', param_name, ' not found'
                print(message)
                services.exception(message)
                raise
        
        return value

    # Try to get component specific config parameter - wraps the exception handling
    def get_component_param(self, services, param_name, optional=False):

        if hasattr(self, param_name):
            value = getattr(self, param_name)
            print(param_name, ' = ', value)
        elif optional:
            print('optional config parameter ', param_name, ' not found')
            value = None
        else:
            message = 'required component config parameter ', param_name, ' not found'
            print(message)
            services.exception(message)
            raise
        
        return value
