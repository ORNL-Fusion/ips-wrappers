#! /usr/bin/env python
"""
driver_torlh_iterate_pwrscale.py  Version 0.0 (Batchelor 9/12/2017)
Driver for torlh/cql3d iteration with logic exposed here rather than buried in the 
components.  Adapted from generic_driver.py. Added coding from previous MIT script to
iterate CQL3D variable pwrscale until CQL3D power matches desired input power.
"""


import sys
import os
import subprocess
import getopt
import shutil
import math
from component import Component
from netCDF4 import *
import Numeric


class generic_driver(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

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
        services.stage_plasma_state()
        services.stage_input_files(self.INPUT_FILES)

      # get list of ports
        ports = services.getGlobalConfigParameter('PORTS')
        port_names = ports['NAMES'].split()
        print 'PORTS =', port_names
        port_dict = {}
        port_id_list = []

      # Instantiate components in port_names list, except DRIVER itself
        print (' ')

        # INIT is already instantiated by the framework. This adds it to port_dict
        if 'INIT' in port_names:
            initComp = services.get_port('INIT') 
            if(initComp == None):
                print 'Error accessing INIT component'
                raise
            port_dict['INIT'] = initComp
            port_id_list.append(initComp)
            print (' ')

        if 'EPA' in port_names:
            epaComp = services.get_port('EPA')
            if(epaComp == None):
                print 'Error accessing EPA component'
                raise
            port_dict['EPA'] = epaComp
            port_id_list.append(epaComp)
            print (' ')
       
        if 'RF_EC' in port_names:
            rf_ecComp = services.get_port('RF_EC')
            if(rf_ecComp == None):
                print 'Error accessing RF_EC component'
                raise
            port_dict['RF_EC'] = rf_ecComp
            port_id_list.append(rf_ecComp)
            print (' ')
        
        if 'RF_IC' in port_names:
            rf_icComp = services.get_port('RF_IC')
            if(rf_icComp == None):
                print 'Error accessing RF_IC component'
                raise
            port_dict['RF_IC'] = rf_icComp
            port_id_list.append(rf_icComp)
            print (' ')
        
        if 'RF_LH' in port_names:
            rf_lhComp = services.get_port('RF_LH')
            if(rf_lhComp == None):
                print 'Error accessing RF_LH component'
                raise
            port_dict['RF_LH'] = rf_lhComp
            port_id_list.append(rf_lhComp)
            print (' ')

        if 'NB' in port_names:
            nbComp = services.get_port('NB')
            if(nbComp == None):
                print 'Error accessing NB component'
                raise
            port_dict['NB'] = nbComp
            port_id_list.append(nbComp)
            print (' ')

        if 'FUS' in port_names:
            fusComp = services.get_port('FUS')
            if(fusComp == None):
                print 'Error accessing FUS component'
                raise
            port_dict['FUS'] = fusComp
            port_id_list.append(fusComp)
            print (' ')
 
        if 'FP' in port_names:
            fpComp = services.get_port('FP')
            if(fpComp == None):
                print 'Error accessing FP component'
                raise
            port_dict['FP'] = fpComp
            port_id_list.append(fpComp)
            print (' ')

        if 'MONITOR' in port_names:
            monitorComp = services.get_port('MONITOR')
            if(monitorComp == None):
                print 'Error accessing MONITOR component'
                raise       
            port_dict['MONITOR'] = monitorComp
            port_id_list.append(monitorComp)
            print (' ')

      # Is this a simulation startup or restart
        sim_mode = services.getGlobalConfigParameter('SIMULATION_MODE')
        print 'SIMULATION_MODE =', sim_mode

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
        services.stage_plasma_state()
        cur_state_file = services.get_config_param('CURRENT_STATE')
        try:
            next_state_file = services.get_config_param('NEXT_STATE')
            shutil.copyfile(cur_state_file, next_state_file)
        except Exception, e:
            print 'generic_driver: No NEXT_STATE file ', e        
        services.update_plasma_state()

       # Get Portal RUNID and save to a file
        run_id = services.get_config_param('PORTAL_RUNID')
        sym_root = services.getGlobalConfigParameter('SIM_ROOT')
        path = os.path.join(sym_root, 'PORTAL_RUNID')
        runid_file = open(path, 'a')
        runid_file.writelines(run_id + '\n')
        runid_file.close()

       # Post init processing: stage plasma state, stage output
        services.stage_output_files(t, self.OUTPUT_FILES)

        print ' init sequence complete'

        INIT_ONLY = self.get_component_param(services, 'INIT_ONLY', optional = True)
        if INIT_ONLY in [True, 'true', 'True', 'TRUE']:   
            message = 'INIT_ONLY: Intentional stop after INIT phase'
            print message
            return

# ------------------------------------------------------------------------------
#
# Start Physics Layer
#
# ------------------------------------------------------------------------------


        print ' \nZeroth step - Run torlh only in toric mode and qldce mode for Maxwellian'
        if sim_mode == 'NORMAL' :   # i.e. not RESTART do Maxwellian runs
            t = tlist_str[0]
            print (' ')
            print '\nDriver: starting iteration ', t
            services.update_time_stamp(t)

            try:
                services.call(rf_lhComp, 'step', t, toric_Mode = 'toric', inumin_Mode = 'Maxwell' , isol_Mode = '1')
            except Exception:
                message = 'RF_LH toric mode step failed'
                print message
                services.exception(message)
                raise 

            try:
                services.call(rf_lhComp, 'step', t, toric_Mode = 'qldce', \
                inumin_Mode = 'Maxwell' , isol_Mode = '1')
            except Exception:
                message = 'RF_LH qldce mode step failed'
                print message
                services.exception(message)
                raise 

        # Iterate through the timeloop, or in this case iteration loop
        for t in tlist_str[1:len(timeloop)]:
            print (' ')
            print '\nDriver: starting iteration ', t
            services.update_time_stamp(t)

        # call pre_step_logic
            services.stage_plasma_state()
            self.pre_step_logic(float(t))
            services.update_plasma_state()
            print (' ')

       # Call step for each component

            print (' ')

            if 'EPA' in port_names:
                self.component_call(services, 'EPA', epaComp, 'step', t)

        #*********************************************************************************
        # Iterate CQL3D pwrscale to match power_lh
        #*********************************************************************************

            icount=0
            running= True
            pwrscale=1.0
            tot_pwr = 0.
            hist_pwrscale=[]
            hist_pwr_result=[]
            
            # Get power_lh <--> goal_pwr from plasma state
            ps = Dataset(cur_state_file, 'r', format = 'NETCDF3_CLASSIC')
            goal_pwr = ps.variables['power_lh'][0]
            while running :
                icount=icount+1
                hist_pwrscale.append(pwrscale)

                comment =  'pwrscale iteration, icount = ' + str(icount) + ' pwrscale = '\
                            + str(pwrscale) + ' goal_pwr = ' + str(goal_pwr) + ' tot_pwr = '\
                            + str(tot_pwr)
                print comment
                services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                  event_comment =  comment)
               
                # Run CQL3D
                try:
                    services.call(fpComp, 'step', t, icount_arg = icount, pwrscale_arg = pwrscale)
                except Exception:
                    message = 'FP step failed'
                    print message
                    services.exception(message)
                    raise 

                tot_pwr = self.services.get_config_param('Pe_LH')
                hist_pwr_result.append(tot_pwr)
                print "history of pwrscale"
                print hist_pwrscale
                print "history of power results"
                print hist_pwr_result
             
                 # Iteration logic
                if ((tot_pwr<(1.05*goal_pwr)) and (tot_pwr>(0.95*goal_pwr))) or icount>10:
                     running= False           
                else:
                     if icount==1:
                         new_pwrscale=pwrscale*goal_pwr/tot_pwr
                     else:
                         new_pwrscale=pwrscale+(hist_pwrscale[icount-1]-hist_pwrscale[icount-2]) \
                                       /(hist_pwr_result[icount-1]-hist_pwr_result[icount-2])*(goal_pwr-tot_pwr)
                         if new_pwrscale<0:
                            new_pwrscale=pwrscale*goal_pwr/tot_pwr 
                     pwrscale=new_pwrscale
                     
            # After convergence record pwrscale history
            f=open('pwrscale.dat','w')
            str_icount=str(icount)
            f.write("\nnumber of feedback\n") 
            f.write(str_icount+"\n")
            f.write("pwrscale\n")
            for j in range(icount):
                str_pwrscale=str(hist_pwrscale[j])+"  "
                f.write(str_pwrscale)
            f.write("\n power in cql3D\n")
            for j in range(icount):
                str_pwr_result=str(hist_pwr_result[j])+"  "
                f.write(str_pwr_result)


        #*********************************************************************************
        # End pwrscale iteration
        #*********************************************************************************

            try:
                services.call(rf_lhComp, 'step', t, toric_Mode = 'toric', \
                inumin_Mode = 'nonMaxwell' , isol_Mode = '1')
            except Exception:
                message = 'RF_LH toric mode step failed'
                print message
                services.exception(message)
                raise 

            try:
                services.call(rf_lhComp, 'step', t, toric_Mode = 'qldce', \
                inumin_Mode = 'nonMaxwell' , isol_Mode = '1')
            except Exception:
                message = 'RF_LH qldce mode step failed'
                print message
                services.exception(message)
                raise 

            if 'MONITOR' in port_names:
                self.component_call(services, 'MONITOR', monitorComp, 'step', t)

            services.stage_plasma_state()

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
        print 'generic_driver.checkpoint() called'
        

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
                print message
                services.exception(message)
                raise 
            
            return
    
    
    # Pre Step Logic
    def pre_step_logic(self, timeStamp):

        cur_state_file = self.services.get_config_param('CURRENT_STATE')

      #  Copy data from next plasma state (if there is one) to current plasma state
        try:
            next_state_file = services.get_config_param('NEXT_STATE')
            shutil.copyfile(next_state_file, cur_state_file)
        except Exception, e:
            pass

      # Update time stamps

        ps = Dataset(cur_state_file, 'r+', format = 'NETCDF3_CLASSIC')
        t1 = ps.variables['t1'].getValue()
        ps.variables['t0'].assignValue(t1)
        ps.variables['t1'].assignValue(timeStamp)
        self.services.log('ps%t1 = ', t1)

        power_ic = 0.0
        if ('power_ic' in ps.variables.keys()):
            power_ic = ps.variables['power_ic'][:]
            print'generic_driver pre_step_logic: power_ic = ', power_ic

        ps.close()
        
    # Copy current plasma state to prior state if there is one
        try:
            prior_state_file = self.services.get_config_param('PRIOR_STATE')
            shutil.copyfile(cur_state_file, prior_state_file)
        except Exception, e:
            pass       
        
        print'generic_driver pre_step_logic: timeStamp = ', timeStamp
        
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
            print param_name, ' = ', value
        except Exception :
            if optional: 
                print 'config parameter ', param_name, ' not found'
                value = None
            else:
                message = 'required config parameter ', param_name, ' not found'
                print message
                services.exception(message)
                raise
        
        return value

    # Try to get component specific config parameter - wraps the exception handling
    def get_component_param(self, services, param_name, optional=False):

        if hasattr(self, param_name):
            value = getattr(self, param_name)
            print param_name, ' = ', value
        elif optional:
            print 'optional config parameter ', param_name, ' not found'
            value = None
        else:
            message = 'required component config parameter ', param_name, ' not found'
            print message
            services.exception(message)
            raise
        
        return value
