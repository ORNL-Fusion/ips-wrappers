#! /usr/bin/env python

# Version 9 (Batchelor 6/2/2011)
# Moved LH component into the list of components that run concurrently


# Version 8 (Batchelor 1/11/2011)
# Adds generic awareness of RF_EC component

# Version 7 (Batchelor 11/19/2010)

# Now writes a file in SIM_ROOT called "PORTAL_RUNID" that contains the 
# PORTAL_RUNID number of the present run.  On restart it appends the new runid
# to the file so there is a record of all portal run id numbers associated with
# this simulation.

# Version 6 (Batchelor 6/30/2010)

# In pre-step_logic eliminated the copy of next_state_file to cur_state_file.  This was
# clobbering the plasma states accumulated from the services.merge_current_plasma_state
# calls.  Besides which TSC was not writing power_ic and power_nbi for the next time step
# into it.  Which was the whole purpose of next_plasma_state.  So now the RF and NBI 
# components don't get the news that the power has changed until the time step is over.
# We need to fix that but with a different method.  However it is possible that some 
# components might still use next_state_file so:
# Put in a test to see if there is a NEXT_STATE before copying between it and
# CURRENT_STATE.  This way it won't crash if the simulation doesn't use NEXT_STATE


# version 5 (Batchelor 5/21/2010)

# This version also has the additions needed to allow checkpoint/restart.  The 
# initialization section of the STEP function below now checks to see if config
# parameter SIMULATION_MODE = 'RESTART'.  If so the restart function is called for
# all components other than the INIT comoponent and the DRIVER component itself.
# At the end of the simulation, services.checkpoint_components is called which
# in turn calls the checkpoint functions for all the components except the DRIVER,
# but including the INIT component.  Since the framework can't call a DRIVER function
# from the DRIVER, the DRIVER's checkpoint function is called internally by the DRIVER.
# This is almost a moot point because as of now the DRIVER checkpoint function doesn't
# do anything, although later the DRIVER certainly could have state that needs to be saved.

# version 4 (Batchelor 2/4/2010)

# Version 4 from 2/4/2010 is more generic than the previous versions in that it picks
# up the components used from the PORTS global variable in the simulation config file.
# Only components listed in the PORTS variables are instantiated, initialized, and
# called during the time loop steps.  Therefore it is no longer necessary to have a
# different driver for each combination of components to be used.  In this version
# the components looked for are: EPA, RF_IC, RF_LH, NB, FP, FUS, and MONITOR. It also
# supports exception handling. When other components are available we will add them. 
#
# The time loop for this version is simple and explicit.  It steps first through all the
# source components (RF, NB, FUS, FP), then steps the EPA and finally the MONITOR.  As 
# with any driver the user can easily customize for any time stepping algorithm he wants.
#
# Note that this version launches the components using the services.call() function.  The
# framework normally launches jobs using MPI = APRUN, which requires an MPP allocation.  
# So either the job must be submitted via a batch script with an MPP repo specification, 
# or one must start an interactive session. It can't be launched from the command line 
# unless the machine configuration file (e.g. franklin.conf) is modified to contain 
# "MPIRUN=eval" instead of "MPIRUN=aprun". The intention is to write a companion to this 
# driver that uses the Python subprocess.call() to launch components.  Then the IPS can be
# started from the command line but the job will be restricted to one processor.
#
# This version also has some of the additions needed to allow checkpoint/restart but
# that is under continued development.

import sys
import os
import subprocess
import getopt
import shutil
import math
from component import Component
from Scientific.IO.NetCDF import *
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

      # Get plasma state files into driver work directory and copy to ps next if there is one
        services.stage_plasma_state()
        cur_state_file = services.get_config_param('CURRENT_STATE')
        try:
            next_state_file = services.get_config_param('NEXT_STATE')
            shutil.copyfile(cur_state_file, next_state_file)
        except Exception, e:
            print 'concurrent_driver: No NEXT_STATE file ', e        
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

        print ' init sequence complete--ready for time loop'

# ------------------------------------------------------------------------------
#
# Start Physics Layer
#
# ------------------------------------------------------------------------------


        # Iterate through the timeloop
        for t in tlist_str[1:len(timeloop)]:
            print (' ')
            print 'Driver: step to time = ', t
            services.update_time_stamp(t)

        # call pre_step_logic
            services.stage_plasma_state()
            self.pre_step_logic(float(t))
            services.update_plasma_state()
            print (' ')

       # Call step for each component

            print (' ')
            
            call_list = []

            if 'EPA' in port_names:
                call_id_epa  = services.call_nonblocking(epaComp, 'step', t)
                call_list.append(call_id_epa)

            if 'RF_EC' in port_names:
                call_id_rf  = services.call_nonblocking(rf_ecComp, 'step', t)
                call_list.append(call_id_rf)

            if 'RF_IC' in port_names:
                call_id_rf  = services.call_nonblocking(rf_icComp, 'step', t)
                call_list.append(call_id_rf)

            if 'RF_LH' in port_names:
                self.component_call_nonblocking(services, 'RF_LH', rf_lhComp, 'step', t)

            if 'NB' in port_names:
                call_id_nb  = services.call_nonblocking(nbComp, 'step', t)
                call_list.append(call_id_nb)
                
            retval_dict = services.wait_call_list(call_list, block=True)
            print 'retval_dict ',retval_dict
            print 'present working directory ', os.getcwd()
            print 'should be file name of partial rf_ic state ', retval_dict[call_id_rf]
#            print 'should be file name of partial nb state ', retval_dict[call_id_nb]

            if 'FUS' in port_names:
                self.component_call(services, 'FUS', fusComp, 'step', t)

            if 'FP' in port_names:
                self.component_call(services, 'FP', fpComp, 'step', t)

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

        if 'FP' in port_names:
            self.component_call(services, 'FP', fpComp, 'finalize', t)

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
                services.exception(comp_mode_string + ' failed')
                raise
            else:
                print comp_mode_string + ' finished'
                print (' ')
            
            return
    
    
    # Pre Step Logic
    def pre_step_logic(self, timeStamp):

        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        prior_state_file = self.services.get_config_param('PRIOR_STATE')

      #  Copy data from next plasma state (if there is one) to current plasma state
        try:
            next_state_file = services.get_config_param('NEXT_STATE')
            shutil.copyfile(next_state_file, cur_state_file)
        except Exception, e:
            pass

      # Update time stamps
        ps = NetCDFFile(cur_state_file, 'r+')
        t1 = ps.variables['t1'].getValue()
        self.services.log('ps%t1 = ', t1)

        power_ic = 0.0
        if ('power_ic' in ps.variables.keys()):
            power_ic = ps.variables['power_ic'].getValue()[0]
            print'concurrent_driver pre_step_logic: power_ic = ', power_ic

        ps.variables['t0'].assignValue(t1)
        ps.variables['t1'].assignValue(timeStamp)

        ps.close()
        shutil.copyfile(cur_state_file, prior_state_file)

        print'concurrent_driver pre_step_logic: timeStamp = ', timeStamp

        return


