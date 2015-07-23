#! /usr/bin/env python
import sys
import os
import subprocess
import getopt
import shutil
import math
from component import Component
from Scientific.IO.NetCDF import *
import Numeric
from plasmastate import *
import traceback

class basic_time_step_driver(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0):
        print 80*'-'
        print 'Driver init'
        print self.__class__.__name__
        print 80*'-'
        return

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

        # instantiate components in port_names list, except DRIVER itself

        for port_name in port_names:
            if port_name in ["DRIVER"]: continue
            port = services.get_port(port_name) 
            if(port == None):
                logMsg = 'Error accessing '+port_name+' component'
                raise Exception(logMsg)
            port_dict[port_name] = port
            port_id_list.append(port)
 
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

        for port_name in port_names:
            if port_name in ['INIT','DRIVER']: continue 
            self.component_call(services,port_name,port_dict[port_name],init_mode,t)
 
        # Get plasma state files into driver work directory and copy to psn if there is one
        services.stage_plasma_state()
        cur_state_file = services.get_config_param('CURRENT_STATE')
        try:
            next_state_file = services.get_config_param('NEXT_STATE')
            shutil.copyfile(cur_state_file, next_state_file)
        except Exception:
            logMsg = 'generic_driver: No NEXT_STATE file '        
            services.exception(logMsg)
            raise
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

        ps = PlasmaState("ips",1)
        ps_work_dir = services.getGlobalConfigParameter('PLASMA_STATE_WORK_DIR')

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

            for port_name in port_names:
                if port_name in ['INIT','DRIVER']: continue 
                self.component_call(services,port_name,port_dict[port_name],'step',t)

                # Santiy check on plasma state density

                this_ps_file = ps_work_dir+'/'+cur_state_file
                ps.read(this_ps_file)
                nS = ps["ns"].shape[0]

                for nn in range(0,nS-1):

                    if np.isnan(np.sum(ps["ns"][nn])):
                        print 'Sanity checking density values for : '+ port_name
                        print os.getcwd()
                        print 'Reading this ps file : ' + this_ps_file
                        logMsg = 'ERROR : NaN detected after running : '+ port_name
                        services.exception(logMsg)
                        raise Exception(logMsg)
 
            services.stage_plasma_state()

            # Post step processing: stage plasma state, checkpoint components and self
            services.stage_output_files(t, self.OUTPUT_FILES)
            services.checkpoint_components(port_id_list, t)
            self.checkpoint(t)

        # Post simulation: call checkpoint on each component

        services.checkpoint_components(port_id_list, t, Force = True)
        self.checkpoint(t)
      
        # Post simulation: call finalize on each component
        print (' ')

        for port_name in port_names:
            if port_name in ['INIT','DRIVER']: continue 
            self.component_call(services, port_name, port_dict[port_name], 'finalize', t)


    def checkpoint(self, timestamp=0.0):
        print 'generic_driver.checkpoint() called'
        

    def finalize(self, timestamp = 0):
        logMsg = 'Finalizing driver - i.e., Simulation End :)'
        self.services.info(logMsg)
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
        except Exception:
            logMsg = 'INFO : No Next Plasma State'
            self.services.info(logMsg)
            pass

        # Update time stamps
        ps = NetCDFFile(cur_state_file, 'r+')
        t1 = ps.variables['t1'].getValue()
        self.services.log('ps%t1 = ', t1)

        power_ic = 0.0
        if ('power_ic' in ps.variables.keys()):
            power_ic = ps.variables['power_ic'].getValue()[0]
            print'generic_driver pre_step_logic: power_ic = ', power_ic

        ps.variables['t0'].assignValue(t1)
        ps.variables['t1'].assignValue(timeStamp)

        ps.close()
        shutil.copyfile(cur_state_file, prior_state_file)
        
        print'generic_driver pre_step_logic: timeStamp = ', timeStamp

        return


