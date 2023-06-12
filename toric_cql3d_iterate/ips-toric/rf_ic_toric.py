#! /usr/bin/env python

"""
SF Cleanup 05/2023: RF Component for running TORIC in IPS
"""
from __future__ import print_function
import sys
import os
import subprocess
import getopt
import shutil
import string
from  component import Component
#Numeric should be replace by numpy, if needed -JCW
from Numeric import *
from Scientific.IO.NetCDF import *
import toric_tools as tt

INIT_Complete = False
CQL_COUPLE_MODE = False

class torlh (Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    # Init Function
    #---------------------------------------------------------------------------------
    def init(self, timeStamp=0):
        print('\ntorlh.init() called')

        services = self.services
        workdir = services.get_working_dir()

        # Get global configuration parameters
        cur_state_file = self.try_get_config_param(services,'CURRENT_STATE')
        cur_eqdsk_file = self.try_get_config_param(services,'CURRENT_EQDSK')
        cur_cql_file = self.try_get_config_param(services,'CURRENT_CQL')
        cur_dql_file = self.try_get_config_param(services,'CURRENT_DQL')

        # Get component-specific configuration parameters. Note: Not all of these are
        # used in 'init' but if any are missing we get an exception now instead of
        # later
        BIN_PATH = self.try_get_component_param(services,'BIN_PATH')
        RESTART_FILES = self.try_get_component_param(services,'RESTART_FILES')
        NPROC = self.try_get_component_param(services,'NPROC')
        global CQL_COUPLE_MODE
        CQL_COUPLE_MODE = self.try_get_component_param(services,'CQL_COUPLE_MODE', \
                                optional = True)

        torlh_log = os.path.join(workdir, 'log.torlh')

        # Copy plasma state files over to working directory
        try:
            services.stage_state()
        except Exception:
            logMsg = 'Error in call to stage_plasma_state()'
            self.services.exception(logMsg)
            raise 

        # Get input files
        try:
            services.stage_input_files(self.INPUT_FILES)
        except:
            logMsg = 'Error in call to stageInputFiles()'
            self.services.exception(logMsg)
            raise

        # Copy machine.inp_<suffix> to generic file name -> machine.inp if there is
        # a suffix
        try:
            suffix = self.INPUT_SUFFIX
            have_suffix = True
            # If suffix is not empty put an underscore in front of it.
            if len(suffix) > 0:
                print('INPUT_SUFFIX = ', suffix)
                suffix = '_' + suffix
            # If suffix is empty you don't really have one
            else:
                have_suffix = False
        except:
            have_suffix = False
            pass

        # If there is a non-empty suffix, copy to generic filename
        if have_suffix:
            try:
                shutil.copyfile('machine.inp' + suffix, 'machine.inp')
            except IOError as xxx_todo_changeme1:
                (errno, strerror) = xxx_todo_changeme1.args
                print('Error copying file %s to %s' % ('machine.inp' + suffix, 'machine.inp', strerror))
                logMsg = 'Error copying machine.inp_<suffix> -> machine.inp'
                services.exception(logMsg)
                raise

        do_init = os.path.join(self.BIN_PATH, 'do_torlh_init_abr')
        retcode = subprocess.call([do_init, cur_state_file])
        if (retcode != 0):
            logMsg = 'Error in call to torlh_init'
            self.services.error(logMsg)
            raise Exception(logMsg)
            
        # Update plasma state files in plasma_state work directory
        try:
            services.update_state()
        except Exception:
            logMsg = 'Error in call to update_plasma_state()'
            self.services.exception(logMsg)
            raise 

        # Archive output files
        # N.B.  do_torlh_init does not produce a complete set of torlh output
        #       files.  This causes an error in stage_output_files().  To
        #       solve this we generate a dummy set of output files here with
        #       system call 'touch'
        for file in self.OUTPUT_FILES.split():
            print('touching ', file)
            subprocess.call(['touch', file])
        # Now stage them
        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception:
            logMsg = 'Error in call to stage_output_files()'
            self.services.exception(logMsg)
            raise 

        if CQL_COUPLE_MODE in [True, 'true', 'True', 'TRUE']:
            self.step(timeStamp)
            
        global INIT_Complete
        INIT_Complete = True

        return 0

    # Restart Function
    #---------------------------------------------------------------------------------
    def restart(self, timeStamp):
        print('\ntorlh.restart() called')

        services = self.services
        workdir = services.get_working_dir()

        # Get restart files listed in config file.
        try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
        except Exception:
            logMsg = 'Error in call to get_restart_files()'
            self.services.exception(logMsg)
            raise

        # Get global configuration parameters
        try:
            self.plasma_state_file = services.get_config_param('CURRENT_STATE')
            self.eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            self.torlh_log = os.path.join(workdir, 'log.torlh')
        except:
            logMsg = 'torlh restart: error in getting config parameters'
            self.services.exception(logMsg)
            raise 
            
        return 0

    # Step Function
    #---------------------------------------------------------------------------------
    def step(self, timeStamp):
        """Take a step for the torlh component.  Really a complete run."""
        print('\ntorlh.step() called')

        if (self.services == None):
            logMsg = 'Error in torlh: step (): No self.services'
            self.services.error(logMsg)
            raise Exception(logMsg)
        services = self.services
        workdir = services.get_working_dir()

      # Copy plasma state files over to working directory
        try:
            services.stage_state()
        except:
            logMsg = 'Error in call to stage_plasma_state()'
            self.services.exception(logMsg)
            raise 

      # Get input files
        try:
            services.stage_input_files(self.INPUT_FILES)
        except:
            logMsg = 'Error in call to stageInputFiles()'
            self.services.exception(logMsg)
            raise 

      # Copy machine.inp_<suffix> to generic file name -> machine.inp if there is
      # a suffix
        try:
            suffix = self.INPUT_SUFFIX
            have_suffix = True
        # If suffix is not empty put an underscore in front of it.
            if len(suffix) > 0:
                print('INPUT_SUFFIX = ', suffix)
                suffix = '_' + suffix
        # If suffix is empty you don't really have one
            else:
                have_suffix = False
        except:
            have_suffix = False
            pass

        # If there is a non-empty suffix, copy to generic filename
        if have_suffix:
            try:
                shutil.copyfile('machine.inp' + suffix, 'machine.inp')
            except IOError as xxx_todo_changeme2:
                (errno, strerror) = xxx_todo_changeme2.args
                print('Error copying file %s to %s' % ('machine.inp' + suffix,
                'machine.inp', strerror))
                logMsg = 'Error copying machine.inp_<suffix> -> machine.inp'
                services.exception(logMsg)
                raise 

        prepare_input = os.path.join(self.BIN_PATH, 'prepare_torlh_input_abr')
        process_output  = os.path.join(self.BIN_PATH, 'process_torlh_output_mcmd')

        zero_RF_LH_power = self.ZERO_LH_POWER_BIN
        torlh_bin = self.TORLH_BIN
        prepare_eqdsk  = self.GEQXPL_BIN

    # Get global configuration parameters
        cur_state_file = self.try_get_config_param(services,'CURRENT_STATE')
        cur_eqdsk_file = self.try_get_config_param(services,'CURRENT_EQDSK')
        cur_cql_file = self.try_get_config_param(services,'CURRENT_CQL')
        cur_dql_file = self.try_get_config_param(services,'CURRENT_DQL')
        # enorm which is used here and in cql3d
        arg_enorm = 'None'
        arg_enorm = self.try_get_config_param(services,'ENORM', optional = True)

        torlh_log = os.path.join(workdir, 'log.torlh')
        cwd = os.getcwd()

# Check if LH power is zero (or effectively zero).  If true don't run torlh just
# run zero_RF_LH_power fortran code
        print('cur_state_file = ', cur_state_file)
        ps = NetCDFFile(cur_state_file, 'r')
        power_lh = ps.variables['power_lh'].getValue()[0]
        ps.close()
        print('power = ', power_lh)
        if(-0.02 < power_lh < 0.02):
            print(zero_RF_LH_power)
            services.send_portal_event(event_type = 'COMPONENT_EVENT',\
              event_comment =  'running ' + zero_RF_LH_power)
            retcode = subprocess.call([zero_RF_LH_power, cur_state_file])
            if (retcode != 0):
                logMsg = 'Error executing ' + prepare_input
                self.services.error(logMsg)
                raise Exception(logMsg)

            # N.B. zero_RF_LH_power does not produce a complete set of torlh output
            #      files.  This causes an error in stage_output_files().  To
            #      solve this we generate a dummy set of output files here with
            #      system call 'touch'
            for file in self.OUTPUT_FILES.split():
                subprocess.call(['touch', file])

# Check if LH power is negative.  If true don't run torlh just
# retain power from previous time step i.e. leave sources untouched in the state.
# However power_lh needs to be reset back to positive

        elif( power_lh < -0.02):
            print('continuing power from previous time step')
            ps.variables['power_lh'].assignValue(-power_lh)
            ps.close()

    # Or actually run torlh

        else:

            if not os.path.isfile(prepare_input):
                logMsg = 'Cannot find torlh prepare_input binary: ' + prepare_input
                self.services.error(logMsg)
                raise Exception(logMsg)

# ------------------------------------------------------------------------------                
            global run_ImChizz
            if CQL_COUPLE_MODE in [True, 'true', 'True', 'TRUE'] and INIT_Complete == True: # Will be False during INIT
                try:
                    print('\nRunning ImChizz')
                    subprocess.call(['cp', cur_cql_file, 'cql3d.cdf' ])
                except Exception:
                    message = 'generic_ps_init: Error copying CURRENT_CQL_FILE to cql3d.cdf'
                    print(message)
                    services.exception(message)
                    raise              
                imchzz_bin = self.ImChizz_BIN
                cmd_imchizz=self.ImChizz_BIN
                try:
                   services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                      event_comment =  'running ' + cmd_imchizz)
                   P=subprocess.Popen(cmd_imchizz,stdin=subprocess.PIPE,stdout=subprocess.PIPE,\
                      stderr=subprocess.STDOUT, bufsize=1)
                #      stderr=subprocess.STDOUT, bufsize=1)
                except :
                   logMsg = "Error executing" + cmd_imchizz
                   self.services.error(logMsg)
                   raise
               # P.stdin.write("b\n")
                print(P.communicate("b\n"))
                #P.wait()
                print('Finished ImChizz')

# ------------------------------------------------------------------------------                
        # Run in toricmode = 'toric'
            # Call toric prepare_input to generate torica.inp

            arg_toric_Mode = 'toric'
            arg_isol_Mode = '1'           
            arg_inumin_Mode = 'Maxwell'
            if INIT_Complete and CQL_COUPLE_MODE:
                arg_inumin_Mode = 'nonMaxwell'                
            
            cmd_prepare_input = [prepare_input, cur_state_file, arg_toric_Mode,\
                      arg_inumin_Mode,arg_isol_Mode, arg_enorm]
                      
            print('running = ', cmd_prepare_input)
            services.send_portal_event(event_type = 'COMPONENT_EVENT',\
              event_comment =  cmd_prepare_input)
            retcode = subprocess.call(cmd_prepare_input)
            if (retcode != 0):
                logMsg = 'Error executing ' + prepare_input
                self.services.error(logMsg)
                raise Exception(logMsg)

            # Call xeqdsk_setup to generate eqdsk.out file
            cmd_eqdsk = [prepare_eqdsk, '@equigs_gen', '/g_filename='+cur_eqdsk_file,\
                                       '/equigs_filename=equigs.data']
            print('running', cmd_eqdsk)
            services.send_portal_event(event_type = 'COMPONENT_EVENT',\
              event_comment =  cmd_eqdsk)
            retcode = subprocess.call(cmd_eqdsk)
            if (retcode != 0):
                logMsg = 'Error in call to prepare_eqdsk'
                self.services.error(logMsg)
                raise Exception(logMsg)

            # Launch toric executable
            print('torlh processors = ', self.NPROC)
            cwd = services.get_working_dir()
            task_id = services.launch_task(self.NPROC, cwd, torlh_bin, logfile=torlh_log)
            retcode = services.wait_task(task_id)
            if (retcode != 0):
                logMsg = 'Error executing command: ' + torlh_bin
                self.services.error(logMsg)
                raise Exception(logMsg)

            
            
            # Preserve torica.out from run in toric mode
            try:
                shutil.copyfile('torica.out', 'torica_toricMode.out')
            except IOError as xxx_todo_changeme3:
                (errno, strerror) = xxx_todo_changeme3.args
                logMsg =  'Error copying file %s to %s' % ('torica.out', 'torica_toricMode.out'\
                        , strerror)
                print(logMsg)
                services.exception(logMsg)
                raise 
            
                
# ------------------------------------------------------------------------------                
        # Run in toricmode = 'qldci' if needed
            if CQL_COUPLE_MODE in [True, 'true', 'True', 'TRUE']:

                arg_toric_Mode = 'qldci'
                arg_isol_Mode = '1'            
                arg_inumin_Mode = 'Maxwell'
                if INIT_Complete:
                    arg_inumin_Mode = 'nonMaxwell'

                cmd_prepare_input = [prepare_input, cur_state_file, arg_toric_Mode,\
                      arg_inumin_Mode,arg_isol_Mode, arg_enorm]
                print('running = ', cmd_prepare_input)
                services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                  event_comment =  cmd_prepare_input)
                retcode = subprocess.call(cmd_prepare_input)
                if (retcode != 0):
                    logMsg = 'Error executing ' + prepare_input
                    self.services.error(logMsg)
                    raise Exception(logMsg)
# 
#                 print '\nRunning torlh in qldce mode, inumin_Mode = ', arg_inumin_Mode
#                 retcode = subprocess.call([prepare_input, cur_state_file, arg_toric_Mode,\
#                       arg_inumin_Mode,arg_isol_Mode])
#                 if (retcode != 0):
#                     logMsg = 'Error executing ' + prepare_input
#                     self.services.error(logMsg)
#                     raise Exception(logMsg)

                # Launch torlh executable
                print('torlh processors = ', self.NPROC)
                cwd = services.get_working_dir()
                task_id = services.launch_task(self.NPROC, cwd, torlh_bin, logfile=torlh_log)
                retcode = services.wait_task(task_id)
                if (retcode != 0):
                    logMsg = 'Error executing command: ' + torlh_bin
                    self.services.error(logMsg)
                    raise Exception(logMsg)

            # Preserve torica.out from run in qldce mode
                try:
                    shutil.copyfile('torica.out', 'torica_qldceMode.out')
                except IOError as xxx_todo_changeme:
                    (errno, strerror) = xxx_todo_changeme.args
                    logMsg =  'Error copying file %s to %s' % ('torica.out', 'torica_qldceMode.out'\
                            , strerror)
                    print(logMsg)
                    services.exception(logMsg)
                    raise
                     
            # Run mapin
                mapin_bin = self.try_get_component_param(services,'MAPIN_BIN')
                print('\nRunning ' + mapin_bin)
                services.send_portal_event(event_type = 'COMPONENT_EVENT', \
                     event_comment = 'running ' + mapin_bin)
                retcode = subprocess.call([mapin_bin])
                if (retcode != 0):
                    logMsg = 'Error executing ' + mapin_bin
                    self.services.error(logMsg)
                    raise Exception(logMsg)
                    
            # Call process_output
            # First rename default fort.* to expected names by component method as of torlh5 r918 from ipp
            #os.rename('fort.9','torlh_cfg.nc')
            #os.rename('fort.21','torlh.nc')
            # No process_output code yet.  And don't find fort.9 or fort.21 in work directory.
            # retcode = subprocess.call([process_output, cur_state_file])
#             if (retcode != 0):
#                 logMsg = 'Error executing' + process_output
#                 self.services.error(logMsg)
#                 raise Exception(logMsg)

# Merge partial plasma state containing updated IC data
        try:
            partial_file = cwd + '/RF_LH_' + cur_state_file
            # No process_output code yet
            #services.merge_current_plasma_state(partial_file, logfile='log.update_state')
            #print 'merged torlh plasma state data ', partial_file
            print('No process_output code yet, so no plasma state merge')
        except:
            logMsg = 'Error in call to merge_current_plasma_state(' + partial_file + ')'
            self.services.exception(logMsg)
            raise 

# Run IDL script if requested
#         do_idl_plots = self.try_get_component_param(services, 'DO_IDL_PLOTS', optional = True)
#         if do_idl_plots != None:
#             if do_idl_plots:
#                 print 'running run_IDL_toricplot()'
#                 self.run_IDL_toricplot()


      # Update plasma state files in plasma_state work directory
        try:
            services.update_state()
        except Exception:
            logMsg = 'Error in call to update_plasma_state()'
            self.services.exception(logMsg)
            raise 

      # Archive output files
        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except:
            logMsg = 'Error in call to stage_output_files()'
            self.services.exception(logMsg)
            raise 

        return 0

# ------------------------------------------------------------------------------
#
# checkpoint function
# Saves plasma state files to restart directory
#
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print('RF_LH_torlh.checkpoint() called')
        services = self.services
        services.save_restart_files(timestamp, self.RESTART_FILES)


# ------------------------------------------------------------------------------
#
# FINALIZE function
# As of now it does nothing
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print('torlh.finalize() called')


# ------------------------------------------------------------------------------
#
# "Private"  methods
#
# ------------------------------------------------------------------------------

    def try_get_config_param(self, services, param_name, optional=False):

        try:
            value = services.get_config_param(param_name)
            print(param_name, ' = ', value)
        except Exception:
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
    def try_get_component_param(self, services, param_name, optional=False):

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

    # IDL plots
    def run_IDL_toricplot(self):

         cmd_toricplot_pro=".r pltoriclhg.pro\n"
         P=subprocess.Popen(["idl"],stdin=subprocess.PIPE, stdout=subprocess.PIPE, \
         stderr=subprocess.STDOUT)
         P.stdin.write(cmd_toricplot_pro)
         P.stdin.write("1\n")
         P.stdin.write("0\n")
         P.stdin.write("1\n")
         print("Make a file(torica.ps)")
         
         return
