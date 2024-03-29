#! /usr/bin/env python

from builtins import str
from builtins import range
import sys
import os
import subprocess
import getopt
import shutil
import string
from  component import Component
from netCDF4 import *
from simple_file_editing_functions import get_lines, lines_to_variable_dict

class cql3d(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        self.firstTime  = True
        self.curTime = -1.0
        self.prevTime = -1.0
        print('Created %s' % (self.__class__))

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp=0):
        print('cql3d.init() called')

        services = self.services

    # Get global configuration parameters
        try:
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            #cur_dql_file = services.get_config_param('CURRENT_DQL')
            cur_cql_file = services.get_config_param('CURRENT_CQL')
           
        except:
            print('fp_cql3d_general: error getting config parameters CURRENT_STATE CURRENT_EQDSK')
            services.error('fp_cql3d_general: error getting config parameters CURRENT_STATE CURRENT_EQDSK')
            raise Exception('rf_ic_gneray: error getting config parameters CURRENT_STATE CURRENT_EQDSK')


        #try:
        #     cur_ImChizz_inp_file = services.get_config_param('CURRENT_ImChizz_inp')
        #except:
        #    print('no ImChizz_inp file found')

        
    # Get component-specific configuration parameters. Note: Not all of these are
    # used in 'init' but if any are missing we get an exception now instead of
    # later
        try:
            NPROC = self.NPROC
            NPPN = self.NPPN
            BIN_PATH = self.BIN_PATH
            INPUT_FILES = self.INPUT_FILES
            OUTPUT_FILES = self.OUTPUT_FILES
            RESTART_FILES = self.RESTART_FILES
            BIN_PATH = self.BIN_PATH
            CQL3D_BIN = self.CQL3D_BIN
            CQL3D_MODE = self.CQL3D_MODE
            CQL3D_OUTPUT = self.CQL3D_OUTPUT
            CQL3D_NML = self.CQL3D_NML
            NSTEPS_STR = self.NSTEPS_STR
            DELTAT_STR = self.DELTAT_STR
            PS_ADD_NML = self.PS_ADD_NML
        except:
            print('fp_cql3d_general init: error getting cql3d-specific config parameters')
            services.error('fp_cql3d_general: error getting cql3d-specific\
            config parameters')
            raise Exception('fp_cql3d_general: error getting cql3d-specific\
            config parameters')

        # enorm which is used here and in cql3d
        arg_enorm = 'None'
        arg_enorm = self.try_get_config_param(services,'ENORM', optional = True)

        arg_eqover = 'n'
        arg_eqover = self.try_get_config_param(services,'EQ_OVERRIDE', optional=True)

        
    # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception as e:
          print('Error in call to stage_plasma_state()' , e)
          services.error('Error in call to stage_plasma_state()')
          raise Exception('Error in call to stage_plasma_state()')
        
    # Get input files  
        try:
          services.stage_input_files(INPUT_FILES)
        except Exception as e:
          print('Error in call to stageInputFiles()' , e)
          services.error('Error in call to stageInputFiles()')
          raise Exception('Error in call to stageInputFiles()')
            
    # Copy cqlinput_<suffix> to generic file name -> cqlinput if there is
    # a suffix
        try:
          suffix = self.INPUT_SUFFIX
          have_suffix = True
        # If suffix is not empty put an underscore in front of it.
          if len(suffix) > 0:
              print('cql3d INPUT_SUFFIX = ', suffix)      
              suffix = '_' + suffix
        # If suffix is empty you don't really have one
          else:
              have_suffix = False
        except:
          have_suffix = False
      
      # If there is a non-empty suffix, copy to generic filename 'cqlinput'
        if have_suffix:      
          try:
              shutil.copyfile('cqlinput' + suffix, 'cqlinput')
          except IOError as xxx_todo_changeme:
              (errno, strerror) = xxx_todo_changeme.args
              print('Error copying file %s to %s' % ('cqlinput' + suffix, 
              'cqlinput', strerror))
              services.error('Error copying cqlinput_<suffix> -> cqlinput')
              raise Exception('Error copying cqlinput_<suffix> -> cqlinput')

    # Copy current plasma state file to generic name -> cur_state.cdf
        try:
            shutil.copyfile(cur_state_file, 'cur_state.cdf')
        except IOError as xxx_todo_changeme3:
            (errno, strerror) = xxx_todo_changeme3.args
            print('Error copying file %s to %s' % (cur_state_file, 'cur_state.cdf', strerror))
            services.error('Error copying cur_state_file -> cur_state.cdf')
            raise Exception('Error copying cur_state_file -> cur_state.cdf')

    # Copy current eqdsk file to generic name -> eqdsk
        if(arg_eqover == 'n'):
            try:
                shutil.copyfile(cur_eqdsk_file, 'eqdsk')
            except IOError as xxx_todo_changeme4:
                (errno, strerror) = xxx_todo_changeme4.args
                print('Error copying file %s to %s' % (cur_eqdsk_file, 'eqdsk', strerror))
                services.error('Error copying cur_eqdsk_file -> eqdsk')
                raise Exception('Error copying cur_eqdsk_file -> eqdsk')

    # Copy current Dql file to generic name -> genray.nc
        #try:
        #    shutil.copyfile(cur_dql_file, 'genray.nc')
        #except IOError as xxx_todo_changeme5:
        #    (errno, strerror) = xxx_todo_changeme5.args
        #    print('Error copying file %s to %s' % (cur_dql_file, 'genray.nc', strerror))
        #    services.error('Error copying cur_dql_file -> genray.nc')
        #    raise Exception('Error copying cur_dql_file -> genray.nc')

        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_cql3d_input')

        cql3d_mode = self.CQL3D_MODE
        cql3d_output = self.CQL3D_OUTPUT
        cql3d_nml = self.CQL3D_NML
        nsteps_str = self.NSTEPS_STR
        deltat_str = self.DELTAT_STR
        ps_add_nml = self.PS_ADD_NML

    # Call prepare_input - init
        print('fp_cql3d: calling prepare_input init')        
        log_file = open('log_prepare_cql3d_input_init', 'w')
        ips_mode = 'init'

    # ptb: Set restart depending on whether the distribution function file from cql3d
    # already exist
        restart = 'disabled'
    #ptb&SS    if os.path.isfile('./cql3d.nc'): restart = 'enabled'
    #ptb&SS    if os.path.isfile('./cql3d.nc'): shutil.copyfile('cql3d.nc','distrfunc.nc')

        if os.path.isfile('./cql3d.nc'): 
            if os.stat('./cql3d.nc')[6] != 0:
                 restart = 'enabled'
                 shutil.copyfile('cql3d.nc','distrfunc.nc')

    # ptb: End of ptb hack

    # ptb:    command = prepare_input_bin + ' ' + ips_mode + ' ' + cql3d_mode  + ' ' +\
    # ptb:    cql3d_output + ' ' + cql3d_nml + ' ' + nsteps_str + ' ' + ps_add_nml
        command = prepare_input_bin + ' ' + ips_mode + ' ' + cql3d_mode  + ' ' +\
        cql3d_output + ' ' + cql3d_nml + ' ' + restart + ' '+ ps_add_nml
        if nsteps_str != None:
            command = command + ' ' + nsteps_str

        if deltat_str != None:
            command = command + ' ' + deltat_str
             
        if arg_enorm != None:
            command = command  + ' ' + arg_enorm
            
        print('running = ', command)
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
          event_comment =  command)
        
        retcode = subprocess.call(command.split(), stdout = log_file,\
                                  stderr = subprocess.STDOUT)
        if (retcode != 0):
            print('Error executing cql3d init ', prepare_input_bin)
            services.error('Error executing cql3d init')
            raise Exception('Error executing cql3d init')

    # Copy generic cur_state.cdf -> current plasma state file
        try:
            shutil.copyfile('cur_state.cdf', cur_state_file)
        except IOError as xxx_todo_changeme6:
            (errno, strerror) = xxx_todo_changeme6.args
            print('Error copying file %s to %s' % ('cur_state.cdf', cur_state_file, strerror))
            services.error('Error copying cur_state.cdf -> cur_state_file')
            raise Exception('Error copying cur_state.cdf -> cur_state_file')

    # Update plasma state files in plasma_state work directory
        try:
          services.update_state()
        except Exception as e:
          print('Error in call to update_plasma_state()', e)
          services.error('Error in call to update_plasma_state()')
          raise Exception('Error in call to update_plasma_state()')
     
    # Archive output files
    # N.B.  prepare_cql3d_input in init mode does not produce a complete set 
    #       of output files.  This causes an error in stage_output_files().
    #       To solve this we generate a dummy set of output files here with 
    #       system call 'touch'
        for file in self.OUTPUT_FILES.split():
          print('touching ', file)
          subprocess.call(['touch', file])
    # Now stage them      
        try:
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception as e:
          print('Error in call to stage_output_files()', e)
          services.error('Error in call to stage_output_files()')
          raise Exception('Error in call to stage_output_files()')

        return 0
# ------------------------------------------------------------------------------
#
# RESTART function
# Gets restart files from restart directory
# Loads the global configuration parameters from the config file
#
# ------------------------------------------------------------------------------
        
    def restart(self, timeStamp):
      print('cql3d.restart() called')

      services = self.services
      if (services == None) :
         print('Error in cql3d: step (): No services')
         services.error('Error in cql3d: step (): No services')
         raise Exception('Error in cql3d: step (): No services')

    # Get restart files listed in config file.        
      try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
      except Exception as e:
            print('Error in call to get_restart_files()' , e)
            services.error('cql3d: error in call to get_restart_files()')
            raise Exception('cql3d: error in call to get_restart_files()')

    # Get global configuration parameters
      try:
            self.plasma_state_file = services.get_config_param('CURRENT_STATE')
            self.eqdsk_file = services.get_config_param('CURRENT_EQDSK')
      except:
            print('cql3d restart: error in getting config parameters')
            services.error('cql3d restart: error in getting config parameters')
            raise Exception('cql3d restart: error in getting config parameters')

      return 0
        
# ------------------------------------------------------------------------------
#
# STEP function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp, **kwargs):
        print('cql3d.step() called')

        if (self.services == None) :
           print('Error in cql3d: step (): No services')
           raise Exception('Error in cql3d: step (): No services')
        services = self.services

        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_cql3d_input')
        process_output_bin  = os.path.join(self.BIN_PATH, 'process_cql3d_output')
        cql3d_bin = os.path.join(self.BIN_PATH, 'cql3d')

    # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception as e:
          print('Error in call to stage_plasma_state()' , e)
          services.error('Error in call to stage_plasma_state()')
          raise Exception('Error in call to stage_plasma_state()')

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        #cur_dql_file = services.get_config_param('CURRENT_DQL')
        cur_cql_file = services.get_config_param('CURRENT_CQL')
        #cur_ImChizz_inp_file = services.get_config_param('CURRENT_ImChizz_inp')
        print('CURRENT_CQL = ', cur_cql_file)
    # Copy current plasma state file to generic name -> cur_state.cdf
        try:
            shutil.copyfile(cur_state_file, 'cur_state.cdf')
        except IOError as xxx_todo_changeme7:
            (errno, strerror) = xxx_todo_changeme7.args
            print('Error copying file %s to %s' % (cur_state_file, 'cur_state.cdf', strerror))
            raise

    # Copy current eqdsk file to generic name -> eqdsk
        arg_eqover = 'n'
        arg_eqover = self.try_get_config_param(services,'EQ_OVERRIDE', optional=True)
        print('eqdsk override: ', arg_eqover)
        if (arg_eqover == 'n'):
            try:
                shutil.copyfile(cur_eqdsk_file, 'eqdsk')
            except IOError as xxx_todo_changeme8:
                (errno, strerror) = xxx_todo_changeme8.args
                print('Error copying file %s to %s' % (cur_eqdsk_file, 'eqdsk', strerror))
                services.error('Error copying cur_eqdsk_file -> eqdsk')
                raise Exception('Error copying cur_eqdsk_file -> eqdsk')

# Check if LHRF power is zero (or effectively zero).  If true don't run Genray

        ps = Dataset(cur_state_file, 'r+', format = 'NETCDF3_CLASSIC')
       
        ps.close()
        #print('power = ', power_lh)
        power_lh = 1.0 #SF fixed so I don't have to deal with this
        if(power_lh > 1.0E-04):

    # Copy current Dql file to generic name -> genray.nc
          #try:
          #    shutil.copyfile(cur_dql_file, 'genray.nc')
          #except IOError as xxx_todo_changeme1:
          #    (errno, strerror) = xxx_todo_changeme1.args
          #    print('Error copying file %s to %s' % (cur_dql_file, 'genray.nc', strerror))
          #    services.error('Error copying cur_dql_file -> genray.nc')
          #    raise Exception('Error copying cur_dql_file -> genray.nc')

          cql3d_mode = self.CQL3D_MODE
          cql3d_output = self.CQL3D_OUTPUT
          cql3d_nml = self.CQL3D_NML
          nsteps_str = self.NSTEPS_STR
          deltat_str = self.DELTAT_STR
          ps_add_nml = self.PS_ADD_NML

        # enorm which is used here and in cql3d
          arg_enorm = 'None'
          arg_enorm = self.try_get_config_param(services,'ENORM', optional = True)

        #if deltat_str .eq. 0 then use timestep ramp from GENRAY/CQL sims
          t0 = float(timeStamp)
          if(deltat_str=='0'):
              if(t0 < 5):
                  deltat_str = '0.00001'
              elif(t0<10):
                  deltat_str = '0.0001'
              elif(t0<30):
                  deltat_str = '0.001'
              elif(t0<40):
                  deltat_str = '0.005'
              else:
                  deltat_str = '0.01'

          
    # Call prepare_input - step
          print('fp_cql3d step: calling prepare_input')
          
          log_file = open('log_prepare_cql3d_input_step', 'w')
          ips_mode = 'step'

    # ptb: Set restart depending on whether the distribution function file from cql3d
    # already exist
          restart = 'disabled'
    #ptb&SS      if os.path.isfile('./cql3d.nc'): restart = 'enabled'
    #ptb&SS      if os.path.isfile('./cql3d.nc'): shutil.copyfile('cql3d.nc','distrfunc.nc')

          if os.path.isfile('./cql3d.nc'): 
              if os.stat('./cql3d.nc')[6] != 0:
                   restart = 'enabled'
                   shutil.copyfile('cql3d.nc','distrfunc.nc')
    # ptb: End of ptb hack

# If this is the first step in a pwrscale iteration, and restart = 'enabled' (e.g. this
# is a restart) save the distrfunc.nc file to initial_distrfunc.nc.   
          if 'icount_arg' in kwargs:
             icount = kwargs.get('icount_arg')
             if (icount == 1) and (restart == 'enabled'):
             	shutil.copyfile('distrfunc.nc', 'initial_distrfunc.nc')
             elif (icount > 1) and (restart == 'enabled'):
             	shutil.copyfile('initial_distrfunc.nc', 'distrfunc.nc')

    # ptb:    command = prepare_input_bin + ' ' + ips_mode + ' ' + cql3d_mode  + ' ' +\
    # ptb:    cql3d_output + ' ' + cql3d_nml + ' ' + nsteps_str + ' ' + ps_add_nml
          command = prepare_input_bin + ' ' + ips_mode + ' ' + cql3d_mode  + ' ' +\
          cql3d_output + ' ' + cql3d_nml + ' ' + restart + ' '+ ps_add_nml
          if nsteps_str != None:
            command = command + ' ' + nsteps_str

          if deltat_str != None:
            command = command + ' ' + deltat_str
             
          if arg_enorm != None:
            command = command  + ' ' + arg_enorm
            
          print('running', command)
          services.send_portal_event(event_type = 'COMPONENT_EVENT',\
              event_comment =  command)

          retcode = subprocess.call(command.split(), stdout = log_file,\
                                    stderr = subprocess.STDOUT)
          if (retcode != 0):
              print('Error executing cql3d: ', prepare_input_bin)
              services.error('Error executing cql3d prepare_input')
              raise Exception('Error executing cql3d prepare_input')
        
# ptb: Need to first copy the cqlinput_new file to cqlinput
          shutil.copyfile('cqlinput_new', 'cqlinput')

# Check if this is a pwrscale iteration such that pwrscale needs to be reset in cqlinput
          if 'pwrscale_arg' in kwargs:
            pwrscale = kwargs.get('pwrscale_arg')
            self.change_cql3d_pwrscale_dbb(pwrscale)

#     Launch cql3d - N.B: Path to executable is in config parameter CQL3D_BIN
          print('fp_cql3d: launching cql3d')
          cwd = services.get_working_dir()
          task_id = services.launch_task(self.NPROC, cwd, self.CQL3D_BIN, task_ppn=self.NPPN, logfile='log.cql3d')
          retcode = services.wait_task(task_id)

#          command = 'srun -N 12 -n 64  -c 6 /project/projectdirs/m77/CompX/CQL3D/master/xcql3d_mpi_intel.cori 2>>log.stdErr 1>>log.stdOut'
#          retcode = subprocess.call(command.split(), stdout = log_file,\
#                                    stderr = subprocess.STDOUT)
          if (retcode != 0):
              print('Error executing command: ', cql3d_bin)
              services.error('Error executing cql3d')
              raise

# If this is the first step in a pwrscale iteration, and restart = 'disabled' (i.e. this
# is first time cql3d has run) save the cql3d.nc file to initial_distrfunc.nc.   
          if 'icount_arg' in kwargs:
             icount = kwargs.get('icount_arg')
             if (icount == 1) and (restart == 'disabled'):
             	shutil.copyfile('cql3d.nc', 'initial_distrfunc.nc')

    # Call process_output - step
          print('fp_cql3d step: calling process_output')          

# Get cql3d_output_file file name <--> mnemonic from cqlinput file
          lines = get_lines('cqlinput')
          cql3d_output_file = lines_to_variable_dict(lines)['MNEMONIC'].strip("'") + ".nc"
          print('cql3d_output_file = ', cql3d_output_file)

          log_file = open('log_process_cql3d_output', 'w')
          mode = 'step'
          command = process_output_bin + ' ' +  cql3d_output+ ' ' +  cql3d_output_file   

          print('running', command)
          services.send_portal_event(event_type = 'COMPONENT_EVENT',\
              event_comment =  command)

          retcode = subprocess.call(command.split(), stdout = log_file,\
                                    stderr = subprocess.STDOUT)                                  
          if (retcode != 0):
              print('Error executing cql3d cql3d process_output ', process_output_bin)
              services.error('Error executing cql3d process_output')
              raise Exception('Error executing cql3d process_output')
          
# Copy generic cql3d partial plasma state file -> FP_CQL3D_cur_state_file  [correct??,
#          try:
#            partial_file = cwd + '/FP_CQL3D_' + cur_state_file
#            print('partial_file ', partial_file)
#            shutil.copyfile('FP_CQL3D_PARTIAL_STATE', partial_file )
#          except IOError as xxx_todo_changeme8:
#            (errno, strerror) = xxx_todo_changeme8.args
#            print('Error copying file %s to %s' % ('cur_state.cdf', cur_state_file, strerror))
#            raise


# Merge partial plasma state containing updated CQL3D data  [BH]
#          try:
#             services.merge_current_plasma_state(partial_file, logfile='log.update_state')
#             print('merged CQL3D plasma state data ', partial_file)
#          except Exception as e:
#             print('Error in call to merge_current_plasma_state(' , partial_file, ')')
#             self.services.error('Error in call to merge_current_plasma_state')
#             raise Exception('Error in call to merge_current_plasma_state')

  # Update plasma state files in plasma_state work directory, but only cur_cql_file
  # and cur_ImChizz_inp_file if there is one.
  # This way it can be used concurrently without overwriting other plasma state files.
          print('CURRENT_CQL = ', cur_cql_file)
          shutil.copyfile(cql3d_output_file,cur_cql_file)
          try:
             if cur_cql_file != None:
               services.update_plasma_state([cur_cql_file])
          except Exception:
             logMsg = 'Error in call to update_plasma_state()'
             self.services.exception(logMsg)
             raise 
            
    # Archive output files
          try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
          except Exception as e:
            print('Error in call to stage_output_files()', e)
            services.error('Error in call to stage_output_files()')
            raise Exception('Error in call to stage_output_files()')
        
#     rename the log file so that it is not appended next step
#        os.rename('log.cql3d', this_logfile)

          return 0

        return 0  # return on "zero" LHRF power condition
# ------------------------------------------------------------------------------
#
# checkpoint function
# Saves plasma state files to restart directory
#
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
      print('cql3d.checkpoint() called')
      if (self.services == None) :
         print('Error in cql3d: checkpoint(): No services')
         services.error('Error in cql3d: checkpoint(): No services')
         raise Exception('Error in cql3d: checkpoint(): No services')
      services = self.services
      services.save_restart_files(timestamp, self.RESTART_FILES)

      return 0
        
# ------------------------------------------------------------------------------
#
# finalize function
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print('cql3d.finalize() called')

# ------------------------------------------------------------------------------
#
# Internal utility functions
#
# ------------------------------------------------------------------------------

    # This doesn't seem to actually be used - ask Lee
    def fix_namelist(self, in_file, out_file):
        fd = open(in_file, 'r')
        line_list = fd.readlines()

        buf = ''.join(line_list)
        out_buf = ''

        in_string = False
        for c in buf:
            string = ''
            if in_string:
                if (c == '\''):
                    out_buf = out_buf + c
                    in_string = False
                elif (c == '\n'):
                    pass
                else :
                    out_buf = out_buf + c
            else:
                if (c == '\''):
                    in_string = True
                out_buf = out_buf + c
        ofd = open(out_file,'w')
        ofd.write(out_buf)
        ofd.close()
        fd.close()
        return

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
# ------------------------------------------------------------------------------

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

    # A utility to do simple editing of txt files on a line by line basis.
    #---------------------------------------------------------------------------------------
    # Open an input file and return the lines
    def get_lines(self, filename):
        file = open(filename, 'r')
        lines = file.readlines()
        file.close()
        return lines

    #---------------------------------------------------------------------------------------
    # Open an output file and write lines into it
    def put_lines(self, filename, lines):
        file = open(filename, 'w')
        file.writelines(lines)
        file.close()


    #---------------------------------------------------------------------------------------
    # Editing utilities
    #---------------------------------------------------------------------------------------

    def change_cql3d_pwrscale_dbb(self, pwrscale):
        print('change_cql3d_pwrscale: pwrscale = ', pwrscale)

        # get lines from namelist file
        inputLines = self.get_lines('cqlinput')
        var = 'PWRSCALE'
        var_line_number = self.find_var_line_number(inputLines, var)
        inputLines[var_line_number] = 'PWRSCALE        = ' + str(pwrscale) + ', 19*1.00000000000000,\n'
        self.put_lines('cqlinput', inputLines)
        return


    def find_var_line_number(self, lines, var):

        # Find the line in the namelist containing 'var = ' and get rid of newline if present.
        var_line_number = -1
        for i in range(len(lines)):
            line = lines[i]
            if '=' in line:
                split_line = line.split('=')
                #print 'split_line = ', split_line
                if (split_line[0].strip()).lower() == var.lower():
                    var_line_number = i

        if var_line_number == -1:
            message = 'Could not find variable ', var, ' in namelist lines'
            print(message)
            raise Exception(message)
        
        #print 'var_line_number = ', var_line_number
        #print 'lines[var_line_number] = ', lines[var_line_number]
        return var_line_number

