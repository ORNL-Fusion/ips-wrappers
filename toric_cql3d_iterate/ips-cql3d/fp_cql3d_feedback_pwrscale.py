#! /usr/bin/env python

from builtins import str
from builtins import range
import sys
import os
import subprocess
import getopt
import shutil
import string
import h5py
import scipy.io as spio
import numpy as np
from  ipsframework import Component
from netCDF4 import *
from simple_file_editing_functions import put_lines, get_lines, lines_to_variable_dict

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
            cur_cql_file = services.get_config_param('CURRENT_CQL')
           
        except:
            print('fp_cql3d_general: error getting config parameters CURRENT_STATE CURRENT_EQDSK')
            services.error('fp_cql3d_general: error getting config parameters CURRENT_STATE CURRENT_EQDSK')
            raise Exception('rf_ic_gneray: error getting config parameters CURRENT_STATE CURRENT_EQDSK')

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
            CQL3D_BIN = self.CQL3D_BIN
            CQL3D_SPECS = services.get_config_param('SPECS') #SF now stored in global config param
            CQL3D_MODE = self.CQL3D_MODE
            CQL3D_NML = self.CQL3D_NML
            NSTEPS_STR = self.NSTEPS_STR
            DELTAT_STR = self.DELTAT_STR
        except:
            print('fp_cql3d_general init: error getting cql3d-specific config parameters')
            services.error('fp_cql3d_general: error getting cql3d-specific\
            config parameters')
            raise Exception('fp_cql3d_general: error getting cql3d-specific\
            config parameters')

        # enorm which is used here and in cql3d
        arg_enorm = self.try_get_config_param(services,'ENORM', optional = True, default = 'None')
        arg_eqover = self.try_get_config_param(services,'EQ_OVERRIDE', optional=True, default = 'n')
        
        # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception as e:
          print('Error in call to stage_state()' , e)
          services.error('Error in call to stage_state()')
          raise Exception('Error in call to stage_state()')
        
        # Get input files  
        try:
            services.stage_input_files(INPUT_FILES)
        except Exception as e:
            print('Error in call to stageInputFiles()' , e)
            services.error('Error in call to stageInputFiles()')
            raise Exception('Error in call to stageInputFiles()')

        try:
            shutil.copyfile(self.CQL3D_NML,'cqlinput')
        except IOError as xxx_todo_changeme4:
            (errno, strerror) = xxx_todo_changeme4.args
            print('Error copying file %s to %s' % (self.CQL3D_NML, 'cqlinput', strerror))
            services.error('Error copying cur_eqdsk_file -> eqdsk')
            raise Exception('Error copying cur_eqdsk_file -> eqdsk')
      
        # Copy current eqdsk file to generic name -> eqdsk
        if(arg_eqover == 'n'):
            try:
                shutil.copyfile(cur_eqdsk_file, 'eqdsk')
            except IOError as xxx_todo_changeme4:
                (errno, strerror) = xxx_todo_changeme4.args
                print('Error copying file %s to %s' % (cur_eqdsk_file, 'eqdsk', strerror))
                services.error('Error copying cur_eqdsk_file -> eqdsk')
                raise Exception('Error copying cur_eqdsk_file -> eqdsk')

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
        workdir = services.get_working_dir()

        # Copy plasma state files over to working directory
        try:
          services.stage_state()
        except Exception as e:
          print('Error in call to stage_state()' , e)
          services.error('Error in call to stage_state()')
          raise Exception('Error in call to stage_state()')

        # Get input files
        try:
            services.stage_input_files(self.INPUT_FILES)
        except:
            logMsg = 'Error in call to stageInputFiles()'
            self.services.exception(logMsg)
            raise 

        # Copy cql3d input namelist to cqlinput
        shutil.copyfile(self.CQL3D_NML,'cqlinput')

        # Get Component specific config params
        CQL3D_MODE = self.CQL3D_MODE
        NSTEPS_STR = self.NSTEPS_STR
        DELTAT_STR = self.DELTAT_STR
        
        # Get global config params
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_cql_file = services.get_config_param('CURRENT_CQL')
        fp_specs = services.get_config_param('SPECS')

        # optional global config params
        arg_enorm = self.try_get_config_param(services,'ENORM', optional = True, default = '0.0')
        arg_nsurf_FP = self.try_get_config_param(services,'NSURF_FP', optional = True, default = '60')
        arg_rhoFPlo = self.try_get_config_param(services,'RHO_FPLO', optional = True, default = '0.01')
        arg_rhoFPhi = self.try_get_config_param(services,'RHO_FPHI', optional = True, default = '0.95')

        # Get kwargs
        arg_rf_code = kwargs.get('rf_code','None')
        pwrscale_on = kwargs.get('pwrscale_on',False)
        icount = kwargs.get('icount',0)
        
        # Get 
        prepare_input_bin = os.path.join(self.BIN_PATH, 'prepare_cql3d_input')
        process_output_bin  = os.path.join(self.BIN_PATH, 'process_cql3d_output')
        cql3d_bin = os.path.join(self.BIN_PATH, 'cql3d')

        # Copy current eqdsk file to generic name -> eqdsk. eqover fixes very niche cql3d bugs sometimes
        arg_eqover = self.try_get_config_param(services,'EQ_OVERRIDE', optional=True, default='n')
        if (arg_eqover == 'n'):
            try:
                shutil.copyfile(cur_eqdsk_file, 'eqdsk') 
            except IOError as xxx_todo_changeme8:
                (errno, strerror) = xxx_todo_changeme8.args
                print('Error copying file %s to %s' % (cur_eqdsk_file, 'eqdsk', strerror))
                services.error('Error copying cur_eqdsk_file -> eqdsk')
                raise Exception('Error copying cur_eqdsk_file -> eqdsk')
        
        # Check if there is an existing cql3d distribution file and if there is
        # put the code in restart mode
        # if first powerscale iteration copy distribution function from first
        # powerscale iteration
        restart = 'disabled'
        if os.path.isfile('./cql3d.nc'):
            if pwrscale_on and icount==0:
                shutil.copyfile('cql3d.nc', 'cql3d.nc_ini')
            elif pwrscale_on:
                shutil.copyfile('cql3d.nc_ini', 'cql3d.nc')

            if os.stat('./cql3d.nc')[6] != 0:
                restart = 'enabled'
                shutil.copyfile('cql3d.nc','distrfunc.nc')
                 
        #if deltat_str .eq. 0 then use timestep ramp from GENRAY/CQL sims
        t0 = float(timeStamp)
        if(DELTAT_STR=='0'):
            if(t0 < 5):
               DELTAT_STR = '0.00001'
            elif(t0<10):
               DELTAT_STR = '0.0001'
            elif(t0<30):
               DELTAT_STR = '0.001'
            elif(t0<40):
               DELTAT_STR = '0.005'
            else:
               DELTAT_STR = '0.01'

        #if pwrscale on set up power rescaling
        print('pwrscale_on ',pwrscale_on)
        if pwrscale_on:
            pwrscale_f = h5py.File('pwrscale.hdf5','r+')
            pwrtarget_dset = pwrscale_f['pwrtarget']
            pwrtarget = np.array(pwrtarget_dset)
            pwrscale_dset = pwrscale_f['pwrscale']
            pwrscale = np.array(pwrscale_dset)
            pwrscale_cql3d_dset = pwrscale_f['pfrac_cql3d']
        else:
            pwrscale = np.full(2,1.0)

        if t0 == 0.0:
            norf_flag = 1
        else:
            norf_flag = 0
            
        #Make input prep namelist
        nml_lines = ['&cql3d_prepare_nml\n']
        nml_lines.append(' cur_state_file = \"' + cur_state_file + '\",\n')
        nml_lines.append(' cql3d_specs = \"' + fp_specs + '\",\n')
        nml_lines.append(' cql3d_mode = \"' + CQL3D_MODE + '\",\n')
        nml_lines.append(' rf_code =\"' + arg_rf_code + '\",\n')
        nml_lines.append(' restart = \"' + restart + '\",\n')
        nml_lines.append(' nsteps = ' + NSTEPS_STR + ',\n')
        nml_lines.append(' arg_deltat = ' + DELTAT_STR + ',\n')
        nml_lines.append(' arg_enorm = ' + arg_enorm + ',\n')
        nml_lines.append(' arg_nsurfFP = ' + arg_nsurf_FP + ',\n')
        nml_lines.append(' arg_rhoFPlo = ' + arg_rhoFPlo + ',\n')
        nml_lines.append(' arg_rhoFPhi = ' + arg_rhoFPhi + ',\n')
        nml_lines.append(' pscale = ' + str(pwrscale[0]) + ', ' + str(pwrscale[1]) + ',\n')
        nml_lines.append(' norf = ' + str(norf_flag) + ',\n')
        nml_lines.append('/\n')
        put_lines('cql3d_prepare.nml',nml_lines)
        
        # Call prepare_input - step
        print('fp_cql3d step: calling prepare_input')

        if not os.path.isfile(prepare_input_bin):
            logMsg = 'Cannot find toric prepare_input binary: ' + prepare_input_bin
            self.services.error(logMsg)
            raise Exception(logMsg)
        
        log_file = open('log_prepare_cql3d_input_step', 'w')
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                                   event_comment =  prepare_input_bin)

        retcode = subprocess.call(prepare_input_bin, stdout = log_file,\
                                    stderr = subprocess.STDOUT)
        if (retcode != 0):
            print('Error executing cql3d: ', prepare_input_bin)
            services.error('Error executing cql3d prepare_input')
            raise Exception('Error executing cql3d prepare_input')
        
        shutil.copyfile('cqlinput_new', 'cqlinput')

        #Actually run CQL3D
        print('fp_cql3d: launching cql3d')
        cwd = services.get_working_dir()
        cpp_cql = int(2 * int(128/int(self.NPPN)))      
        task_id = services.launch_task(self.NPROC, cwd, self.CQL3D_BIN, whole_nodes=True, \
                                         task_cpp = cpp_cql, task_ppn=self.NPPN, logfile='log.cql3d')
        retcode = services.wait_task(task_id, timeout=600.0,delay=15.0)
        if (retcode != 0):
           print('Error executing command: ', cql3d_bin)
           services.error('Error executing cql3d')
           raise

        #get the name of the cql3d output file
        lines = get_lines('cqlinput')
        
        cql3d_output_file_tmp = lines_to_variable_dict(lines)['MNEMONIC'].replace(',','').replace("'","")
        cql3d_output_file_tmp = cql3d_output_file_tmp.rstrip()
        cql3d_output_file = cql3d_output_file_tmp + ".nc"

        #if power rescaling read rf powers from cql3d file and save
        #fraction of power to pwrscale hdf5 file
        if pwrscale_on:
            cql_nc = spio.netcdf_file(cql3d_output_file,'r')
            cql_powers_int = np.copy(cql_nc.variables['powers_int'].data)
            cql_nc.close()
            
            if fp_specs == 'MIN':
                power_ic_min = cql_powers_int[-1,0,4]
                pwrscale_cql3d_dset[0] = power_ic_min/pwrtarget
            elif fp_specs == 'MIN+':
                power_ic_min = cql_powers_int[-1,0,4]
                power_ic_bulk = cql_powers_int[-1,1,4]
                pwrscale_cql3d_dset[0] = power_ic_min/pwrtarget
                pwrscale_cql3d_dset[1] = power_ic_bulk/pwrtarget
            print('pfrac cql3d',np.array(pwrscale_cql3d_dset))
            pwrscale_f.close()
            #update the state with the latest power scale
            try:
                services.update_state(['pwrscale.hdf5'])
            except Exception:
                logMsg = 'Error in call to update_state()'
                self.services.exception(logMsg)
                raise 
                
        # Call process_output - step
        print("running process output")
        nml_lines = ['&cql3d_process_nml\n']
        nml_lines.append(' cur_state_file = \"' + cur_state_file + '\",\n')
        nml_lines.append(' cql3d_specs = \"' + fp_specs + '\",\n')
        nml_lines.append(' cql3d_mode = \"' + CQL3D_MODE + '\",\n')
        nml_lines.append(' cql3d_output_file = \"' + cql3d_output_file + '\",\n')
        nml_lines.append(' /\n')
        put_lines('cql3d_process.nml',nml_lines)
        
        services.send_portal_event(event_type = 'COMPONENT_EVENT',\
                                   event_comment =  process_output_bin)

        log_file = open('log_process_cql3d_output_step', 'w')
        retcode = subprocess.call(process_output_bin, stdout = log_file,\
                                   stderr = subprocess.STDOUT)                                  
        if (retcode != 0):
            print('Error executing cql3d cql3d process_output ', process_output_bin)
            services.error('Error executing cql3d process_output')
            raise Exception('Error executing cql3d process_output')

        # Update plasma state files in plasma_state work directory
        # This way it can be used concurrently without overwriting other plasma state files.
        print('CURRENT_CQL = ', cur_cql_file)
        shutil.copyfile(cql3d_output_file,cur_cql_file)
        try:
           if cur_cql_file != None:
             services.update_state([cur_cql_file])
        except Exception:
             logMsg = 'Error in call to update_state()'
             self.services.exception(logMsg)
             raise 

        
        # Archive output files
        try:
          services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception as e:
          print('Error in call to stage_output_files()', e)
          services.error('Error in call to stage_output_files()')
          raise Exception('Error in call to stage_output_files()')

        return 0

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

    def try_get_config_param(self, services, param_name, optional=False, default=None):

        try:
            value = services.get_config_param(param_name)
            print(param_name, ' = ', value)
        except Exception:
            if optional:
                print('config parameter ', param_name, ' not found')
                value = default
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

