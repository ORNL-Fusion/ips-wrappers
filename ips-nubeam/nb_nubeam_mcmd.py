#! /usr/bin/env python
"""
MCMD version of NUBEAM component.  Adapted from Long-Poe Ku's nb_nubeam_pll.py

"""

# Working notes Version 3:  DBB 5/24/2010
# This version has the additions needed to allow checkpoint/restart.

# Working notes:  DBB 5/14/2010
#
# Note: This component assumes that the first two files listed in [NUBEAM] config parameter
# INPUT_FILES are nubeam_init_files.dat and nubeam_step_files.dat respectively although the
# file names can be simulation specific.  These are then copied to generic names 
# nubeam_init_files.dat and nubeam_step_files.dat respectively. This constrains the order that 
# the files are listed but not the names themselves.  In other components we have opted to
# put input files for different simulations in different branches of the directory tree and to
# use generic names for input files.  We might want to go back and do that here.
#
# Like we did for the TSC component I am changing the convention on expected naming of input 
# plasma state files.  Instead of requiring the plasma state files names to be of the form 
# <sim_name>_ps.cdf, I get the plasma state file names from the global config parameter 
# CURRENT_STATE etc.  These files are then copied to the generic name input_state.cdf which is
# the name nubeam_comp_exec expects. 

import sys
import os
import subprocess
import getopt
import shutil
import string
from  component import Component
from Scientific.IO.NetCDF import *
import time

class nubeam(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        self.firstTime  = True
        print 'Created %s' % (self.__class__)

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------


    def init(self, timeStamp):
        print 'nubeam.init() called'

        services = self.services
        workdir = services.get_working_dir()

    # Get global configuration parameters
        try:
            self.plasma_state_file = services.get_config_param('CURRENT_STATE')
        except:
            print 'error in getting global config parameters'
            self.services.error('error in getting global config parameters')
            raise Exception, 'error in getting global config parameters'
        cur_state_file = self.plasma_state_file

    # Get [NUBEAM] configuration parameters
        try:
            files = self.INPUT_FILES.split()
            suffix = self.INPUT_SUFFIX
            nubeam_exec = self.NUBEAM_BIN

        except:
            print 'error in getting NUBEAM config parameters'
            self.services.error('error in getting NUBEAM config parameter')
            raise Exception, 'error in getting NUBEAM config parameter'

    # Get input files  
        try:
            services.stage_input_files(self.INPUT_FILES)
        except Exception, e:
            print 'Error in call to stageInputFiles()' , e
            self.services.error('Error in call to stageInputFiles()')
            raise Exception, 'Error in call to stageInputFiles()'

    # Copy init_files and step_files lists to generic names.  Unless the filenames
    # are already generic. NB: assumes nubeam_init_files.dat and nubeam_step_files.dat
    # are the first two files in the NUBEAM INPUT_FILES list
        nubeam_init_files = files[0]
        nubeam_step_files = files[1]
        generic_file_name = "nubeam_init_files.dat"
        if nubeam_init_files != generic_file_name:
            shutil.copyfile(nubeam_init_files, generic_file_name)
        generic_file_name = "nubeam_step_files.dat"
        if nubeam_step_files != generic_file_name:
            shutil.copyfile(nubeam_step_files, generic_file_name)

    # copy PREACT from phys/nubeam/share
    # need phys_bin_root
        try:
            adas = self.ADAS
            print 'adas = ', adas
            os.environ['ADASDIR'] = adas
        except Exception:
            services.exeception('need an ADAS parameter in nubeam conf--use PREACT model')
        preact =  self.PREACT
        print "preact = ", preact
    # copy into working directory
        try:
            shutil.copytree(preact, os.path.join(workdir, "PREACT"))
        except:
            print 'PREACT directory not copied or already there.'

    # Copy plasma state files over to working directory
        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'Error in call to stage_plasma_state()' , e
            self.services.error('Error in call to stage_plasma_state()')
            raise Exception, 'Error in call to stage_plasma_state()'

    # Copy plasma state files over to generic file names
        try:
            shutil.copyfile(cur_state_file, "input_state.cdf")
        except IOError, (errno, strerror):
            print 'Error copying file %s in nubeam' % ("input_state.cdf", strerror)
            self.services.error('Error copying plasma state files over to generic file names')
            raise Exception, 'Error copying plasma state files over to generic file names'

        os.environ['NUBEAM_ACTION'] = 'init'
        os.environ['PREACTDIR'] = 'PREACT'
        os.environ['FRANTIC_INIT'] = '50'
        try:
            del os.environ['FRANTIC_ACTION']
        except:
            pass

        task_id = services.launch_task(1, workdir, nubeam_exec, logfile = 'log.init_nubeam')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            print 'Error executing command:  mpi_nubeam_comp_exec: init '
            raise Exception('Error executing command:  mpi_nubeam_comp_exec: init ')

    # Update plasma state
        services.merge_current_plasma_state("state_changes.cdf", logfile='log.update_state')

    # Archive output files
        try:
            shutil.copyfile('log.init_nubeam', 'log.nubeam') # fmwk expects log.nubeam as output
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
            print 'Error in call to stage_output_files()', e
            self.services.error('Error in call to stage_output_files()')
            raise Exception, 'Error in call to stage_output_files()'

        return

# ------------------------------------------------------------------------------
#
# RESTART function
# Gets restart files from restart directory
# Loads the global configuration parameters from the config file
#
# ------------------------------------------------------------------------------
        
    def restart(self, timeStamp):
        print 'nubeam.restart() called'

        services = self.services
        workdir = services.get_working_dir()
        
    # Get restart files listed in config file.        
        try:
            restart_root = services.get_config_param('RESTART_ROOT')
            restart_time = services.get_config_param('RESTART_TIME')
            services.get_restart_files(restart_root, restart_time, self.RESTART_FILES)
        except Exception, e:
            print 'Error in call to get_restart_files()' , e
            raise

    # Get global configuration parameters
        try:
            self.plasma_state_file = services.get_config_param('CURRENT_STATE')
        except:
            print 'nubeam restart: error in getting config parameters'
            self.services.error('error in getting config parameters')
            raise Exception, 'error in getting config parameters'

        self.nubeam_log = os.path.join(workdir, 'log.nubeam')

    # copy PREACT from phys/nubeam/share
    # need phys_bin_root
        try:
            adas = self.ADAS
            print 'adas = ', adas
            os.environ['ADASDIR'] = adas
        except Exception:
            services.exeception('need an ADAS parameter in nubeam conf--use PREACT model')
        preact =  self.PREACT
        print "preact = ", preact
    # copy into working directory
        try:
            shutil.copytree(preact, os.path.join(workdir, "PREACT"))
        except:
            print 'PREACT directory not copied or already there.'
        os.environ['PREACTDIR'] = 'PREACT'

        return 0
        
# ------------------------------------------------------------------------------
#
# STEP function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'nubeam.step() called'

        if (self.services == None) :
           print 'Error in nubeam: step (): No self.services'
           self.services.error('Error in nubeam: step (): No self.services')
           raise Exception, 'Error in nubeam: step (): No self.services'
        services = self.services
        cur_state_file = self.plasma_state_file
      
    # Copy plasma state files over to working directory
        try:
           services.stage_plasma_state()
        except Exception, e:
           print 'Error in call to stage_plasma_state()' , e
           self.services.error('Error in call to stage_plasma_state()')
           raise Exception, 'Error in call to stage_plasma_state()'

    # Copy plasma state files over to generic file names, unless file names are
    # already generic
        if cur_state_file != "input_state.cdf" :
            try:
                shutil.copyfile(cur_state_file, "input_state.cdf")
            except IOError, (errno, strerror):
                print 'Error copying file %s in nubeam' % ("input_state.cdf", strerror)
                self.services.error('Error copying plasma state files to generic file names')
                raise Exception, 'Error copying plasma state files to generic file names'

        ps = NetCDFFile(cur_state_file, 'r')
        t1 = ps.variables['t1'].getValue()
        t0 = ps.variables['t0'].getValue()
        print 't1 = ', t1, 't0 = ',t0
        print  self.NSTEP_INT + 'x' + str((t1-t0)/int(self.NSTEP_INT))
        print  'step value = ', (t1-t0)/int(self.NSTEP_INT)
        print 'current state ', cur_state_file

        ps.close()

        workdir = services.get_working_dir()
        nubeam_binary = os.path.join(self.BIN_PATH, 'nubeam_comp_exec')
        update_binary = os.path.join(self.BIN_PATH, 'update_state')
        sim_name = services.get_config_param("SIM_NAME")

    # remove stale state_changes file 
        cmd = "rm " + workdir + "/state_changes.cdf"
        os.popen(cmd)

    # Set nubeam envirnment variables
        os.environ['NUBEAM_ACTION'] = 'step'
        os.environ['NUBEAM_REPEAT_COUNT'] = self.NSTEP_INT + 'x' + str((t1-t0)/int(self.NSTEP_INT))
        os.environ['NUBEAM_POSTPROC'] = 'fbm_write'
        os.environ['STEPFLAG'] = 'TRUE'
        os.environ['NUBEAM_POSTPROC'] = 'summary_test'
        try:
            del os.environ['FRANTIC_INIT']
        except:
            pass

        os.environ['FRANTIC_ACTION'] = 'execute'

        nubeam_exec = self.NUBEAM_BIN

    # On Franklin and jaguar we have observed NUBEAM to fail randomly but run ok when the exact
    # same job is resubmitted.  The logic below will resubmit a failed NUBEAM job up to max_attempt
    # times to try to avoid this.
        done = False
        run_count = 0
        max_attempt = 2
        while not done:
            run_count += 1
            task_id = services.launch_task(self.NPROC, workdir, nubeam_exec,
                                      logfile = 'log.nubeam')
            retcode = services.wait_task(task_id)
            if (retcode != 0):
                print 'Error executing command:  nubeam step ', nubeam_exec
                if (run_count >= max_attempt):
                    raise Exception('Error executing command:  nubeam step : %s' %
                                               (nubeam_exec))
                else:
                    self.services.error('Failed execution %d of %s', 
                                        run_count, nubeam_exec)
                    time.sleep(60)
            else:
                done = True

    # Merge partial plasma state containing updated nubeam data
        try:
           partial_file = "state_changes.cdf"
           services.merge_current_plasma_state(partial_file, logfile='log.update_state')
           print 'merged NUBEAM plasma state data ', partial_file
        except Exception, e:
           print 'Error in call to merge_current_plasma_state(' , partial_file, ')'
           self.services.error('Error in call to merge_current_plasma_state')
           raise Exception, 'Error in call to merge_current_plasma_state'

    # Archive output files
        try:
           services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
           print 'Error in call to stage_output_files()', e
           self.services.error('Error in call to stage_output_files()')
           raise Exception, 'Error in call to stage_output_files()'

        return 0

# ------------------------------------------------------------------------------
#
# checkpoint function
#
# Saves restart files to restart directory
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print 'nubeam.checkpoint() called'
        services = self.services
        services.save_restart_files(timestamp, self.RESTART_FILES)
        

    def finalize(self, timestamp=0.0):
        print 'nubeam.finalize() called'

