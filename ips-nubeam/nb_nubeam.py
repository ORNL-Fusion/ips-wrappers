#! /usr/bin/env python
"""
Version 4: Batchelor (5/27/2013)

This version allows the power for multiple beams to be set and changed over time from
within the simulation config file.  If this option is used, any programming of the 
beam power from the EPA component is over-ridden.  This behavior is triggered by the
presence a variable N_BEAMS_PROGRAMMED in the [NUBEAM] section of the simulation config 
file.  If N_BEAMS_PROGRAMMED = 0, or is absent, the beam powers are taken from the current
plasma state file.  If > 0 the component looks for space delimited lists of time points
for power changes and the associated powers.  The lists should be named BEAM1_TIMES,
BEAM1_POWERS_MW, BEAM2_TIMES, BEAM2_POWERS_MW, etc.  N_BEAMS_PROGRAMMED must match the number
of beams allocated in the plasma state (ps%nbeams).  The first time the simulation time at
the beginning of a time step (ps%t0) is >= BEAMx_TIME then BEAMx_POWER is changed.
BEAMx_TIMES don't have to match simulation time step beginnings, but the changes won't  
take place until the beginning of the next time step.

Note: for ease of typing enter powers in config file in MW, they get converted to Watts in
this component.

The mechanism for changing the power is to modify the plasma state just before launching
the nubeam step.  The modification is done by a small FORTRAN code called
set_ps_beam_power.f90.

This version was obtained by extending the MCMD version of NUBEAM component, which was 
itself adapted from Long-Poe Ku's nb_nubeam_pll.py

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

#----------------------------------------------------------------------------------------------        
# Utility functions
#----------------------------------------------------------------------------------------------


    def  piecewise_constant(self, x, x_points, y_points):
    
    # Piecewise constant function. 
    # x = function argument
    # x_points = monotonically increasing list of x values
    # y_points = list of function values at the x_points, must be same length as x_points
    # 
    # If x < x_points[0] returns zero
    # If x >= x_points[i] and x < x_points[i+1] returns y_points[i]
    # If x > x_points[-1] returns y_points[-1]

        eps = 10.0**(-10)
        if len(x_points) != len(y_points):
            raise Exception('piecewise_constant: len(x_points != len(y_points)')
              
        if x < x_points[0]:
            return 0.0
        elif x > x_points[-1] - eps:
            return y_points[-1]
        else:
            j_up = 1
            while x > x_points[j_up] - eps:  # Find the index of the point in x_points just above x 
                j_up = j_up + 1
            return y_points[j_up -1]      
        
#----------------------------------------------------------------------------------------------

    def time_points_and_powers(self, n_beams):
    
    # Does initializations for set_ps_power() defined below.
    # Gets the BEAMx_TIME and BEAMx_POWER_MW from the [NUBEAM] section of the config file,
    # does a little checking for consistency and returns a 2 element list [times_list,
    # powers_list].  Where times_list and  powers_list are length = n_beams and contain
    # lists with the time points and associated powers for each of the beams.  Thus the
    # shape of the returned list is (2,n_beams, len(BEAMx_TIME)) for x = 1 to n_beams.
    # It also converts powers to Watts so as to be consistent with plasma state.
    
        times_list = []
        powers_list = []
        for i in range(n_beams):
            beam_name = 'BEAM' + str(i+1)

            # get time points for this beam
            name = beam_name + '_TIMES'
            time_expr = 'self.' + name + '.split()'
            try:
                str_time_points = eval(time_expr)
            except:
                message = 'error in getting NUBEAM config parameters ' + name
                print message
                self.services.error(message)
                raise Exception, message
            print name, ' = ', str_time_points

            # get powers for this beam
            name = beam_name + '_POWERS_MW'
            power_expr = 'self.' + name + '.split()'
            try:
                str_powers = eval(power_expr)
            except:
                message = 'error in getting NUBEAM config parameters ' + name
                print message
                self.services.error(message)
                raise Exception, message
            print name, ' = ', str_powers
                
            # check that times and powers have same length
            if len(str_time_points) != len(str_powers):
                message = 'error: ', beam_name, \
                ' time points and powers have different lengths' 
                print message
                self.services.error(message)
                raise Exception, message
                
            # Convert to float, convert powers to Watts, and append to full lists
            times_list.append([float(x) for x in str_time_points])
            powers_list.append([10.0**6*float(x) for x in str_powers])
            
        return [times_list, powers_list]


# ------------------------------------------------------------------------------
            
    def set_ps_power(self, ps_file, BIN_PATH, t, time_points_and_powers_list):
        services = self.services    
        times_list = time_points_and_powers_list[0]
        powers_list = time_points_and_powers_list[1]
        
        n_beams = len(times_list)
        set_power_list = []
        for i in range(n_beams):
            power_beam_i = self.piecewise_constant(t, times_list[i], powers_list[i])
            set_power_list.append(str(power_beam_i))
        
        set_power_bin = os.path.join(BIN_PATH, 'set_ps_beam_power')
        args = [set_power_bin, ps_file] + set_power_list
        print 'args = ', args

        try:
            retcode = subprocess.call(args) 
        except Exception:
            services.error(' error executing set_power_bin')
            raise Exception(' error executing set_power_bin')

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------


    def init(self, timeStamp):
        print 'nubeam.init() called'

        services = self.services
        workdir = services.get_working_dir()
        global power_programming, times_powers_list

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
            BIN_PATH = self.BIN_PATH
            files = self.INPUT_FILES.split()
            suffix = self.INPUT_SUFFIX
            nubeam_exec = self.NUBEAM_BIN

        except:
            print 'error in getting NUBEAM config parameters'
            self.services.error('error in getting NUBEAM config parameter')
            raise Exception, 'error in getting NUBEAM config parameter'

        # Get [NUBEAM] power programming configuration parameters, if present
        n_beams = 0
        power_programming = False
        try:
            n_beams = int(self.N_BEAMS_PROGRAMMED)
        except:
            print '\nCould not get beam programming from config file'

        if n_beams > 0: 
            power_programming = True
            times_powers_list = self.time_points_and_powers(n_beams)
        else:
            print '\nUsing beam power in plasma state'

    # Get input files  
        try:
            services.stage_input_files(self.INPUT_FILES)
        except Exception, e:
            print 'Error in call to stage_input_files()' , e
            self.services.error('Error in call to stage_input_files()')
            raise Exception, 'Error in call to stage_input_files()'

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
            services.stage_state()
        except Exception, e:
            print 'Error in call to stage_state()' , e
            self.services.error('Error in call to stage_state()')
            raise Exception, 'Error in call to stage_state()'

    # Copy plasma state files over to generic file names
        try:
            shutil.copyfile(cur_state_file, "input_state.cdf")
        except IOError, (errno, strerror):
            print 'Error copying file %s in nubeam' % ("input_state.cdf", strerror)
            self.services.error('Error copying plasma state files over to generic file names')
            raise Exception, 'Error copying plasma state files over to generic file names'

    # If power is programmed from config file set powers in plasma state
        if power_programming == True:
            # Get t0 from plasma state
            ps = NetCDFFile(cur_state_file, 'r')
            t0 = ps.variables['t0'].getValue()
            ps.close()
            # set powers those for time = t0
            self.set_ps_power("input_state.cdf", BIN_PATH, t0, times_powers_list)

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
        services.merge_current_state("state_changes.cdf", logfile='log.update_state')

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
        global power_programming, times_powers_list
        
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

    # Get [NUBEAM] power programming configuration parameters, if present.
    # Regenerate times_powers_list
        n_beams = 0
        power_programming = False
        try:
            n_beams = int(self.N_BEAMS_PROGRAMMED)
        except:
            print '\nCould not get beam programming from config file'

        if n_beams > 0: 
            power_programming = True
            times_powers_list = self.time_points_and_powers(n_beams)
        else:
            print '\nUsing beam power in plasma state'


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
        global power_programming, times_powers_list
      
    # Copy plasma state files over to working directory
        try:
           services.stage_state()
        except Exception, e:
           print 'Error in call to stage_state()' , e
           self.services.error('Error in call to stage_state()')
           raise Exception, 'Error in call to stage_state()'

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

    # Get component-specific configuration parameters. 
        try:
            BIN_PATH = self.BIN_PATH
            NPROC = self.NPROC
            nubeam_exec = self.NUBEAM_BIN
        except:
            print 'model_NB init: error getting component-specific config parameters'
            services.error('nb step: error getting component-specific\
            config parameters')
            raise Exception, 'nb step: error getting component-specific\
            config parameters'

        workdir = services.get_working_dir()
        #nubeam_binary = os.path.join(self.BIN_PATH, 'nubeam_comp_exec')
        #update_binary = os.path.join(self.BIN_PATH, 'update_state')
        sim_name = services.get_config_param("SIM_NAME")

    # remove stale state_changes file 
        cmd = "rm " + workdir + "/state_changes.cdf"
        os.popen(cmd)

    # Set nubeam environment variables
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

    # If power is programmed from config file set powers for time = t0 in plasma state
        if power_programming == True:            # set powers those for time = t0
            self.set_ps_power("input_state.cdf", BIN_PATH, t0, times_powers_list)



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
           services.merge_current_state(partial_file, logfile='log.update_state')
           print 'merged NUBEAM plasma state data ', partial_file
        except Exception, e:
           print 'Error in call to merge_current_state(' , partial_file, ')'
           self.services.error('Error in call to merge_current_state')
           raise Exception, 'Error in call to merge_current_state'

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

