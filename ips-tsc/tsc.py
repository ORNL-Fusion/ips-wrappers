#! /usr/bin/env python

"""

There are 3 use cases for this component:

1) Start a simulation from scratch.  TSC is the first component to 
"init", that is after the simulation "init" component, "minimal_state_init".  
It allocates all EPA arrays including all  plasma species, even the ones that technically
belong to other components - like minority or beam ions.  It writes initial EPA data 
only - initial equilibrium, initial geqdsk, initial plasma profiles.  This only requires TSC
input file.  It does not need a binary restart file "sprsin". The other components follow with
their 'init', allocating their arrays and writing their initial plasma state data.  Then step
through the time loop.

2) Start a new simulation from TSC input and binary restart files of a previous run.
TSC does "init" but takes the data from previous "sprsin file".  On "init" TSC 
writes only EPA data to plasma state.  This includes things like minority or beam species and
dimensions of RF sources and neutral beams (that technically belong to other components).   
Other components do their standard "init", possibly with grids or initial conditions different
from the simulation used to "init" from. Then step through the time loop.

3) Restart an existing simulation.  Don't do "init" at all but continue time 
stepping from a valid, previous time step. 

Note: In all of these cases when executing a "step" function TSC is in restart mode.  That is
      TSC takes all it's inputs from the binary restart file "sprsin" and the only thing
      it gets from the "input" file is the end time of the time step.

The interaction between this component script and the TSC ascii input file is now done by a
helper python script called TSC_input_file.py rather than by hard coding in this component 
script.  The helper script exports a group of functions that make it easy to manipulate a 
standard TSC input file.  The TSC_input_file.py module is imported and used to read and write  
TSC input files, set the time step end time, insert a number 99 card needed on restart, etc.
It also does some consistency checks during the 'init' phase and gives warnings if anything 
seems weird, such as start time of simulation != TSC start time, bootstrap current calculation 
turned off, energy transport calculation turned off, etc      

"""

# Working notes 8/13/2010:
# TSC now writes two partial states for merging by the concurrent framework.  These are:
# ps_update_state_eq.cdf and ps_update_state_pa.cdf.  This conponent now does
# services.merge_current_plasma_state() on these rather than services.update_plasma_state.
# However it still does an update_plasma_state with the reduced_ps_list as mentioned just
# below.
# Working notes 6/30/10:
# In 'step' added services.update_plasma_state(reduced_ps_list) where reduced_ps_list is
# the list of plasma state files written by epa but excluding the current plasma state so as
# not to step on the merged state. This list presently contains eqdsk_file and jso_file.
# The prior_state_file has been eliminated.  It was an error to be there.
# It should never have been updated except by the
# driver at the begnnning of the time step.  Also next_plasma_state has been eliminated. TSC
# is not writing the icrf and nbi power for the next time step in there, so it serves no 
# purpose.  Later we should communicate the power to the other components a different way.
#
# Working notes 5/24/10:
# This version has the additions needed to allow checkpoint/restart.
#
# Working notes:
#
#
# TSC interface: (Steve, Jin or Long-Poe please correct anything that's wrong)
# 
# command line arguments
# 1) sim_name - defaults to "I1234"
# 2) tokamak_name - defaults to "ITER"
# 3) year - defaults to "09"
# 4) shot_number - defaults to "10001"
# 5) mode - If there are 5 args and the 5th one = "init" it generates a plasma state and
#           quits. Therefore one cannot use the defaults when in a SWIM component.  As far as
#           I can tell it does not test for mode = "step" or "finalize"
#
# Input file name conventions:
# input.<sim_name> - Needed both for normal startup and restart
# config_nbi_iter.dat - It is hard coded into TSC that if the machine name is ITER it looks 
#                       for this file, independent of whether an NB component is being used.
#                       This probably should be generalized to look for other machines that 
#                       might have neutral beams, like DIII-D, NSTX and JET for example.
#
# sprsin.<sim_name> - Needed for restart. However, this needs to be listed as an input file  
#                     in the config file since it's needed for the step.  The framework will 
#                     look for it in the INPUT_DIR and will crash if it's not there.  
#                     Therefore a (perhaps dummy) file of this name needs to be in INPUT_DIR 
#                     even when starting from t = 0.
# wall_data - no prefix or suffix on this file
#
# Output file name conventions:
#
# sprsou.<sim_name>, output.<sim_name>, fort66.<sim_name>, eqflou.<sim_name>, 
# divhis.<sim_name>, osprsou.<sim_name>
# <sim_name>_ps.cdf, <sim_name>_psp.cdf, <sim_name>_psn.cdf, <sim_name>_ps.jso,
# <sim_name>_ps.geq
# 
#
# N.B. In the above <sim_name> is the first tsc command line argument.  I think its only 
#      function is to tell tsc what extension to expect on input.<sim_name> and
#      sprsin.<sim_name> and what extension to tack onto the sprsou.<sim_name> file.  In this 
#      component it could be anything, not necessarily the config file parameter SIM_NAME.  
#      We want the external names of input files to have a suffix that doesn't have to 
#      be SIM_NAME.  So for now I'm keeping it simple: using the TSC config parameter 
#      INPUT_SUFFIX for the <sim_name> command line argument to TSC. 
#
#      So:
#      TSC expects input file name = input.<INPUT_SUFFIX>
#      On restart TSC expects restart file name = sprsin.<INPUT_SUFFIX>
#      TSC writes binary output file sprsou.<INPUT_SUFFIX> which gets copied by this 
#      component to sprsin.<INPUT_SUFFIX>.  Similarly TSC writes output files 
#      output.<INPUT_SUFFIX>, fort66.<INPUT_SUFFIX>, eqflou.<INPUT_SUFFIX>, divhis.<INPUT_SUFFIX>,
#      osprsou.<INPUT_SUFFIX> but as of now these aren't copied to any new name.
#
#      TSC writes plasma state file <INPUT_SUFFIX>_ps.cdf which gets copied to <CURRENT_STATE>
#      where <CURRENT_STATE> is the global parameter in the simulation config file. The plasma 
#      state file is also copied to  the <PRIOR_STATE> and <NEXT_STATE> files as defined in the 
#      config file. Similarly TSC writes additional plasma state files <INPUT_SUFFIX>_ps.geq, 
#      <INPUT_SUFFIX>_ps.jso which are copied to <CURRENT_EQDSK> and <CURRENT_JSDSK> respectively.
#      Note: the framework doesn't require you to define <CURRENT_STATE> etc in the config  
#      file, it's just a convenience.  But this component will break if you don't define them 
#      there. 

import sys
import os
import subprocess
import getopt
import shutil
import string
import time
from  component import Component
from TSC_input_file import *

class tsc(Component):
    def __init__(self, services=None, config=None):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

# ------------------------------------------------------------------------------
#
# init function
#
# ------------------------------------------------------------------------------

    def init(self, timeStamp):
        print 'epa.init() called'

        services = self.services
        workdir = services.get_working_dir()

    # Get global configuration parameters
        try:
            mode = services.get_config_param('SIMULATION_MODE')
            self.run_id = services.get_config_param('RUN_ID')
            self.tokamak =services.get_config_param('TOKAMAK_ID')  
            self.shot_number = services.get_config_param('SHOT_NUMBER')
            self.sim_name = services.get_config_param('SIM_NAME')
            self.plasma_state_file = services.get_config_param('CURRENT_STATE')
            self.prior_state_file = services.get_config_param('PRIOR_STATE')
            self.next_plasma_state_file = services.get_config_param('NEXT_STATE')
            self.eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            self.jso_file = services.get_config_param('CURRENT_JSDSK')
            self.tsc_log = os.path.join(workdir, 'log.tsc')
        except:
            print 'epa_tsc: error in getting config parametefs'
            raise 

        sim_name = self.sim_name
        state_file = self.plasma_state_file
        prior_state_file = self.prior_state_file
        next_state_file = self.next_plasma_state_file
        eqdsk_file = self.eqdsk_file
        jso_file = self.jso_file
        reduced_ps_list = [eqdsk_file, jso_file]
        
        # get YEAR from TSC section of config file, if present, else use default
        self.year = '9999'
        try:
            num_year = self.YEAR
        except Exception, e:
            print 'Error getting config param YEAR, using default = ', self.year, ' ', e
        else: 
            self.year = str(num_year)
        
        try:
            services.stage_input_files(self.INPUT_FILES)
        except Exception, e:
            print 'Error in call to stageInputFiles()' , e

        # TSC input and binary restart files assumed of form
        # input<.suffix> and sprsin<.suffix>.  Get suffix
        suffix = self.INPUT_SUFFIX
        input_file_name = 'input.' + suffix
        print 'input_file_name = ', input_file_name

        # Read the tsc input file
        input_file = TSC_input_file(input_file_name)
        input_file.get_lines()
    
        # check for unusual stuff in input file
        input_file.check_settings(timeStamp)
 
    # Copy current, prior and next state over to working directory
        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'Error in call to services.stage_plasma_state()', e

        # Make backups of the starting point plasma state files
        try:
            shutil.copyfile(state_file, 'state-0.nc')
        except IOError, (errno, strerror):
            print 'Error creating backup of initial state file' % (strerror)
            raise


    # Launch TSC init
        try:
            cwd = services.get_working_dir()
            task_id = services.launch_task(self.NPROC, cwd, self.TSC_BIN, suffix,
                      self.tokamak, self.year, self.shot_number, 'init', logfile=self.tsc_log)
            ret_val = services.wait_task(task_id)
            if ret_val != 0:
                print 'epa_tsc: init: abnormal return from tsc'
                raise
        except:
            print 'error in calling epa/tsc comp init'
            raise 

    # So far we can't get TSC to return a non-zero stop value if it fails.  So look in the 
    # log.tsc file to try to determine if there was an abnormal exit.

        tsc_failed = any([True for l in open('log.tsc').readlines() if \
                          'abnormal exit' in l])
        if (tsc_failed):
           print 'epa_tsc: init: abnormal return from tsc'
           raise Exception('Error in tsc invocation')

    # Rename output plasma state file -> CURRENT_STATE. Also eqdsk_file and jso_file. Unless
    # state files already have the same prefix as INPUT_SUFFIX
        print 'suffix + "_ps.cdf" = ', suffix + "_ps.cdf"
        print 'state_file = ', state_file
        if suffix + "_ps.cdf" != state_file:
            try:
                subprocess.call(['mv', suffix + "_ps.cdf", state_file ])
                subprocess.call(['mv', suffix + "_ps.geq", eqdsk_file ])
                subprocess.call(['mv', suffix + "_ps.jso", jso_file ])
            except Exception, e:
                print 'epa_tsc: Error renaming %s ' %(state_file), e
                raise

    # Copy plasma state file to prior and next states
        try:
            shutil.copyfile(state_file, prior_state_file)
            shutil.copyfile(state_file, next_state_file)
        except IOError, (errno, strerror):
            print 'Error copying previous and next state files' % (strerror)
            raise

    # Update plasma state files in plasma_state work directory
        try:
            services.update_plasma_state()
        except Exception, e:
            print 'epa_tsc: Error in call to updatePlasmaState()', e
            raise
            
    # copy binary output file sprsou<.suffix> -> binary restart file sprsin<.suffix>
        try:
            shutil.copyfile('sprsou.' + suffix, 'sprsin.' + suffix)
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s' % ('sprsou.' + suffix,\
            'sprsin.' + suffix, strerror)
            raise

    # In TSC input file set restart flag to 1.0 and add 99 card after card 11 with flag 29
    # (i.e. the one that sets stop time)
        input_file.set_IRST1_restart()
 
    # write modified input file to input<.suffix>
        output_file = TSC_input_file(input_file_name)
        output_file.lines = input_file.lines
        output_file.put_lines()
        
    # Archive output files
        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
            print 'epa_tsc: Error in call to stage_output_files()', e
            raise

    # Clean up
        cmd = "rm " + "/ographa" + " " + "/osprsou"
        try:
            os.popen(cmd)
        except:
            print 'not an error for tsc init if osprsou cannnot be deleted???'

        return 0
        
# ------------------------------------------------------------------------------
#
# RESTART function
# Gets restart files from restart directory
# Loads the global configuration parameters from the config file
#
# ------------------------------------------------------------------------------
        
    def restart(self, timeStamp):
        print 'epa.restart() called'

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
            mode = services.get_config_param('SIMULATION_MODE')
            self.run_id = services.get_config_param('RUN_ID')
            self.tokamak =services.get_config_param('TOKAMAK_ID')  
            self.shot_number = services.get_config_param('SHOT_NUMBER')
            self.sim_name = services.get_config_param('SIM_NAME')
            self.plasma_state_file = services.get_config_param('CURRENT_STATE')
            self.prior_state_file = services.get_config_param('PRIOR_STATE')
            self.next_plasma_state_file = services.get_config_param('NEXT_STATE')
            self.eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            self.jso_file = services.get_config_param('CURRENT_JSDSK')
            self.tsc_log = os.path.join(workdir, 'log.tsc')
        except:
            print 'epa_tsc restart: error in getting config parametefs'
            raise 

        # get YEAR from TSC section of config file, if present, else use default
        self.year = '9999'
        try:
            num_year = self.YEAR
        except Exception, e:
            print 'Error getting config param YEAR, using default = ', self.year, ' ', e
        else: 
            self.year = str(num_year)

        return 0
        
# ------------------------------------------------------------------------------
#
# STEP function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        print 'epa.step() calld'

        services = self.services
        workdir = services.get_working_dir()
        
        sim_name = self.sim_name
        state_file = self.plasma_state_file
        prior_state_file = self.prior_state_file
        next_state_file = self.next_plasma_state_file
        eqdsk_file = self.eqdsk_file
        jso_file = self.jso_file
        suffix = self.INPUT_SUFFIX
        reduced_ps_list = [eqdsk_file, jso_file]

        # Copy current and prior state over to working directory
        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'epa_tsc: Error in call to stageCurrentPlasmaState()', e
            raise

        # Copy CURRENT_STATE -> <suffix> + _ps.cdf. Unless state file already 
        # has the same prefix as INPUT_SUFFIX
        if suffix + "_ps.cdf" != state_file:
            try:
                subprocess.call(['mv',state_file, suffix + "_ps.cdf",  ])
            except Exception, e:
                print 'epa_tsc: Error renaming %s ' %(state_file), e
                services.error('epa_tsc: Error renaming state file')
                raise Exception, 'epa_tsc: Error renaming state file'

    # OPEN TSC input file and set stop time = timeStamp
        
        # TSC input and binary restart files assumed of form
        # input<.suffix> and sprsin<.suffix>.
        input_file_name = 'input.' + suffix
        input_file = TSC_input_file(input_file_name)
        input_file.get_lines()
        input_file.set_end_time(float(timeStamp))
            
        # write modified input file to input<.suffix>
        output_file = TSC_input_file(input_file_name)
        output_file.lines = input_file.lines
        output_file.put_lines()

    # Launch TSC step
        try:
            cwd = services.get_working_dir()
            task_id = services.launch_task(self.NPROC, cwd, self.TSC_BIN, suffix,
                      self.tokamak, self.year, self.shot_number, 'step', logfile=self.tsc_log)
            ret_val = services.wait_task(task_id)
            if ret_val != 0:
                print 'epa_tsc: init: abnormal return from tsc'
                raise
        except:
            print 'error in calling epa/tsc step'
            raise            

    # So far we can't get TSC to return a non-zero stop value if it fails.  So look in the 
    # log.tsc file to try to determine if there was an abnormal exit.

        tsc_failed = any([True for l in open('log.tsc').readlines() if \
                          'abnormal exit' in l])
        if (tsc_failed):
           print 'epa_tsc: init: abnormal return from tsc'
           raise Exception('Error in tsc invocation')

    # Rename output plasma state file -> CURRENT_STATE. Also eqdsk_file and jso_file. Unless
    # state files already have the same prefix as INPUT_SUFFIX
        if suffix + "_ps.cdf" != state_file:
            try:
                subprocess.call(['mv', suffix + "_ps.cdf", state_file ])
                subprocess.call(['mv', suffix + "_ps.geq", eqdsk_file ])
                subprocess.call(['mv', suffix + "_ps.jso", jso_file ])
            except Exception, e:
                print 'epa_tsc: Error renaming %s ' %(state_file), e
                services.error('epa_tsc: Error renaming state files')
                raise Exception, 'epa_tsc: Error renaming state files'

    # Merge partial plasma state containing updated tsc data

        try:
           partial_file = "ps_update_state_eq.cdf"
           services.merge_current_plasma_state(partial_file, logfile='log.update_state')
           print 'merged TSC plasma state data ', partial_file
        except Exception, e:
           print 'Error in call to merge_current_plasma_state(' , partial_file, ')'
           self.services.error('Error in call to merge_current_plasma_state')
           raise Exception, 'Error in call to merge_current_plasma_state'

        try:
           partial_file = "ps_update_state_pa.cdf"
           services.merge_current_plasma_state(partial_file, logfile='log.update_state')
           print 'merged TSC plasma state data ', partial_file
        except Exception, e:
           print 'Error in call to merge_current_plasma_state(' , partial_file, ')'
           self.services.error('Error in call to merge_current_plasma_state')
           raise Exception, 'Error in call to merge_current_plasma_state'

    # Update plasma state files in plasma_state work directory, but excluding current_plasma_state
        try:
            services.update_plasma_state(reduced_ps_list)
        except Exception, e:
            print 'epa_tsc: Error in call to updatePlasmaState()', e
            raise
            
    # copy binary output file: sprsou<.suffix> -> binary restart file: sprsin<.suffix>
        try:
            shutil.copyfile('sprsou.' + suffix, 'sprsin.' + suffix)
        except IOError, (errno, strerror):
            print 'Error copying file %s to %s: %s' % ('sprsou' + suffix,\
            'sprsin' + suffix, strerror)
            raise

    # Archive output files
        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
            print 'Error in call to stage_output_files()', e
            raise

    #
    # remove all signal files and clean up
#     cmd = "rm " + workdir + "/*.old" + " " + workdir + "/transp_*.dat"
#     os.popen(cmd)
        cmd = "rm " + workdir + "/ographa" + " " + workdir + "/osprsou"
        os.popen(cmd)
        return 0


# ------------------------------------------------------------------------------
#
# checkpoint function
#
# Saves plasma state files to restart directory
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print 'epa_tsc.checkpoint() called'
        services = self.services
        services.save_restart_files(timestamp, self.RESTART_FILES)
        

# ------------------------------------------------------------------------------
#
# FINALIZE function
# As of now it does nothing
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print 'epa.finalize() called'
        return 0

