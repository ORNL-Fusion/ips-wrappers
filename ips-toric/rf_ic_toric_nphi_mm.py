#! /usr/bin/env python

"""
rf_ic_nphi_mm.py  2/26/2013

Multiple mode MCMD version of TORIC component. This version launches parallel toric runs
for multiple toroidal modes for each of multiple icrf sources.  The information on number
of sources, their frequencies, and number of toroidal modes to be run for each source is
obtained from the Plasma State variables freq_ic(), num_nphi(), and nphi().

INIT:

A modified version of do_toric_init_abr.f90 has been written called do_toric_init_mm.f90.
The modified code reads a plasma state file which has EPA data already initialized.  
Then it reads Plasma State namelist files containing ICRF machine description and shot
configuration namelists using the plasma state functions ps_mdescr_read() and 
ps_sconfig_read() and copies the new ICRF data to the current plasma state and saves it.  
The original section of do_toric_init does further initialization 
to allocate and fill the profile arrays as is normally done in the component init.
do_toric_init_mm.f90 runs once per simulation and generates the initial plasma state.

Also during init a  run sub-directory is created for each toric run and a mapping file
(run_map.txt) is generated which describes the mapping between icrf source number, nphi 
and the path to the sub directory.  The format of the mapping file is one field per line:

line 1 -> num_runs = total number of runs and sub-directories

followed by num_runs groups of 5 lines each:
 index_src = fortran index of icrf source
 src_freq  = frequency of that source
 index_nphi = fortran index of nphi
 nphi(index_nphi, index_src) = nphi value for that index_nphi
 path name to corresponding directory
 
STEP:

For each time step prepare_toric_input_abr.f90 runs once and produces the torica.inp file.
Then for each freq/nphi run the torica.inp file is modified to contain the proper
values of toric variables FREQCY and NPHI.  All input files are copied to the appropriate
sub-directory.  The toric run is added to the task pool 'rf_pool' 
using services.add_task() and the loop goes on to deal with the next freq/nphi run.  When
all tasks in rf_pool are complete execution continues with process_output.

"""

# Working notes:

# 2/26/2013 – There is presently a raise Exception after the rf_pool finishes because we
#             don't yet have multimode process output.  we need to remember to get rid of
#             this later.


import os
import subprocess
import shutil
from  component import Component
from netCDF4 import Dataset

run_map_filename = 'run_map.txt'
input_inp_filename = 'torica.inp'
input_file_list = ['torica.inp', 'equidt.data', 'equigs.data']

class toric(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        self.run_map = {}
        self.map_file = 'run_map.txt'

    def init(self, timeStamp=0):
        print 'toric.init() called'

        services = self.services
        workdir = services.get_working_dir()

        # Get global configuration parameters
        try:
            self.plasma_state_file = services.get_config_param('CURRENT_STATE')
            self.eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            self.toric_log = os.path.join(workdir, 'log.toric')
        except:
            print 'rf_ic_toric_mcmd: error in getting config parameters'
            self.services.error('rf_ic_toric_mcmd: error in getting config parameters')
            raise Exception('rf_ic_toric_mcmd: error in getting config parameters')

        cur_state_file = self.plasma_state_file
        toric_log = self.toric_log


        # Copy plasma state files over to working directory
        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'Error in call to stage_plasma_state()', e
            self.services.error('Error in call to stage_plasma_state()')
            raise Exception('Error in call to stage_plasma_state()')

        # Get input files
        try:
            services.stage_input_files(self.INPUT_FILES)
        except Exception, e:
            print 'Error in call to stageInputFiles()', e
            self.services.error('Error in call to stageInputFiles()')
            raise Exception('Error in call to stageInputFiles()')

        # Copy machine.inp_<suffix> to generic file name -> machine.inp if
        # there is a suffix
        try:
            suffix = self.INPUT_SUFFIX
            have_suffix = True
        # If suffix is not empty put an underscore in front of it.
            if len(suffix) > 0:
                print 'INPUT_SUFFIX = ', suffix
                suffix = '_' + suffix
        # If suffix is empty you don't really have one
            else:
                have_suffix = False
        except:
            have_suffix = False

        # If there is a non-empty suffix, copy to generic filename
        if have_suffix:
            try:
                shutil.copyfile('machine.inp' + suffix, 'machine.inp')
            except IOError, (errno, strerror):
                print 'Error copying file %s to %s : %s' % \
                ('machine.inp' + suffix, 'machine.inp', strerror)
                services.error('Error copying machine.inp_<suffix> -> machine.inp')
                raise Exception('Error copying machine.inp_<suffix> -> machine.inp')

        # run TORIC init
        do_input = os.path.join(self.BIN_PATH, 'do_toric_init_mm')
        retcode = subprocess.call([do_input, cur_state_file])
        if (retcode != 0):
            print 'Error in call to toric_init'
            self.services.error('Error in call to toric_init')
            raise Exception('Error in call to toric_init')

        # Debugging
#         ps = Dataset(cur_state_file)
#         power_ic = ps.variables['power_ic'][:]
#         print 'power_ic = ', power_ic
#         raise Exception, 'Intentional stop in rf_ic_nphi_mm.py after do_toric_init'


        ps = Dataset(cur_state_file)
        num_nphi = ps.variables['num_nphi'][:]
        nphi = ps.variables['nphi'][:]
        freq = ps.variables['freq_ic'][:]
        self.run_map = {}
        text = ''
        for i in range(freq.size):
            src_num_nphi = num_nphi[i]
            src_freq = freq[i]
            src_nphi = nphi[i]
            for j in range(src_num_nphi):
                nphi_j = src_nphi[j]
                dir_name = os.path.join(os.getcwd(), 'run_%d_%d' % (i+1, nphi_j))
                try:
                    os.mkdir(dir_name)
                except OSError, e:
                    if not os.path.isdir(dir_name):
                        self.services.error('Error creating directory %s : %s' % 
                                            (dir_name, str(e)))
                        raise
                text += '%d\n%g\n%d\n%d\n%s\n' % (i+1, src_freq, j+1, nphi_j, dir_name)
                self.run_map[i+1, nphi_j] = dir_name
        num_runs = len(self.run_map.keys())
        out_text = '%d\n' % (num_runs) + text
        open(self.map_file, 'w').write(out_text)
        
        # Update plasma state files in plasma_state work directory
        try:
            services.update_plasma_state()
        except Exception, e:
            print 'Error in call to update_plasma_state()', e
            self.services.error('Error in call to update_plasma_state()')
            raise Exception('Error in call to update_plasma_state()')

        # Archive output files
        # N.B.  do_toric_init does not produce a complete set of TORIC output
        #       files.  This causes an error in stage_output_files().  To
        #       solve this we generate a dummy set of output files here with
        #       system call 'touch'
        for file in self.OUTPUT_FILES.split():
            print 'touching ', file
            subprocess.call(['touch', file])
      # Now stage them
        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
            print 'Error in call to stage_output_files()', e
            self.services.error('Error in call to stage_output_files()')
            raise Exception, 'Error in call to stage_output_files()'

        return 0

# ------------------------------------------------------------------------------
#
# RESTART function
# Gets restart files from restart directory
# Loads the global configuration parameters from the config file
#
# ------------------------------------------------------------------------------

    def restart(self, timeStamp):
        print 'toric.restart() called'

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
            self.eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            self.toric_log = os.path.join(workdir, 'log.toric')
        except:
            print 'toric restart: error in getting config parameters'
            self.services.error('error in getting config parameters')
            raise Exception, 'error in getting config parameters'

        return 0

# ------------------------------------------------------------------------------
#
# STEP function
#
# ------------------------------------------------------------------------------

    def step(self, timeStamp):
        """Take a step for the toric component.  Really a complete run."""
        print 'toric.step() called'

        if (self.services == None) :
            print 'Error in toric: step (): No self.services'
            self.services.error('Error in toric: step (): No self.services')
            raise Exception, 'Error in toric: step (): No self.services'
        services = self.services

      # Copy plasma state files over to working directory
        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'Error in call to stage_plasma_state()' , e
            self.services.error('Error in call to stage_plasma_state()')
            raise Exception, 'Error in call to stage_plasma_state()'

#      # Get input files
#        try:
#            services.stage_input_files(self.INPUT_FILES)
#        except Exception, e:
#            print 'Error in call to stageInputFiles()' , e
#            self.services.error('Error in call to stageInputFiles()')
#            raise Exception, 'Error in call to stageInputFiles()'

      # Copy machine.inp_<suffix> to generic file name -> machine.inp if there is
      # a suffix
        try:
            suffix = self.INPUT_SUFFIX
            have_suffix = True
        # If suffix is not empty put an underscore in front of it.
            if len(suffix) > 0:
                print 'INPUT_SUFFIX = ', suffix
                suffix = '_' + suffix
        # If suffix is empty you don't really have one
            else:
                have_suffix = False
        except:
            have_suffix = False

        # If there is a non-empty suffix, copy to generic filename
        if have_suffix:
            try:
                shutil.copyfile('machine.inp' + suffix, 'machine.inp')
            except IOError, (errno, strerror):
                print 'Error copying file %s to %s' % ('machine.inp' + suffix,
                'machine.inp', strerror)
                services.error('Error copying machine.inp_<suffix> -> machine.inp')
                raise Exception, 'Error copying machine.inp_<suffix> -> machine.inp'

        prepare_input = os.path.join(self.BIN_PATH, 'prepare_toric_input_abr')
        process_output  = os.path.join(self.BIN_PATH, 'process_toric_output_mcmd_mm')
        zero_RF_IC_power  = os.path.join(self.BIN_PATH, 'zero_RF_IC_power')
        toric_bin = self.TORIC_BIN
        prepare_eqdsk  = self.GEQXPL_BIN

        cur_state_file = self.plasma_state_file
        cur_eqdsk_file = self.eqdsk_file
        toric_log = self.toric_log
        cwd = os.getcwd()

# Check if ICRF power is zero (or effectively zero).  If true don't run toric just
# run zero_RF_IC_power fortran code
        print 'cur_state_file = ', cur_state_file
        ps = Dataset(cur_state_file)
        power_ic = ps.variables['power_ic'][:]
        ps.close()

#         ps = NetCDFFile(cur_state_file, 'r')
#         power_ic = ps.variables['power_ic'].getValue()[0]
#         ps.close()
        print 'power = ', power_ic
        if(-0.02 < sum(power_ic) < 0.02):
            retcode = subprocess.call([zero_RF_IC_power, cur_state_file])
            if (retcode != 0):
                print 'Error executing ', prepare_input
                self.services.error('Error executing zero_RF_IC_power')
                raise Exception, 'Error executing zero_RF_IC_power'

            # N.B. zero_RF_IC_power does not produce a complete set of TORIC output
            #      files.  This causes an error in stage_output_files().  To
            #      solve this we generate a dummy set of output files here with
            #      system call 'touch'
            for file in self.OUTPUT_FILES.split():
                subprocess.call(['touch', file])

# Check if ICRF power_ic[0] is negative.  If true don't run toric just
# retain power from previous time step i.e. leave sources untouched in the state.
# However power_ic[0] needs to be reset back to positive

        elif( power_ic[0] < -0.02):
            print 'continuing power from previous time step'
            power_ic[0] = -power_ic[0]
            ps.variables['power_ic'].assignValue(power_ic)
            ps.close()

    # Or actually run TORIC

        else:

            # Call TORIC prepare_input to generate torica.inp
            retcode = subprocess.call([prepare_input, cur_state_file]) #, cur_eqdsk_file])
            if (retcode != 0):
                print 'Error executing ', prepare_input
                self.services.error('Error executing TORIC prepare_input')
                raise Exception, 'Error executing TORIC prepare_input'

            # Call xeqdsk_setup to generate eqdsk.out file
            print 'prepare_eqdsk', prepare_eqdsk, cur_eqdsk_file

            retcode = subprocess.call([prepare_eqdsk, \
                                       '@equigs_gen', '/g_filename='+cur_eqdsk_file,\
                                       '/equigs_filename=equigs.data'])
            if (retcode != 0):
                print 'Error in call to prepare_eqdsk'
                self.services.error('Error executing TORIC prepare_eqdsk')
                raise Exception, 'Error executing TORIC prepare_eqdsk'

            #-----------------------------------------------------------------------------
            # Modify torica.inp and copy input files to run directories
            
            # Get run_map and initial machine.inp data
            file = open(run_map_filename, 'r')
            run_map = file.readlines()
            file.close()
            num_runs = int(run_map[0])

            file = open(input_inp_filename, 'r')
            torica_inp = file.readlines()
            file.close()
            
            # Find the lines in torica.inp containing FREQCY and NPHI
            for i in range(len(torica_inp)):
                split_line = torica_inp[i].split()
                if len(split_line) > 0:
                    if split_line[0].lower() == 'freqcy':
                        freq_line_number = i
                    if split_line[0].lower() == 'nphi':
                        nphi_line_number = i

            tasks = {}
            cwd = self.services.get_working_dir()
            pool = self.services.create_task_pool('rf_pool')
            run_dirs=[]
            for i in range(num_runs):
                n_src = int(run_map[5*i+1])
                freq = float(run_map[5*i+2])
                index_nphi = int(run_map[5*i+3])
                nphi = int(run_map[5*i+4])
                path = run_map[5*i+5].rstrip()
                # Change FREQCY and NPHI in lines
                torica_inp[freq_line_number] =' FREQCY = ' + str(freq) + ',\n'
                torica_inp[nphi_line_number] =' NPHI = ' + str(nphi) + ',\n'

                # update machine.inp file
                file = open(input_inp_filename, 'w')
                file.writelines(torica_inp)
                file.close()
        
                # copy files to run directories
                for file in input_file_list:
                    try:
                        shutil.copy(file, path)
                    except OSError, e:
                        self.services.error('Error copying file %s : %s' %  (file, str(e)))
                        raise  
                run_dir = os.path.join(cwd, path)
                run_dirs.append(run_dir)
                task_log_file = os.path.join(cwd, toric_log + '.%d.%d' % (n_src, nphi))
                self.services.add_task('rf_pool', 'task_%d_%d' %(n_src, nphi), 
                                       self.NPROC, run_dir, toric_bin, 
                                       logfile = task_log_file )
                
            ret_val = self.services.submit_tasks('rf_pool')
            print 'ret_val = ', ret_val
            exit_status = self.services.get_finished_tasks('rf_pool')
            print exit_status
            self.services.remove_task_pool('rf_pool')

#uncommment to run process_toric output
#            raise Exception, 'Intentional stop in rf_ic_nphi_mm.py'
             
            #-----------------------------------------------------------------------------
            # Launch TORIC executable
#             print 'toric processors = ', self.NPROC
#             cwd = services.get_working_dir()
#             task_id = services.launch_task(self.NPROC, cwd, toric_bin, logfile=toric_log)
#             retcode = services.wait_task(task_id)
#             if (retcode != 0):
#                 print 'Error executing command: ', toric_bin
#                 self.services.error('Error executing TORIC')
#                 raise Exception, 'Error executing TORIC'

            # Call process_output
            # First rename default fort.* to expected names by component method as of toric5 r918 from ipp

            for run_dir in run_dirs:
                os.chdir(run_dir)
                os.rename('fort.9','toric_cfg.nc')
                os.rename('fort.21','toric.nc')
                os.chdir('../')

            retcode = subprocess.call([process_output, cur_state_file])
            if (retcode != 0):
                print 'Error executing',  process_output
                self.services.error('Error executing TORIC process_output')
                raise Exception, 'Error executing TORIC process_output'


# Merge partial plasma state containing updated IC data
        try:
            partial_file = cwd + '/RF_IC_' + cur_state_file
            services.merge_current_plasma_state(partial_file, logfile='log.update_state')
            print 'merged TORIC plasma state data ', partial_file
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
# Saves plasma state files to restart directory
#
# ------------------------------------------------------------------------------

    def checkpoint(self, timestamp=0.0):
        print 'rf_ic_toric.checkpoint() called'
        services = self.services
        services.save_restart_files(timestamp, self.RESTART_FILES)


# ------------------------------------------------------------------------------
#
# FINALIZE function
# As of now it does nothing
#
# ------------------------------------------------------------------------------

    def finalize(self, timestamp=0.0):
        print 'toric.finalize() called'
