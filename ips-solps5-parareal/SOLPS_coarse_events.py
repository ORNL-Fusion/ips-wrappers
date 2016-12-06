#! /usr/bin/env python

from  component import Component
import subprocess
import shutil
import glob
import os
from ipsExceptions import InsufficientResourcesException
import time

class CoarseTask(object):
    def __init__(self, 
                 comp, 
                 iteration, 
                 slice, 
                 in_files = None, 
                 out_files = None, 
                 new_slice = False,
                 working_dir = None,
                 cmd_args = None):
        self.comp = comp
        self.services = comp.services
        self.iter = iteration
        self.slice = slice
        self.new_slice = new_slice
        if in_files:
            self.in_files = [f for f in in_files]
        self.out_files = {}
        self.dependencies = {}   # To be filled when dependencies are satisfied 
        self.status = None
        self.task_id = None
        self.working_dir = working_dir
        self.cmd_args = cmd_args
        self.suffix = '%04d.%04d' % (self.iter, self.slice)
        self.submit_args = []
        self.first_in_slice = False

        print 'coarse cmd args are ', self.cmd_args
                
#        if (self.new_slice):
#            self.satisfy_dependency('CONVERGE', iter - 1, slice, False)
#           self.satisfy_dependency('COARSE', iter - 1, slice - 1, 'DUMMY')
#           self.satisfy_dependency('FINE', iter - 1, slice - 1, 'DUMMY')
        
        
    def satisfy_dependency(self, kind, iteration, slice, value):
        if kind not in ('COARSE', 'FINE', 'CONVERGE'):
            raise Exception('Error: Task %d:%d received wrong dependency type %s' % 
                            (self.iter, self.slice, kind))
        self.dependencies[kind, iteration, slice] = value
        print 'COARSETASK %d %d : Satisfied %s %d %d with %s' % \
            (self.iter, self.slice, str(kind), int(iteration), int(slice), str(value))
        
        return
    
    def set_first_in_slice(self):
        self.first_in_slice = True
    
    def prepare_input(self, c2f_bin = None, f2c_bin = None, correction_bin = None):
        prev_slice_coarse_out = self.dependencies['COARSE', self.iter, self.slice - 1]

        values_input = 'coarse_in_values.%s' % (self.suffix)
        deltas_input = 'coarse_in_deltas.%s' % (self.suffix)
        wall_input = 'coarse_in_wall.%s' % (self.suffix)
        
# We're skipping the call to prepare_input() for coarse(1,1) and doing the copying in the component
        if (self.iter == 1 or self.first_in_slice):
            print 'COARSETASK using coarse out', self.suffix, self.first_in_slice 
            shutil.copyfile(prev_slice_coarse_out['COARSE_OUT_VALUES_CORRECTED'], values_input)
            shutil.copyfile(prev_slice_coarse_out['COARSE_OUT_DELTAS_CORRECTED'], deltas_input)
            shutil.copyfile(prev_slice_coarse_out['COARSE_OUT_WALL_CORRECTED'], wall_input)
        else:
            converge_prev_slice = self.dependencies['CONVERGE', self.iter-1, self.slice-1]
            if converge_prev_slice:
                print 'COARSETASK using fine out', self.suffix, 'CONV_PREV = ', converge_prev_slice 
                fine_out_prev_slice = self.dependencies['FINE', self.iter-1, self.slice-1]
                shutil.copyfile(fine_out_prev_slice['FINE_OUT_VALUES'], values_input)
                shutil.copyfile(fine_out_prev_slice['FINE_OUT_DELTAS'], deltas_input)
                shutil.copyfile(fine_out_prev_slice['FINE_OUT_WALL'], wall_input)
#                values_input = fine_out_prev_slice['FINE_OUT_VALUES']
#                deltas_input = fine_out_prev_slice['FINE_OUT_DELTAS']
            else:
                print 'COARSETASK using fine out', self.suffix, 'CONV_PREV = ',converge_prev_slice 
                coarse_out_prv_slice = self.dependencies['COARSE', self.iter, self.slice-1]
                shutil.copyfile(prev_slice_coarse_out['COARSE_OUT_VALUES_CORRECTED'], values_input)
                shutil.copyfile(prev_slice_coarse_out['COARSE_OUT_DELTAS_CORRECTED'], deltas_input)
                shutil.copyfile(prev_slice_coarse_out['COARSE_OUT_WALL_CORRECTED'], wall_input)
#                values_input =  prev_slice_coarse_out['COARSE_OUT_VALUES_CORRECTED']
#                deltas_input =  prev_slice_coarse_out['COARSE_OUT_DELTAS_CORRECTED']
            #this file now needs to be coarsened
#            subprocess.call([f2c_bin, prev_slice_output_values, my_input,
#                                 '1', '192', '192', '100', '100'])
        
        values_out = 'coarse_out_values.' + self.suffix
        deltas_out = 'coarse_out_deltas.' + self.suffix
        wall_out = 'coarse_out_wall.' + self.suffix
        plot_out = 'coarse_err_profile.%s' % (self.suffix)
#        pr_slice = str(self.slice)
        pr_path=str(self.suffix)
        args = [pr_path, values_input, values_out, deltas_input, deltas_out, wall_input, wall_out]
#  plot_out, self.comp.TIME_INIT, self.comp.TIME_FINAL, pr_slice, self.comp.STEP_NUM]
        self.out_files['COARSE_OUT_VALUES'] = os.path.join(self.working_dir, values_out)
        self.out_files['COARSE_OUT_DELTAS'] = os.path.join(self.working_dir, deltas_out)
        self.out_files['COARSE_OUT_WALL'] = os.path.join(self.working_dir, wall_out)
        self.submit_args = args
        return
    
    def process_output(self, c2f_bin = None, f2c_bin = None, correction_bin = None):
        
        #  fill high-k modes with zeros for ps
#        coarse_zero_pad = 'coarse_zero_pad.'+ self.suffix 
        my_out = self.out_files
#        subprocess.call([c2f_bin, my_out, 
#                         coarse_zero_pad,
#                         '2', '192', '192', '100', '100'])
#        self.out_files['COARSE_ZERO_PAD'] = os.path.join(self.working_dir, coarse_zero_pad)
        
#        coarse_advance_padded = 'coarse_advance_padded.' + self.suffix
        coarse_values_corrected = 'coarse_out_values_corrected.' + self.suffix
        coarse_deltas_corrected = 'coarse_out_deltas_corrected.' + self.suffix
        coarse_wall_corrected = 'coarse_out_wall_corrected.' + self.suffix
        if (self.iter == 1 or self.first_in_slice):
#            next_slice_coarse_in = 'coarse_in.' + '%04d.%04d' % (self.iter, self.slice + 1)
#            shutil.copyfile(my_out, next_slice_coarse_in)
#            try:
#                print 'coarse to fine in first iteration'
#                call_list = [c2f_bin, my_out, coarse_advance_padded, '1', '192',
#                 '192', '100', '100']
#                print 'call list = ', call_list
#                subprocess.call(call_list)
#            except:
#                self.services.exception('subprocess call failed in c2f_bin(1)')
#                raise
            shutil.copyfile(my_out['COARSE_OUT_VALUES'], coarse_values_corrected)
            shutil.copyfile(my_out['COARSE_OUT_DELTAS'], coarse_deltas_corrected)
            shutil.copyfile(my_out['COARSE_OUT_WALL'], coarse_wall_corrected)
            self.out_files['COARSE_OUT_VALUES_CORRECTED'] = os.path.join(self.working_dir, 
                                                                   coarse_values_corrected)
            self.out_files['COARSE_OUT_DELTAS_CORRECTED'] = os.path.join(self.working_dir, 
                                                                   coarse_deltas_corrected)
            self.out_files['COARSE_OUT_WALL_CORRECTED'] = os.path.join(self.working_dir, 
                                                                   coarse_wall_corrected)
            return
        
        #now with the three files, call the parareal advance
        prior_iter_coarse_out = self.dependencies['COARSE', self.iter-1, self.slice]
        prior_iter_fine_out = self.dependencies['FINE', self.iter-1, self.slice]
        # Need to add an output file spec to the command line arguments. 
        run_num=str(self.suffix)
        #Calculating number of lines-SD
#        numlines = len(open(my_out['COARSE_OUT_VALUES']).readlines())
        cmd = [correction_bin, 
                         prior_iter_fine_out['FINE_OUT_VALUES'],
                         my_out['COARSE_OUT_VALUES'],
                         prior_iter_coarse_out['COARSE_OUT_VALUES'], 
#                         coarse_values_corrected,run_num]
#                         str(numlines)]
                         prior_iter_fine_out['FINE_OUT_DELTAS'],
                         my_out['COARSE_OUT_DELTAS'],
                         prior_iter_coarse_out['COARSE_OUT_DELTAS'],
                         prior_iter_fine_out['FINE_OUT_WALL'],
                         coarse_values_corrected,
                         coarse_deltas_corrected,
                         coarse_wall_corrected,
                         run_num]
#                         str(numlines)]
#       Note, for PR correction, b2.wall_save.parameters is not corrected. The fine output is simply read in.
        print 'Executing CORRECTION command : ', ' '.join(cmd)
        subprocess.call(cmd, stdout = open('coarse_correct.%s.log' % (self.suffix), 'w'), stderr=subprocess.STDOUT)
        self.out_files['COARSE_OUT_VALUES_CORRECTED'] = os.path.join(self.working_dir, 
                                                                   coarse_values_corrected)
        self.out_files['COARSE_OUT_DELTAS_CORRECTED'] = os.path.join(self.working_dir, 
                                                                   coarse_deltas_corrected)
        self.out_files['COARSE_OUT_WALL_CORRECTED'] = os.path.join(self.working_dir, 
                                                                   coarse_wall_corrected)
        return
            
    def get_cmd_args(self):
        return self.submit_args
    
    def ready_to_run(self):
        print self.iter, self.comp.max_iter
        if self.iter > self.comp.max_iter:
            print self.suffix, 'Exceeded max_iter - Will not run'
            return False
        
        if len(self.dependencies) == 6:
            print self.suffix, 'len(self.dependencies) == 6'
            try:
                prior_slice_converge = self.dependencies['CONVERGE', self.iter-1, self.slice - 1]
                current_slice_converge = self.dependencies['CONVERGE', self.iter-1, self.slice]
            except KeyError:
                pass
            else:
                if prior_slice_converge and not current_slice_converge:
                    self.satisfy_dependency('COARSE', self.iter, self.slice-1, None)
                    
        if len(self.dependencies) == 7:
            print self.suffix, 'len(self.dependencies) == 7'
            print self.suffix, self.dependencies
            if self.dependencies['CONVERGE', self.iter-1, self.slice] == False:
                print 'COARSETASK %d %d ready to run' % (self.iter, self.slice)
                return True
        print 'COARSETASK %d %d NOT ready to run - dep = %s' % \
                (self.iter, self.slice, str(self.dependencies))
        return False
        
    def get_external_outputs(self):
        return self.out_files
    
class Coarse(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        self.event_arrived = False
        self.events_received = []
        self.ready_tasks = []
        self.task_table = {}    # Keyed by (iteration, slice)
        self.active_tasks = {}  # Keyed by taskid
        self.pipeline_width = int(self.services.get_config_param('NT_SLICE'))
        self.max_iter = int(self.services.get_config_param('MAX_ITERATION'))
        self.total_steps = int(self.TOTAL_STEPS)
        self.step_width = float(self.STEP_WIDTH)
#        self.time_init = float(self.TIME_INIT)
#        self.time_final = float(self.TIME_FINAL)
        self.solver_executable = self.EXECUTABLE
        self.total_steps = 0
        self.step_width = 0
        self.solver_executable = None

        try:
            self.max_slice = int(self.services.get_config_param('MAX_SLICE'))
        except KeyError:
            self.max_slice = self.pipeline_width
        self.nt_slice= {}
        for iter in range(1, self.max_iter+1):
            self.nt_slice[iter] = self.pipeline_width
        self.cmd_args = []

    def init(self, timeStamp = 0.0):
        self.working_dir = self.services.get_working_dir()
        self.total_steps = int(self.TOTAL_STEPS)
        self.step_width = float(self.STEP_WIDTH) 
#        self.time_init = int(self.TIME_INIT)
#        self.time_final = float(self.TIME_FINAL)
        self.solver_executable = self.EXECUTABLE
        self.services.subscribe('PARAREAL_FINE', self.handle_fine_events)
        self.services.subscribe('PARAREAL_CONVERGE', self.handle_converge_events)
        # Build task table, satisfying "pre-existing" dependencies 
        # for itertaion 0 coarse tasks
        c2f_bin = self.C2F_BIN
        f2c_bin = self.F2C_BIN
        correction_bin = self.CORRECTION_BIN
        iteration = 1
        slice = 1
        self.services.stage_input_files(self.INPUT_FILES)
                        
        coarse_out_dict={}
        coarse_out_dict['COARSE_OUT_VALUES'] = os.path.join(self.working_dir, self.INPUT_VALUES)
        coarse_out_dict['COARSE_OUT_DELTAS'] = os.path.join(self.working_dir, self.INPUT_DELTAS)
        coarse_out_dict['COARSE_OUT_WALL'] = os.path.join(self.working_dir, self.INPUT_WALL)
        coarse_out_dict['COARSE_OUT_VALUES_CORRECTED'] = os.path.join(self.working_dir, self.INPUT_VALUES)
        coarse_out_dict['COARSE_OUT_DELTAS_CORRECTED'] = os.path.join(self.working_dir, self.INPUT_DELTAS)
        coarse_out_dict['COARSE_OUT_WALL_CORRECTED'] = os.path.join(self.working_dir, self.INPUT_WALL)
        
        task_11 = CoarseTask(self, iteration, slice, 
                             working_dir=self.working_dir, 
                             cmd_args = self.cmd_args)
        task_11.satisfy_dependency('COARSE', iteration, slice-1, coarse_out_dict)
        task_11.satisfy_dependency('COARSE', iteration-1 , slice-1, None)
        task_11.satisfy_dependency('COARSE', iteration-1 , slice, None)
        task_11.satisfy_dependency('FINE', iteration-1, slice-1, None)
        task_11.satisfy_dependency('FINE', iteration-1, slice, None)
        task_11.satisfy_dependency('CONVERGE', iteration-1, slice, False)
        task_11.satisfy_dependency('CONVERGE', iteration-1, slice-1, False)
        self.task_table[iteration, slice] = task_11
        self.ready_tasks.append(task_11)
        return

    def step(self, timeStamp = 0.0):
        
        coarse_bin = self.EXECUTABLE 
        done = False
        services = self.services
        curdir = self.services.get_working_dir()
        converged_slices = []
        c2f_bin = self.C2F_BIN
        f2c_bin = self.F2C_BIN
        correction_bin = self.CORRECTION_BIN
        while not done:
            while len(self.ready_tasks) > 0:
                new_task = self.ready_tasks.pop(0)
                new_task.prepare_input(c2f_bin, f2c_bin, correction_bin)
                args = new_task.get_cmd_args()
                log_file_name =  'coarse.' + new_task.suffix + '.log'
                try:
                    task_tag = new_task.suffix
                    task_id = services.launch_task(self.NPROC, curdir, coarse_bin, \
                                          *args, block=False, logfile=log_file_name, tag = task_tag)
                except InsufficientResourcesException:
                    self.ready_tasks.insert(0, new_task)
                    break
                except:
                    services.exception('Error launching task COARSE %04d.%04d'  , 
                                       new_task.iter, new_task.slice)
                    raise
                else:
                    self.active_tasks[task_id] = new_task
                
            self.event_arrived = False
            self.events_received = []
            services.process_events()
            for event in self.events_received:
                if event['KIND'] == 'FINE_FINISHED':
                    event_slice = event['SLICE']
                    event_iter = event['ITER']
                    event_fine_outfiles = event['OUT_FILES']
                    self.update_create_task(task_iter = event_iter + 1,
                                            task_slice = event_slice,
                                            dep_kind = 'FINE', 
                                            dep_iter = event_iter,
                                            dep_slice = event_slice, 
                                            dep_value = event_fine_outfiles)
#                    if event_slice == self.nt_slice:
#                        continue
                    if (event_slice + 1 in converged_slices):
                        continue
                    self.update_create_task(task_iter = event_iter + 1,
                                            task_slice = event_slice + 1,
                                            dep_kind = 'FINE', 
                                            dep_iter = event_iter,
                                            dep_slice = event_slice, 
                                            dep_value = event_fine_outfiles)
                elif event['KIND'] == 'SLICE_CONVERGED':
                    event_slice = event['SLICE']
                    event_iter = event['ITER']
                    converged_slices.append(event_slice)
                    self.update_create_task(task_iter = event_iter + 1, 
                                            task_slice = event_slice + 1, 
                                            dep_kind = 'CONVERGE', 
                                            dep_iter = event_iter,
                                            dep_slice = event_slice, 
                                            dep_value = True)
                    self.update_create_task(task_iter = event_iter + 1, 
                                            task_slice = event_slice, 
                                            dep_kind = 'CONVERGE', 
                                            dep_iter = event_iter,
                                            dep_slice = event_slice, 
                                            dep_value = True)
                    for iter in range(event_iter+1, self.max_iter + 1):
                        self.nt_slice[iter] = min(self.nt_slice[iter] + 1, 
                                                  self.max_slice)
                        print 'Set COARSE NT_SLICE[%d] to %d' % (iter, self.nt_slice[iter])
                elif event['KIND'] == 'SLICE_NOT_CONVERGED':
                    target_slice = event['SLICE']
                    target_iter = event['ITER']
                    self.update_create_task(task_iter = target_iter + 1, 
                                            task_slice = target_slice, 
                                            dep_kind = 'CONVERGE', 
                                            dep_iter = target_iter,
                                            dep_slice = target_slice, 
                                            dep_value = False)
                    self.update_create_task(task_iter = target_iter + 1, 
                                            task_slice = target_slice + 1, 
                                            dep_kind = 'CONVERGE', 
                                            dep_iter = target_iter,
                                            dep_slice = target_slice, 
                                            dep_value = False)
                elif event['KIND'] == "ALL_CONVERGED" or event['KIND'] == "ALL_ABORT" :
                    done = True
                    break
            if done:
                break
            
            finished_tasks = self.services.wait_tasklist(self.active_tasks.keys(), block = False)
            while len(finished_tasks) > 0:
                task_id, retval = finished_tasks.popitem()
                task = self.active_tasks[task_id]
                del self.active_tasks[task_id]
                task.process_output(c2f_bin, f2c_bin, correction_bin)
                new_event = {}
                new_event['KIND'] = 'COARSE_FINISHED'
                new_event['SLICE'] = task.slice
                new_event['ITER'] = task.iter
                out_files = task.get_external_outputs()
                new_event['OUT_FILES'] = out_files
                new_event['EXIT_STATUS'] = retval
                self.services.publish('PARAREAL_COARSE', 'COARSE_FINISHED', new_event)
                
                self.update_create_task(task_iter = task.iter + 1, 
                                        task_slice = task.slice, 
                                        dep_kind = 'COARSE', 
                                        dep_iter = task.iter,
                                        dep_slice = task.slice, 
                                        dep_value = out_files)
                
                #if task.slice == self.nt_slice:
                #    continue
                
                self.update_create_task(task_iter = task.iter, 
                                        task_slice = task.slice + 1, 
                                        dep_kind = 'COARSE', 
                                        dep_iter = task.iter,
                                        dep_slice = task.slice, 
                                        dep_value = out_files)
                
                self.update_create_task(task_iter = task.iter + 1, 
                                        task_slice = task.slice + 1, 
                                        dep_kind = 'COARSE', 
                                        dep_iter = task.iter,
                                        dep_slice = task.slice, 
                                        dep_value = out_files)
                                                    
            time.sleep(0.1)
        return
    
    def finalize(self, timeStamp = 0.0):
        pass
    
    def update_create_task(self, task_iter = 0, 
                                 task_slice = 0, 
                                 dep_kind = None, 
                                 dep_iter = 0, 
                                 dep_slice = 0, 
                                 dep_value= None):
        try:
            affected_task = self.task_table[task_iter, task_slice]
        except KeyError:
            affected_task = CoarseTask(self, task_iter, task_slice,
                                       working_dir=self.working_dir, 
                                       cmd_args = self.cmd_args)
            self.task_table[task_iter, task_slice] = affected_task
            if (task_iter == 1 or task_slice > self.nt_slice[task_iter - 1]):
                affected_task.satisfy_dependency('CONVERGE', task_iter-1, task_slice, False)
                affected_task.satisfy_dependency('CONVERGE', task_iter-1, task_slice-1, False)            
                affected_task.satisfy_dependency('FINE', task_iter-1, task_slice-1, None)            
                affected_task.satisfy_dependency('FINE', task_iter-1, task_slice, None)            
                affected_task.satisfy_dependency('COARSE', task_iter-1, task_slice-1, None)            
                affected_task.satisfy_dependency('COARSE', task_iter-1, task_slice, None)            
                affected_task.set_first_in_slice()
                
        affected_task.satisfy_dependency(dep_kind, dep_iter, 
                                         dep_slice, dep_value)
        if task_iter > self.max_iter:
            return
        
        if task_slice <= self.nt_slice[task_iter]:
            if affected_task.ready_to_run():
                print 'Adding COARSETASK %d %d to ready queue' % (task_iter, task_slice)
                self.ready_tasks.append(affected_task)
        return
    
    def handle_fine_events(self, topicName, theEvent):
        self.event_arrived = True
        self.events_received.append(theEvent.getBody())
        return
    
    def handle_converge_events(self, topicName, theEvent):
        return self.handle_fine_events(topicName, theEvent)
                    
