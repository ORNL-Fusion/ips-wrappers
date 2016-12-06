#! /usr/bin/env python

from  component import Component
import os
import glob
import subprocess
import shutil
from ipsExceptions import InsufficientResourcesException
import time

class FineTask(object):
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
        self.working_dir = working_dir
        self.cmd_args = cmd_args
        if in_files:
            self.in_files = [f for f in in_files]
        self.out_files = {}
        self.dependencies = {}   # To be filled when dependencies are satisfied 
        self.status = None
        self.task_id = None
        self.suffix = '%04d.%04d' % (self.iter, self.slice)
        self.ran_prepare_input = False
        self.submit_args = []
        self.first_in_slice = False
        if (self.new_slice):
            self.satisfy_dependency('CONVERGE', iter - 1, slice, False)
            self.satisfy_dependency('COARSE', iter - 1, slice - 1, 'DUMMY')
            self.satisfy_dependency('FINE', iter - 1, slice - 1, 'DUMMY')
        
    def satisfy_dependency(self, kind, iteration, slice, value):
        if kind not in ('COARSE', 'FINE', 'CONVERGE'):
            raise Exception('Error: FINE Task %d:%d received wrong dependency type %s' % 
                            (self.iter, self.slice, kind))
        self.dependencies[kind, iteration, slice] = value
        print 'FINETASK %d %d : Satisfied %s %d %d with %s' % \
            (self.iter, self.slice, str(kind), int(iteration), int(slice), str(value))
        return
    
    def set_first_in_slice(self):
        self.first_in_slice = True
    
    def prepare_input(self, c2f_bin = None, f2c_bin = None, correction_bin = None):
        
        if self.ran_prepare_input:
            return
        
        coarse_out_prv_slice = self.dependencies['COARSE', self.iter, self.slice-1]
        fine_in_values = 'fine_in_values_solps.' + self.suffix
        fine_in_deltas = 'fine_in_deltas_solps.' + self.suffix
        fine_in_wall = 'fine_in_wall_solps.' + self.suffix
        if (self.iter == 1 or self.first_in_slice):
            fine_in_values = 'fine_in_values_solps.' + self.suffix
            fine_in_deltas = 'fine_in_deltas.' + self.suffix
            fine_in_wall = 'fine_in_wall.' + self.suffix
            print 'coarse to fine in first invocation of fine for the current slice'
            if (self.slice == 1):
                shutil.copyfile(coarse_out_prv_slice['COARSE_OUT_VALUES_CORRECTED'], fine_in_values)
                shutil.copyfile(coarse_out_prv_slice['COARSE_OUT_DELTAS_CORRECTED'], fine_in_deltas)
                shutil.copyfile(coarse_out_prv_slice['COARSE_OUT_WALL_CORRECTED'], fine_in_wall)
            else: # Same case here, could be different of coarse out needs to be converted to fine in
                shutil.copyfile(coarse_out_prv_slice['COARSE_OUT_VALUES_CORRECTED'], fine_in_values)
                shutil.copyfile(coarse_out_prv_slice['COARSE_OUT_DELTAS_CORRECTED'], fine_in_deltas)
                shutil.copyfile(coarse_out_prv_slice['COARSE_OUT_WALL_CORRECTED'], fine_in_wall)
        else:
            converge_prev_slice = self.dependencies['CONVERGE', self.iter-1, self.slice-1]
            if converge_prev_slice:
                fine_out_prev_slice = self.dependencies['FINE', self.iter-1, self.slice-1]
                shutil.copyfile(fine_out_prev_slice['FINE_OUT_VALUES'],fine_in_values)
                shutil.copyfile(fine_out_prev_slice['FINE_OUT_DELTAS'],fine_in_deltas)
                shutil.copyfile(fine_out_prev_slice['FINE_OUT_WALL'],fine_in_wall)
#                fine_in_values = fine_out_prev_slice['FINE_OUT_VALUES']
#                fine_in_deltas = fine_out_prev_slice['FINE_OUT_DELTAS']
            else:
                shutil.copyfile(coarse_out_prv_slice['COARSE_OUT_VALUES_CORRECTED'],fine_in_values)
                shutil.copyfile(coarse_out_prv_slice['COARSE_OUT_DELTAS_CORRECTED'],fine_in_deltas)
                shutil.copyfile(coarse_out_prv_slice['COARSE_OUT_WALL_CORRECTED'],fine_in_wall)
#                fine_in_values = coarse_out_prv_slice['COARSE_OUT_VALUES_CORRECTED']
#                fine_in_deltas = coarse_out_prv_slice['COARSE_OUT_DELTAS_CORRECTED']
            
        fine_out_values = 'fine_out_values.' + self.suffix
        fine_out_nc = 'b2time_fine.' + self.suffix + '.nc'
        fine_out_deltas = 'fine_out_deltas.' + self.suffix
        fine_out_wall = 'fine_out_wall.' + self.suffix
        plot_out = 'fine_err_profile.%s' % (self.suffix)
#        pr_slice = str(self.slice)
        pr_path = str(self.suffix)
        args = [pr_path,fine_in_values, fine_out_values,fine_in_deltas,fine_out_deltas,fine_in_wall,fine_out_wall] 
#                plot_out, self.comp.TIME_INIT, self.comp.TIME_FINAL, pr_slice, self.comp.STEP_NUM]
        self.out_files = {'FINE_OUT_VALUES' : os.path.join(self.working_dir, fine_out_values),
                          'FINE_OUT_NC' : os.path.join(self.working_dir, fine_out_nc),
                          'FINE_OUT_DELTAS' : os.path.join(self.working_dir, fine_out_deltas),
                          'FINE_OUT_WALL' : os.path.join(self.working_dir, fine_out_wall)}
        print 'fine command args = \n', args
        self.submit_args = args
        self.ran_prepare_input = True
        return 
    
    def process_output(self, c2f_bin = None, f2c_bin = None, correction_bin = None):
        return
    
    def get_cmd_args(self):
        return self.submit_args
        
    def ready_to_run(self):
        
        if self.iter > self.comp.max_iter:
            return False
        if len(self.dependencies) == 3:
            try:
                prior_slice_converge = self.dependencies['CONVERGE', self.iter-1, self.slice - 1]
                current_slice_converge = self.dependencies['CONVERGE', self.iter-1, self.slice]
            except KeyError:
                pass
            else:
                if prior_slice_converge and not current_slice_converge:
                    self.satisfy_dependency('COARSE', self.iter, self.slice-1, None)
                    
        if len(self.dependencies) == 4:
            if self.dependencies['CONVERGE', self.iter-1, self.slice] == False:
                print 'FINETASK %d %d ready to run' % (self.iter, self.slice)
                return True
        print 'FINETASK %d %d NOT ready to run - dep = %s' % \
                (self.iter, self.slice, str(self.dependencies))
        return False
    
    def get_external_outputs(self):
        return self.out_files
    
class Fine(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        self.event_arrived = False
        self.events_received = []
        self.ready_tasks = []
        self.task_table = {}    # Keyed by (iteration, slice)
        self.active_tasks = {}  # Keyed by taskid
        self.max_iter = int(self.services.get_config_param('MAX_ITERATION'))
        self.pipeline_width = int(self.services.get_config_param('NT_SLICE'))
        try:
            self.max_slice = int(self.services.get_config_param('MAX_SLICE'))
        except KeyError:
            self.max_slice = self.pipeline_width
        self.nt_slice= {}
        for iter in range(1, self.max_iter+1):
            self.nt_slice[iter] = self.pipeline_width

        self.cmd_args = []

    def init(self, timeStamp = 0.0):
        self.services.subscribe('PARAREAL_COARSE', self.handle_coarse_events)
        self.services.subscribe('PARAREAL_CONVERGE', self.handle_converge_events)
        self.working_dir = self.services.get_working_dir()
        # Build task table, satisfying "pre-existing" dependencies 
        # for itertaion 0 coarse tasks
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
                        
        task_11 = FineTask(self, iteration, slice, 
                           working_dir = self.working_dir,
                           cmd_args = self.cmd_args)
                
        task_11.satisfy_dependency('COARSE', iteration, slice-1, coarse_out_dict)
        task_11.satisfy_dependency('FINE', iteration-1, slice-1, None)
        task_11.satisfy_dependency('CONVERGE', iteration-1, slice, False)
        task_11.satisfy_dependency('CONVERGE', iteration-1, slice-1, False)
        
        self.task_table[iteration, slice] = task_11
        self.ready_tasks.append(task_11)
        return

    def step(self, timeStamp = 0.0):
        
        c2f_bin = self.C2F_BIN
        f2c_bin = self.F2C_BIN
        correction_bin = self.CORRECTION_BIN
        fine_bin = self.EXECUTABLE
        
        done = False
        services = self.services
        converged_slices = []
        while not done:
            if len(self.ready_tasks) > 0:
                new_task = self.ready_tasks.pop(0)
                new_task.prepare_input(c2f_bin, f2c_bin, correction_bin)
                args = new_task.get_cmd_args()
                i_suffix = new_task.suffix
                log_file_name =  'fine' + i_suffix + '.log'
                task_tag = '%04d.%04d' % (new_task.iter, new_task.slice)
                try:
                    task_id = services.launch_task(self.NPROC, self.working_dir, fine_bin, *args, \
                                          block=False, logfile=log_file_name, tag = task_tag)
                except InsufficientResourcesException:
                    self.ready_tasks.insert(0, new_task)
                except:
                    services.exception('Error launching task FINE %04d.%04d'  , 
                                       new_task.iter, new_task.slice)
                    raise
                else:
                    self.active_tasks[task_id] = new_task
                
            self.event_arrived = False
            self.events_received = []
            services.process_events()
            for event in self.events_received:
                if event['KIND'] == 'COARSE_FINISHED':
                    finished_slice = event['SLICE']
                    finished_iter = event['ITER']
                    event_coarse_outfiles = event['OUT_FILES']
                    
#                    if finished_slice == self.nt_slice:
#                        continue
                    self.update_create_task(task_iter = finished_iter,
                                            task_slice = finished_slice + 1, 
                                            dep_kind = 'COARSE', 
                                            dep_iter = finished_iter, 
                                            dep_slice = finished_slice, 
                                            dep_value = event_coarse_outfiles)
                elif event['KIND'] == 'SLICE_CONVERGED':
                    finished_slice = event['SLICE']
                    finished_iter = event['ITER']
                    converged_slices.append(finished_slice)
                    self.update_create_task(task_iter = finished_iter + 1, 
                                            task_slice = finished_slice, 
                                            dep_kind = 'CONVERGE', 
                                            dep_iter = finished_iter,
                                            dep_slice = finished_slice, 
                                            dep_value = True)
                    self.update_create_task(task_iter = finished_iter + 1, 
                                            task_slice = finished_slice + 1, 
                                            dep_kind = 'CONVERGE', 
                                            dep_iter = finished_iter,
                                            dep_slice = finished_slice, 
                                            dep_value = True)
                    for iter in range(finished_iter+1, self.max_iter + 1):
                        self.nt_slice[iter] = min(self.nt_slice[iter] + 1, self.max_slice)
                        print 'Set FINE NT_SLICE[%d] to %d' %(iter, self.nt_slice[iter])
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
                new_event = {}
                new_event['KIND'] = 'FINE_FINISHED'
                new_event['SLICE'] = task.slice
                new_event['ITER'] = task.iter
                out_files = task.get_external_outputs()
                new_event['OUT_FILES'] = out_files
                new_event['EXIT_STATUS'] = retval
                self.services.publish('PARAREAL_FINE', 'FINE_FINISHED', new_event)
                self.update_create_task(task_iter = task.iter + 1, 
                                        task_slice = task.slice + 1, 
                                        dep_kind = 'FINE', 
                                        dep_iter = task.iter,
                                        dep_slice = task.slice, 
                                        dep_value = out_files)
                
            time.sleep(0.1)
        return
    
    def finalize(self, timeStamp = 0.0):
        pass
    
    def handle_coarse_events(self, topicName, theEvent):
        self.event_arrived = True
        self.events_received.append(theEvent.getBody())
        return
    
    def handle_converge_events(self, topicName, theEvent):
        return self.handle_coarse_events(topicName, theEvent)

    def update_create_task(self, task_iter = 0, 
                                 task_slice = 0, 
                                 dep_kind = None, 
                                 dep_iter = 0, 
                                 dep_slice = 0, 
                                 dep_value= None):
        try:
            affected_task = self.task_table[task_iter, task_slice]
        except KeyError:
            affected_task = FineTask(self, task_iter, task_slice, 
                                     working_dir = self.working_dir, cmd_args = self.cmd_args)
            self.task_table[task_iter, task_slice] = affected_task
            if (task_iter == 1 or task_slice > self.nt_slice[task_iter - 1]):
                affected_task.satisfy_dependency('CONVERGE', task_iter-1, task_slice, False)            
                affected_task.satisfy_dependency('CONVERGE', task_iter-1, task_slice-1, False)            
                affected_task.satisfy_dependency('FINE', task_iter-1, task_slice-1, None)         
                affected_task.set_first_in_slice()
        affected_task.satisfy_dependency(dep_kind, dep_iter, 
                                         dep_slice, dep_value)
        
        if task_iter > self.max_iter:
            return        

        if task_slice <= self.nt_slice[task_iter]:
            if affected_task.ready_to_run():
                self.ready_tasks.append(affected_task)
                print 'Adding FINETASK %d %d to ready queue' % (task_iter, task_slice)
        return    
