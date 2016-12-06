#! /usr/bin/env python

from  component import Component
import os
import glob
import shutil
import subprocess
from math import sqrt
from ipsExceptions import InsufficientResourcesException
import time

class ConvergeTask(object):
    def __init__(self, comp, iteration, slice, in_files = None, out_files = None, working_dir = None):
        self.comp = comp
        self.services = comp.services
        self.iter = iteration
        self.slice = slice
        if in_files:
            self.in_files = [f for f in in_files]
        self.out_files = {}
        self.dependencies = {}   # To be filled when dependencies are satisfied 
        self.status = None
        self.task_id = None
        self.submit_args = []
        self.suffix = '%04d.%04d' % (self.iter, self.slice)
        self.working_dir = working_dir
        
    def satisfy_dependency(self, kind, iteration, slice, value):
        if kind  != 'FINE':
            raise Exception('Error: CONVERGE Task %d:%d received wrong dependency type %s' % 
                            (self.iter, self.slice, kind))
        self.dependencies[kind, iteration, slice] = value
        print 'CONVERGETASK %d %d : Satisfied %s %d %d with %s' % \
            (self.iter, self.slice, str(kind), int(iteration), int(slice), str(value))
        return
    
    def prepare_input(self):
        fine_in_cur_iter = self.dependencies['FINE', self.iter, self.slice]
        fine_in_prev_iter = self.dependencies['FINE', self.iter-1, self.slice]
        if (fine_in_cur_iter and fine_in_prev_iter):
            prev_iter_values = fine_in_prev_iter['FINE_OUT_NC']
            cur_iter_values = fine_in_cur_iter['FINE_OUT_NC']
            converge_out = 'converge_out.' + self.suffix
#            num_lines = len(open(cur_iter_values).readlines())
            self.submit_args = [prev_iter_values, cur_iter_values, converge_out]
            self.out_files = {'CONVERGE_OUT' : os.path.join(self.working_dir, converge_out)}
        else:
            self.submit_args = None
        return
    
    def process_output(self):
        pass
    
    def get_cmd_args(self):
        return self.submit_args
    
    def get_cmd(self):
        return 
    
    def ready_to_run(self):
        if self.iter > self.comp.max_iter:
            return False
        if len(self.dependencies) == 2:
            print 'CONVERGETASK %d %d ready to run' % (self.iter, self.slice)
            return True
        print 'CONVERGETASK %d %d NOT ready to run - dep = %s' % \
                (self.iter, self.slice, str(self.dependencies))
        return False    
        
    def get_external_outputs(self):
        return self.out_files
    
class Converge(Component):
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
                    
        self.first_non_converge = {}


    def init(self, timeStamp = 0.0):
        self.working_dir = self.services.get_working_dir()
        self.services.subscribe('PARAREAL_FINE', self.handle_fine_events)
        iter = 1
        for slice in range(1,self.max_slice+1):
            new_task = self.update_create_task(iter, slice, 'FINE', iter-1, slice, None)
        return

    def step(self, timeStamp = 0.0):
        
        converge_bin = self.CONV_BIN
        done = False
        services = self.services
        curdir = self.services.get_working_dir()
        converged_slices = {}
        last_true_converge = {}
        last_true_converge[1] = 1
        error = {}
        num_converged_slices = 0
        num_converge_tests = {}
        for i_iter in range(1, self.max_iter + 1):
            num_converge_tests[i_iter] = 0
        
        while not done:
            if len(self.ready_tasks) > 0:
                new_task = self.ready_tasks.pop(0)
                new_task.prepare_input()
                args = new_task.get_cmd_args()
                task_iter = new_task.iter
                task_slice = new_task.slice
                slice_converged = False
                num_converge_tests[task_iter] += 1
                if not args:
                    conv_str = '0'
                else:
                    cmd_list = [converge_bin] + args + [self.TOL]
                    print 'Converge Command ', ' '.join(cmd_list)
                    subprocess.call(cmd_list)
                    conv_file = new_task.get_external_outputs()['CONVERGE_OUT'] 
                    conv_str = open(conv_file, 'r').readline()
                    print 'Convergence Str =',  conv_str
                    (converge, actual_tolerance) =conv_str.strip().split()
                    if converge == '1':
                        slice_converged = True
                    err = abs(float(actual_tolerance))
                    err_str = '%10.3e' % err
                    error[task_iter, task_slice] = err
                    print 'convergence string = ', converge, actual_tolerance, 'iteration = ', new_task.iter
                    self.services.send_portal_event('converge_out', str(new_task.iter) + ' ' + 
                                                    str(new_task.slice) +  ' ' +  err_str)
                if (task_iter > 1 ):
                    if (task_slice == self.first_non_converge[task_iter-1]):
                        slice_converged = True
                else:
                    if (task_slice == 1):
                        slice_converged = True
                    
                print 'diag converge ', task_iter, task_slice, self.first_non_converge, slice_converged
                print 'diag converge ', task_iter, task_slice, num_converged_slices
                if not slice_converged:
                    if task_slice < self.first_non_converge[task_iter]:
                        new_event = {}
                        new_event['KIND'] = 'SLICE_NOT_CONVERGED'
                        new_event['ITER'] = task_iter
                        new_event['OUT_FILES'] = None
                        new_event['EXIT_STATUS'] = None
                        print 'diag1a converge ', task_iter,  task_slice,  self.first_non_converge[task_iter]
                        for i_slice in range(task_slice, self.first_non_converge[task_iter] + 1):
                            new_event['SLICE'] = i_slice
                            self.services.publish('PARAREAL_CONVERGE', 'CONVERGE_FINISHED', new_event)
                        self.first_non_converge[task_iter] = min (task_slice, self.first_non_converge[task_iter])
                else:
                    print 'made it into the else loop ', task_iter, task_slice, converged_slices, last_true_converge
                    
                    try:
                        converged_list = converged_slices[task_iter]
                    except KeyError:
                        print 'in the KeyError except', task_iter, task_slice
                        converged_slices[task_iter] = []
                        converged_list = converged_slices[task_iter]
                        if (task_iter > 1):
                            last_true_converge[task_iter] = self.first_non_converge[task_iter - 1] - 1
                        else:
                            last_true_converge[task_iter] = 0
                    
                    #print 'diag1 converge ', task_iter, task_slice, self.first_non_converge, last_true_converge, converged_slices
                    if task_slice <= self.first_non_converge[task_iter]:
                        converged_list.append(task_slice)
                        converged_slices[task_iter] = sorted(converged_list)
                        last_published_converge = last_true_converge[task_iter]
                        new_event = {}
                        new_event['KIND'] = 'SLICE_CONVERGED'
                        new_event['ITER'] = task_iter
                        new_event['OUT_FILES'] = None
                        new_event['EXIT_STATUS'] = None
                        #print 'before for loop ', task_iter, task_slice, range( last_published_converge + 1, max(converged_list) + 1)
                        for candidate_slice in range( last_published_converge + 1, 
                                                     max(converged_list) + 1):
                            print 'diag2 converge', candidate_slice, converged_list
                            if candidate_slice in converged_list:
                                new_event['SLICE'] = candidate_slice
                                self.services.publish('PARAREAL_CONVERGE', 'CONVERGE_FINISHED', new_event)
                                last_true_converge[task_iter] = candidate_slice
                                num_converged_slices += 1
                                for iter in range(task_iter + 1, self.max_iter):
                                    self.nt_slice[iter] = min(self.nt_slice[iter]+ 1, self.max_slice)
                                #print 'diag converge ', task_iter, task_slice, num_converged_slices
                            else:
                                break
                        if task_iter > 1 and last_true_converge[task_iter] == self.nt_slice[task_iter-1]:
                            self.first_non_converge[task_iter] = last_true_converge[task_iter] + 1

                            new_event = {}
                            new_event['KIND'] = 'SLICE_NOT_CONVERGED'
                            new_event['ITER'] = task_iter
                            new_event['OUT_FILES'] = None
                            new_event['EXIT_STATUS'] = None
                            for slice in range(last_true_converge[task_iter]+ 1, self.nt_slice[iter]+1):
                                new_event['SLICE'] = slice
                                self.services.publish('PARAREAL_CONVERGE', 'CONVERGE_FINISHED', new_event)
                            
            #print 'diag converge for exit test ', self.nt_slice, num_converged_slices
            if (num_converged_slices == self.max_slice):
                done = True
                new_event = {}
                new_event['KIND'] = 'ALL_CONVERGED'
                self.services.publish('PARAREAL_CONVERGE', 'CONVERGE_DONE', new_event)
                break
            if task_iter == self.max_iter : 
                print 'max iter test ', task_iter, self.max_iter, num_converge_tests[self.max_iter], \
                      last_true_converge[self.max_iter - 1]
                if (num_converge_tests[self.max_iter] == self.max_slice - \
                                       last_true_converge[self.max_iter - 1]):
                    done = True
                    new_event = {}
                    new_event['KIND'] = 'ALL_ABORT'
                    self.services.publish('PARAREAL_CONVERGE', 'CONVERGE_DONE', new_event)
                    break
                            
            self.event_arrived = False
            self.events_received = []
            services.process_events()
            for event in self.events_received:
                if event['KIND'] == 'FINE_FINISHED':
                    finished_slice = event['SLICE']
                    finished_iter = event['ITER']
                    event_outfiles = event['OUT_FILES']
                    
                    affected_slice = finished_slice
                    for affected_iter in (finished_iter, finished_iter + 1):
                        self.update_create_task(task_iter = affected_iter,
                                                task_slice = affected_slice, 
                                                dep_kind = 'FINE', 
                                                dep_iter = finished_iter, 
                                                dep_slice = finished_slice, 
                                                dep_value = event_outfiles)
            if done:
                break
            time.sleep(0.1)
        return
    
    def finalize(self, timeStamp = 0.0):
        pass
    
    def handle_fine_events(self, topicName, theEvent):
        self.event_arrived = True
        self.events_received.append(theEvent.getBody())
        return
    
    def update_create_task(self, task_iter = 0, 
                                 task_slice = 0, 
                                 dep_kind = None, 
                                 dep_iter = 0, 
                                 dep_slice = 0, 
                                 dep_value= None):
        try:
            affected_task = self.task_table[task_iter, task_slice]
        except KeyError:
            affected_task = ConvergeTask(self, task_iter, task_slice, working_dir = self.working_dir)
            self.task_table[task_iter, task_slice] = affected_task
            if task_iter not in self.first_non_converge.keys():
                self.first_non_converge[task_iter] = self.max_slice + 1
            if task_iter == 1:
                affected_task.satisfy_dependency('FINE', task_iter-1, task_slice, None)            
                affected_task.satisfy_dependency('FINE', task_iter, task_slice, None)
        affected_task.satisfy_dependency(dep_kind, dep_iter, 
                                         dep_slice, dep_value)
        if affected_task.ready_to_run():
            self.ready_tasks.append(affected_task)
