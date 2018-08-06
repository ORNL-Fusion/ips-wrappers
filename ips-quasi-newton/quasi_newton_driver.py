#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for QUASI-NEWTON Driver component. This driver only implements the
#  quasi newton algorthium for optimization.
#
#-------------------------------------------------------------------------------

from component import Component
from utilities import ZipState
import os
import shutil
import json
import numpy

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver Constructor
#
#-------------------------------------------------------------------------------
class quasi_newton_driver(Component):
    def __init__(self, services, config):
        print('quasi_newton_driver: Construct')
        Component.__init__(self, services, config)
        self.model_workers = []

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        print('quasi_newton_driver: init')
    
#  Get config filenames.
        self.current_model_state = self.services.get_config_param('MODEL_INPUT')
        self.quasi_newton_config_file = self.services.get_config_param('QUASI_NEWTON_CONFIG')
        self.current_quasi_newton_state = self.services.get_config_param('CURRENT_QUASI_NEWTON_STATE')
        ips_model_config = self.services.get_config_param('MODEL_SIM_CONFIG')

#  Stage plasma state and extract all files.
        self.services.stage_plasma_state()
        with ZipState.ZipState(self.current_quasi_newton_state, 'a') as zip_ref:
            zip_ref.extractall()

#  Load the quasi-newton json file.
        with open(self.quasi_newton_config_file, 'r') as config_file:
            quasi_newton_config = json.load(config_file)
            self.signal_sigma = numpy.array(quasi_newton_config['signal_sigma'])
            self.signal_observed = numpy.array(quasi_newton_config['signal_observed'])
            self.signal_weights = numpy.sqrt(numpy.array(quasi_newton_config['signal_weights']))
            self.dchi2_tol = quasi_newton_config['dchi2_tol']

#  Set keys for the subworkflows.
            keys = {'PWD'              : self.services.get_config_param('PWD'),
                    'USER_INPUT_FILES' : self.current_model_state}

#  Copy the model state to the input file staging directory. Since all the
#  model instances start from the same state, we only need this in one place.
            if os.path.exists('model_inputs'):
                shutil.rmtree('model_inputs')
            os.mkdir('model_inputs')
            shutil.copy2(self.current_model_state, 'model_inputs')

            keywords = {}

            for i, param in enumerate(quasi_newton_config['params']):
            
                self.model_workers.append({'sim_name': None, 'init': None, 'driver': None,
                                           'result'  : '{}_result.json'.format(param['name']),
                                           'output'  : '{}_model_state.zip'.format(param['name']),
                                           'name'    : param['name'],
                                           'vrnc'    : param['vrnc'],
                                           'value'   : param['init'],
                                           'scale'   : (float(i) + 1.0)/float(len(quasi_newton_config['params']))})

                keys['SIM_NAME'] = param['name']
                keys['LOG_FILE'] = 'log.{}'.format(param['name'])
                keys['USER_OUTPUT_FILES'] = self.model_workers[i]['output']

                (self.model_workers[i]['sim_name'],
                 self.model_workers[i]['init'],
                 self.model_workers[i]['driver']) = self.services.create_sub_workflow(param['name'], ips_model_config,
                                                                                      keys, 'model_inputs')
                
                self.model_workers[i]['wait'] = self.services.call_nonblocking(self.model_workers[i]['init'],
                                                                               'init', timeStamp)

                keywords[param['name']] = param['init']

#  Run the inital convergence. Only one model is needed to run but the staging
#  expects values from each sub_work_flow so we need to launch them all.
            for worker in self.model_workers:
                self.services.wait_call(worker['wait'], True)
                worker['wait'] = self.services.call_nonblocking(worker['driver'], 'init', timeStamp,
                                                                **keywords)

            for worker in self.model_workers:
                self.services.wait_call(worker['wait'], True)
                worker['wait'] = self.services.call_nonblocking(worker['driver'], 'step', timeStamp,
                                                                result_file=worker['result'])

            for worker in self.model_workers:
                self.services.wait_call(worker['wait'], True)
                            
            self.services.stage_subflow_output_files()
                
            with ZipState.ZipState(self.model_workers[0]['output'], 'a') as zip_ref:
                zip_ref.extract(self.model_workers[0]['result'])
                with open(self.model_workers[0]['result'], 'r') as result_file:
                    self.signal_model = numpy.array(json.load(result_file)['signal_model'])
                
#  FIXME: DEBUGING REMOVE
                os.rename(self.model_workers[0]['result'], 'saved.json')

#  The initial model state may have changed reset it from the one of the
#  workers.
            os.rename(self.model_workers[0]['output'], self.current_model_state)

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver step method. This runs the vmec component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        print('quasi_newton_driver: step')

#  Compute chi^2. Set the inital change in chi^2 higher than the tolarance to
#  ensure at least one iteration of the while loop is performed.
        e = self.signal_weights*((self.signal_model - self.signal_observed)/self.signal_sigma)
        chi2 = numpy.sum(e*e)
        dchi2 = self.dchi2_tol + 1.0

        jacobian = numpy.empty((len(self.model_workers), e.size), dtype=float)
#        hessian = numpy.empty((len(self.model_workers), len(self.model_workers)), dtype=float)
#        gradient = numpy.empty(len(self.model_workers), dtype=float)
        etry = numpy.empty((len(self.model_workers), e.size), dtype=float)
        chi2try = numpy.empty(len(self.model_workers), dtype=float)

#  Perform a quasi-newton minimization.
        while dchi2 > self.dchi2_tol:
            timeStamp += 1.0
    
#  Set the model to a known state.
            shutil.copy2(self.current_model_state, 'model_inputs')
            for worker in self.model_workers:
                worker['wait'] = self.services.call_nonblocking(worker['init'],
                                                                'init', timeStamp)

#  Perturb the parameters.
            for worker in self.model_workers:
                self.services.wait_call(worker['wait'], True)
                keywords = {worker['name'] : worker['value'] + worker['vrnc']}
                print(keywords[worker['name']])
                worker['wait'] = self.services.call_nonblocking(worker['driver'],
                                                                'init', timeStamp,
                                                                **keywords)
            
#  Recompute the model.
            for worker in self.model_workers:
                self.services.wait_call(worker['wait'], True)
                worker['wait'] = self.services.call_nonblocking(worker['driver'], 'step', timeStamp,
                                                                result_file=worker['result'])

#  Collect the results and reset the model.
            for worker in self.model_workers:
                self.services.wait_call(worker['wait'], True)
                worker['wait'] = self.services.call_nonblocking(worker['init'],
                                                                'init', timeStamp)
                        
            self.services.stage_subflow_output_files()

#  Compute the normalized jacobian A.
#
#    A_ij = d e_i/d a_j                                                      (1)
#
#  Where e is the error vector.
#
#    e_i = W_i*((S_i - M_i)/sigma_i)^2                                       (2)
#
#  Note due to the what the memory is laid out the Jacobian is transposed.
            for i, worker in enumerate(self.model_workers):
                with ZipState.ZipState(worker['output'], 'a') as zip_ref:
                    zip_ref.extract(worker['result'])
                    with open(worker['result'], 'r') as result_file:
                        self.signal_model = numpy.array(json.load(result_file)['signal_model'])

                ep = self.signal_weights*((self.signal_model - self.signal_observed)/self.signal_sigma)
                jacobian[i] = e - ep
#  Approximate the hessian.
#
#    alpha = A^T * A
#
#  The gradient
#
#    beta = A^T * e
#
#            hessian = numpy.dot(jacobian, numpy.transpose(jacobian))
#            gradient = numpy.dot(jacobian, e)
#            delta_a = gradient*numpy.dot(gradient, gradient)
#            delta_a /= numpy.dot(gradient, numpy.dot(hessian, gradient))
            j_svd_u, j_svd_w, j_svd_vt = numpy.linalg.svd(numpy.transpose(jacobian))

            print(numpy.transpose(jacobian))
            print('------------------------------------------------------------')
            print(j_svd_u)
            print('------------------------------------------------------------')
            print(j_svd_w)
            print('------------------------------------------------------------')
            print(j_svd_vt)

            delta_a = numpy.empty((len(j_svd_w), len(self.model_workers)), dtype=float)
            
            temp_work1 = numpy.where(j_svd_w > 0, 1.0/j_svd_w, 0)
            temp_work2 = numpy.dot(numpy.transpose(j_svd_u), e)
            temp_work3 = numpy.zeros(len(self.model_workers))
            for i in range(len(j_svd_w)):
                temp_work3[i] = temp_work1[i]*temp_work2[i]
                delta_a[i] = numpy.dot(numpy.transpose(j_svd_vt), temp_work3)

            new_values = []
            for i, worker in enumerate(self.model_workers):
                new_values.append(delta_a[i]*worker['vrnc'])

            for worker in self.model_workers:
                keywords = {}
                for new_value in new_values:
                    keywords[worker['name']] = worker['value'] + new_value*(worker['scale']*50.0 + 50.0)
                self.services.wait_call(worker['wait'], True)
                worker['wait'] = self.services.call_nonblocking(worker['driver'],
                                                                'init', timeStamp,
                                                                **keywords)

            for worker in self.model_workers:
                self.services.wait_call(worker['wait'], True)
                worker['wait'] = self.services.call_nonblocking(worker['driver'], 'step', timeStamp,
                                                                result_file=worker['result'])

#  Collect the results and reset the model.
            for worker in self.model_workers:
                self.services.wait_call(worker['wait'], True)
            self.services.stage_subflow_output_files()

            for i, worker in enumerate(self.model_workers):
                with ZipState.ZipState(worker['output'], 'a') as zip_ref:
                    zip_ref.extract(worker['result'])
                    with open(worker['result'], 'r') as result_file:
                        self.signal_model = numpy.array(json.load(result_file)['signal_model'])
            
                etry[i] = self.signal_weights*((self.signal_observed - self.signal_model)/self.signal_sigma)
                chi2try[i] = numpy.dot(etry[i], etry[i])

            min_i = numpy.argmin(chi2try)
            dchi2 = chi2 - chi2try[min_i]
            chi2 = chi2try[min_i]

            for i, worker in enumerate(self.model_workers):
                worker['value'] += new_values[i]*(self.model_workers[min_i]['scale']*50.0 + 50.0)

            os.rename(self.model_workers[min_i]['output'], self.current_model_state)
            print('--------------------------------------------------------------------------------')
            print('--------------------------------------------------------------------------------')
            print('--------------------------------------------------------------------------------')
            print(timeStamp, min_i, chi2, dchi2)
            print('--------------------------------------------------------------------------------')
            print('--------------------------------------------------------------------------------')
            print('--------------------------------------------------------------------------------')

#  FIXME: Force the loop to end.
            dchi2 = 0

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver finalize method. This cleans up afterwards.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        print('quasi_newton_driver: finalize')

        for worker in self.model_workers:
            worker['wait'] = [
                              self.services.call_nonblocking(worker['init'], 'finalize', timeStamp),
                              self.services.call_nonblocking(worker['driver'], 'finalize', timeStamp)
                             ]

        for worker in self.model_workers:
            self.services.wait_call_list(worker['wait'], True)

