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
import math

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
            self.signal_sigma = numpy.absolute(numpy.array(quasi_newton_config['signal_sigma']))
            self.signal_observed = numpy.array(quasi_newton_config['signal_observed'])
            self.signal_weights = numpy.sqrt(numpy.array(quasi_newton_config['signal_weights']))
            self.dchi2_tol = quasi_newton_config['dchi2_tol']
        
#  Singular value step controls.
#  Cutoff value for relative singular values.
            if 'cut_svd' in quasi_newton_config:
                self.cut_svd = quasi_newton_config['cut_svd']
            else:
                self.cut_svd = 0.0
#  Cutoff value for expected step efficiency.
            if 'cut_eff' in quasi_newton_config:
                self.cut_eff = quasi_newton_config['cut_eff']
            else:
                self.cut_eff = 0.0
#  Cutoff value for expected marginal step efficiency.
            if 'cut_marg_eff' in quasi_newton_config:
                self.cut_marg_eff = quasi_newton_config['cut_marg_eff']
            else:
                self.cut_marg_eff = 0.0
#  Cutoff value for expected step size.
            if 'cut_delta_a' in quasi_newton_config:
                self.cut_delta_a = quasi_newton_config['cut_delta_a']
            else:
                self.cut_delta_a = 0.0
#  Cutoff value for expected change in g^2.
            if 'cut_dg2' in quasi_newton_config:
                self.cut_dg2 = quasi_newton_config['cut_dg2']
            else:
                self.cut_dg2 = 0.0
            
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
        self.e = self.signal_weights*((self.signal_model - self.signal_observed)/self.signal_sigma)
        chi2 = numpy.dot(self.e, self.e)
        dchi2 = self.dchi2_tol + 1.0

        self.jacobian = numpy.empty((len(self.model_workers), self.e.size), dtype=float)
        self.hessian = numpy.empty((len(self.model_workers), len(self.model_workers)), dtype=float)
        self.gradient = numpy.empty(len(self.model_workers), dtype=float)
        self.delta_a = numpy.empty((min(len(self.e), len(self.model_workers)) + 1, len(self.model_workers)), dtype=float)

        etry = numpy.empty((len(self.model_workers), self.e.size), dtype=float)
        chi2try = numpy.empty(len(self.model_workers), dtype=float)

#  Perform a quasi-newton minimization.
        while dchi2 > self.dchi2_tol:
            timeStamp += 1.0
            
            self.eval_jacobian(timeStamp)
            num_sv = self.get_k_svd()

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

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver eval_jacobian method. This perturbs each parameter
#  independently in parallel.
#
#-------------------------------------------------------------------------------
    def eval_jacobian(self, timeStamp=0.0):
        print('quasi_newton_driver: eval_jacobian')

#  Set the model to a known state.
        shutil.copy2(self.current_model_state, 'model_inputs')
        for worker in self.model_workers:
            worker['wait'] = self.services.call_nonblocking(worker['init'],
                                                            'init', timeStamp)
                
#  Perturb the parameters.
        for worker in self.model_workers:
            self.services.wait_call(worker['wait'], True)
            keywords = {worker['name'] : worker['value'] + worker['vrnc']}
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
                                                                                            
            self.jacobian[i] = self.e - self.signal_weights*((self.signal_model - self.signal_observed)/self.signal_sigma)

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver get_k_svd method. Performs a singular value decomposition
#  of the jacobian. Based on various cuttoffs, the maximum number of sigular
#  values retained is determined.
#
#-------------------------------------------------------------------------------
    def get_k_svd(self):
        print('quasi_newton_driver: get_k_svd')

#  Approximate the hessian.
#
#    alpha = A^T * A                                                         (1)
#
#  The gradient
#
#    beta = A^T * e                                                          (2)
#
        self.hessian = numpy.matmul(self.jacobian, numpy.transpose(self.jacobian))
        self.gradient = numpy.dot(self.jacobian, self.e)
        
#  Compute the steepest descent setp size in normalized parameter space.
#
#    da = beta * beta * beta / (beta * alpha * beta)                         (3)
        self.delta_a[0,:] = self.gradient*numpy.dot(self.gradient, self.gradient)/numpy.dot(self.gradient, numpy.matmul(self.hessian, self.gradient))

#  Singular value decomposition of the Jacobian.
        j_svd_u, j_svd_w, j_svd_vt = numpy.linalg.svd(numpy.transpose(self.jacobian))
        
#  Define Inverse singular values.
        temp_work1 = numpy.where(j_svd_w > 0, 1.0/j_svd_w, 0)
#    U^T * e                                                                 (4)
        temp_work2 = numpy.dot(numpy.transpose(j_svd_u), self.e)

#  Perform pseudo-inversion for successive numbers of singular values retained.
        temp_work3 = numpy.zeros(len(self.model_workers))
        for i in range(0, len(j_svd_w)):
            temp_work3[i] = temp_work1[i]*temp_work2[i]
            self.delta_a[i + 1,:] = numpy.matmul(numpy.transpose(j_svd_vt), temp_work3)

#  Estimate the expected changes in g^2. Equation 22 in Hanson et. al.
#  doi: 10.1088/0029-5515/49/7/075031
        exp_dg2 = numpy.empty(len(j_svd_w) + 1, dtype=float)
        delta_a_len = numpy.empty(len(j_svd_w) + 1, dtype=float)
        exp_eff = numpy.empty(len(j_svd_w) + 1, dtype=float)
        for i in range(0, len(exp_dg2)):
            exp_dg2[i] = self.get_exp_dg2(self.delta_a[i,:])
            delta_a_len[i] = math.sqrt(numpy.dot(self.delta_a[i,:], self.delta_a[i,:]))
            
#  Equation 23 in Hanson et. al. doi: 10.1088/0029-5515/49/7/075031
            exp_eff[i] = abs(exp_dg2[i])/delta_a_len[i]


#  Although marginal efficiencies are not strictly defined for the Steepest
#  Descent (index 0) and just one singular value cases, for convenience define
#  them here.
        marg_exp_eff = numpy.empty(len(j_svd_w) + 1, dtype=float)
        marg_exp_eff[0:1] = exp_eff[0:1]
        for i in range(2, len(exp_dg2)):
            d_len = max(delta_a_len[i] - delta_a_len[i - 1], 1.0e-10)
            marg_exp_eff[i] = abs(exp_dg2[i]) - abs(exp_dg2[i - 1])/d_len
    
#  Check for cutoffs.
        largest_w = 1.0
        if j_svd_w[0] > 0:
            largest_w = j_svd_w[0]

#  Find the largest number of singular values to use.
        for i in range(len(j_svd_w), 1, -1):
            meets_cut = (j_svd_w[i - 1]/largest_w >= self.cut_svd)
            meets_cut = meets_cut and (exp_eff[i] >= self.cut_eff)
            meets_cut = meets_cut and (marg_exp_eff[i] >= self.cut_marg_eff)
            meets_cut = meets_cut and (delta_a_len[i] >= self.cut_delta_a)
            meets_cut = meets_cut and (numpy.abs(exp_dg2[i]) >= self.cut_dg2)
            if meets_cut:
                return i

#  All the selected criteria failed use no singular values.
        return 0
        
#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver get_exp_dg2 method. Obtains the expected change in g^2.
#
#-------------------------------------------------------------------------------
    def get_exp_dg2(self, delta):
        print('quasi_newton_driver: get_exp_dg2')
        
#  Linear part.
#
#    2e * A * da                                                             (1)
        exp_dg2_lin = 2.0*numpy.dot(self.e, numpy.matmul(numpy.transpose(self.jacobian), delta))

#  Qaudratic part.
#
#    da * A^T * A * da = da * alpha * da                                     (2)
#
#  Alpha is the hessian matrix. Equation 14 in Hanson et. al.
#  doi: 10.1088/0029-5515/49/7/075031
        exp_dg2_quad = numpy.dot(delta, numpy.matmul(self.hessian, delta))

        return exp_dg2_lin - exp_dg2_quad
