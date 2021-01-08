#! /usr/bin/env python

#-------------------------------------------------------------------------------
#
#  IPS driver for QUASI-NEWTON Driver component. This driver only implements the
#  quasi newton algorthium for optimization.
#
#-------------------------------------------------------------------------------

from ipsframework import Component
from utilities import ZipState
from utilities import ScreenWriter
import os
import shutil
import json
import numpy
import math
from scipy import optimize
import time

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver Constructor
#
#-------------------------------------------------------------------------------
class quasi_newton_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        self.model_workers = []
        self.history = []

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver init method.
#
#-------------------------------------------------------------------------------
    def init(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_driver: init')
    
#  Get config filenames.
        self.current_model_state = self.services.get_config_param('MODEL_INPUT')
        self.quasi_newton_config_file = self.services.get_config_param('QUASI_NEWTON_CONFIG')
        self.current_quasi_newton_state = self.services.get_config_param('CURRENT_QUASI_NEWTON_STATE')
        ips_model_config = self.services.get_config_param('MODEL_SIM_CONFIG')

#  Stage state and extract all files.
        self.services.stage_state()
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

#  Minimization controls
            if 'max_step' in quasi_newton_config:
                self.max_step = quasi_newton_config['max_step']
            else:
                self.max_step = 100.0

#  Maximum reconstruction steps
            if 'max_recon_steps' in quasi_newton_config:
                self.max_recon_steps = quasi_newton_config['max_recon_steps']
            else:
                self.max_recon_steps = 20

#  Maximum number of trys a step can take to reduce g^2
            if 'max_step_try' in quasi_newton_config:
                self.max_step_try = quasi_newton_config['max_step_try']
            else:
                self.max_step_try = 10

#  Set keys for the subworkflows.
            keys = {'PWD'              : self.services.get_config_param('PWD'),
                    'USER_INPUT_FILES' : self.current_model_state,
                    'OUTPUT_LEVEL'     : self.services.get_config_param('OUTPUT_LEVEL')}

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
                    self.e = self.get_e(result_file)

#  The initial model state may have changed reset it from the one of the
#  workers.
            os.rename(self.model_workers[0]['output'], self.current_model_state)

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver step method. This runs the vmec component.
#
#-------------------------------------------------------------------------------
    def step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_driver: step')

#  Compute chi^2. Set the inital change in chi^2 higher than the tolarance to
#  ensure at least one iteration of the while loop is performed.
        self.chi2 = numpy.dot(self.e, self.e)
        dchi2 = self.dchi2_tol + 1.0
        ScreenWriter.screen_output(self, 'quiet',   '-------------------------------------------------------------------------------------------------')
        ScreenWriter.screen_output(self, 'quiet',   'Step {:>4.0f} : chi^2 = {:12.5e}'.format(timeStamp, self.chi2))
        ScreenWriter.screen_output(self, 'verbose', '-------------------------------------------------------------------------------------------------')

        self.jacobian = numpy.empty((len(self.model_workers), self.e.size), dtype=float)
        self.hessian = numpy.empty((len(self.model_workers), len(self.model_workers)), dtype=float)
        self.gradient = numpy.empty(len(self.model_workers), dtype=float)
        self.delta_a = numpy.empty((min(len(self.e), len(self.model_workers)) + 1, len(self.model_workers)), dtype=float)

        self.history.append({'Time' : timeStamp, 'Chi2' : self.chi2})

#  Perform a quasi-newton minimization.
        while dchi2 > self.dchi2_tol and timeStamp < self.max_recon_steps:
            timeStamp += 1.0

            self.eval_jacobian(timeStamp)
            if self.try_step(timeStamp):
                new_chi2 = numpy.dot(self.e, self.e)
                dchi2 = self.chi2 - new_chi2
                self.chi2 = new_chi2
                self.history.append({'Time' : timeStamp, 'Chi2' : self.chi2, 'DChi2' : dchi2, 'Num SV' : self.k_use, 'Size' : self.norm_len})
                ScreenWriter.screen_output(self, 'verbose', '-------------------------------------------------------------------------------------------------')
                ScreenWriter.screen_output(self, 'quiet',   'Step {:>4.0f} : chi^2 = {:12.5e} : dchi^2 = {:12.5e} : Num SV = {:4} : Norm Size {:12.5e}'.format(timeStamp, self.chi2, dchi2, self.k_use, self.norm_len))
                ScreenWriter.screen_output(self, 'verbose', '-------------------------------------------------------------------------------------------------')
            else:
                break

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver finalize method. This cleans up afterwards.
#
#-------------------------------------------------------------------------------
    def finalize(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_driver: finalize')

        for worker in self.model_workers:
            worker['wait'] = [
                              self.services.call_nonblocking(worker['init'], 'finalize', timeStamp),
                              self.services.call_nonblocking(worker['driver'], 'finalize', timeStamp)
                             ]

        for worker in self.model_workers:
            self.services.wait_call_list(worker['wait'], True)

        ScreenWriter.screen_output(self, 'quiet', '-------------------------------------------------------------------------------------------------')
        ScreenWriter.screen_output(self, 'quiet', '{:<4} : {:<12} : {:<12} : {:<6} : {:<12}'.format('Step', 'chi^2', 'dchi^2', 'Num SV', 'Norm Size'))
        ScreenWriter.screen_output(self, 'quiet', '')
        for step in self.history:
            if "DChi2" in step:
                ScreenWriter.screen_output(self, 'quiet', '{Time:>4.0f} : {Chi2:>12.5e} : {DChi2:>12.5e} : {Num SV:>6} : {Size:>12.5e}'.format(**step))
            else:
                ScreenWriter.screen_output(self, 'quiet', '{Time:>4.0f} : {Chi2:>12.5e}'.format(**step))
        ScreenWriter.screen_output(self, 'quiet', '-------------------------------------------------------------------------------------------------')

        correlation_matrix = numpy.linalg.inv(self.hessian)
        ScreenWriter.screen_output(self, 'quiet', '{:<50} {:<12} {:<12}'.format('Parameter', 'Value', 'Sigma'))
        ScreenWriter.screen_output(self, 'quiet', '')
        for i, worker in enumerate(self.model_workers):
            worker['sigma'] = math.sqrt(correlation_matrix[i,i]*worker['vrnc']**2.0)
            ScreenWriter.screen_output(self, 'quiet', '{name:<50} {value:>12.5e} {sigma:>12.5e}'.format(**worker))
        ScreenWriter.screen_output(self, 'quiet', '-------------------------------------------------------------------------------------------------')

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver eval_jacobian method. This perturbs each parameter
#  independently in parallel.
#
#-------------------------------------------------------------------------------
    def eval_jacobian(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_driver: eval_jacobian')

#  Set the model to a known state.
        shutil.copy2(self.current_model_state, 'model_inputs')
        for worker in self.model_workers:
            worker['wait'] = self.services.call_nonblocking(worker['init'],
                                                            'init', timeStamp)
                
#  Perturb the parameters.
        for worker in self.model_workers:
            keywords = {worker['name'] : worker['value'] + worker['vrnc']}
            self.services.wait_call(worker['wait'], True)
            worker['wait'] = self.services.call_nonblocking(worker['driver'],
                                                            'init', timeStamp,
                                                            **keywords)
                                                            
#  Recompute the model.
        for worker in self.model_workers:
            self.services.wait_call(worker['wait'], True)
            worker['wait'] = self.services.call_nonblocking(worker['driver'], 'step', timeStamp,
                                                            result_file=worker['result'])
                                                                    
#  Collect the results.
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
                    self.jacobian[i] = self.e - self.get_e(result_file)

        with open('jacobian.log', 'a') as jacobian_ref:
            jacobian_ref.write('Jacobian step {}\n'.format(timeStamp));
            for j in range(len(self.e)):
                self.jacobian[:,j].tofile(jacobian_ref, sep=',', format='%12.5e');
                jacobian_ref.write('\n')
            jacobian_ref.write('\n')

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver get_e method. This computes the non squared error.
#
#-------------------------------------------------------------------------------
    def get_e(self, result_file):
        return self.signal_weights*((numpy.array(json.load(result_file)['signal_model']) - self.signal_observed)/self.signal_sigma)

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver try_step method. Trys different step sizes to move to the
#  minimum in parameter space. Since there are potentially up to the number of
#  parameter parallel instances, try that many different step sizes.
#
#-------------------------------------------------------------------------------
    def try_step(self, timeStamp=0.0):
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_driver: try_step')

        self.k_use = self.get_k_svd()

#  Try different Levenberg-Marquardt step sizes.
        new_max = min(self.delta_a_len[self.k_use], self.max_step)
        step_use = numpy.empty(len(self.model_workers), dtype=float)
        delta_try = numpy.empty((len(self.model_workers), len(self.model_workers)), dtype=float)
        e_try = numpy.empty((len(self.model_workers), len(self.signal_observed)), dtype=float)
        chi2try = numpy.empty(len(self.model_workers), dtype=float)
        
        num_trys = 0
        
        while num_trys < self.max_step_try:
            num_trys += 1
            
            for i, worker in enumerate(self.model_workers):
                step_use[i] = new_max - i*new_max/(2.0*len(self.model_workers))
                delta_try[i] = self.lm_step(step_use[i])

#  Set new parameters.
                keywords = {}
                for j, worker2 in enumerate(self.model_workers):
                    keywords[worker2['name']] = worker2['value'] + delta_try[i,j]*worker2['vrnc']
                self.services.wait_call(worker['wait'], True)
                worker['wait'] = self.services.call_nonblocking(worker['driver'],
                                                                'init', timeStamp,
                                                                **keywords)

#  Recompute the model.
            for worker in self.model_workers:
                self.services.wait_call(worker['wait'], True)
                worker['wait'] = self.services.call_nonblocking(worker['driver'], 'step', timeStamp,
                                                                result_file=worker['result'])

#  Collect the results.
            for worker in self.model_workers:
                self.services.wait_call(worker['wait'], True)

            self.services.stage_subflow_output_files()

#  Compute chi^2 for each attempted step. And keep the largest.
            with open('chi.log', 'a') as chi_ref:
                chi_ref.write('Chi step {}\n'.format(timeStamp));
                for i, worker in enumerate(self.model_workers):
                    with ZipState.ZipState(worker['output'], 'a') as zip_ref:
                        zip_ref.extract(worker['result'])
                    with open(worker['result'], 'r') as result_file:
                        e_try[i] = self.get_e(result_file)
                        chi2try[i] = numpy.dot(e_try[i], e_try[i])

                    chi_ref.write('chi2 = {} : '.format(chi2try[i]))
                    e_try[i].tofile(chi_ref, sep=',', format='%12.5e')
                    chi_ref.write('\n')
                chi_ref.write('\n')

            i_min = numpy.argmin(chi2try)
            if chi2try[i_min] <= self.chi2:
#  Chi^2 decreased. Set the best case to the current model.
                os.rename(self.model_workers[i_min]['output'], self.current_model_state)
                shutil.copy2(self.current_model_state, 'model_inputs')
                self.e = e_try[i_min]
                
#  Set the new parameter values.
                current_values = {}
                for i, worker in enumerate(self.model_workers):
                    worker['value'] += delta_try[i_min,i]*worker['vrnc']
                    current_values[worker['name']] = worker['value']
                
#  Dump current values to a json file. This can be used for restarting a
#  reconstruction.
                with open('current_values.json', 'w') as current_values_file:
                    json.dump(current_values, current_values_file)

                self.norm_len = numpy.sqrt(numpy.dot(delta_try[i_min], delta_try[i_min]))

                return True
            else:
#  Cut the step size in half and reset the model.
                new_max /= 2.0
                for worker in self.model_workers:
                    worker['wait'] = self.services.call_nonblocking(worker['init'],
                                                                    'init', timeStamp)

        return False

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver get_k_svd method. Performs a singular value decomposition
#  of the jacobian. Based on various cuttoffs, the maximum number of sigular
#  values retained is determined.
#
#-------------------------------------------------------------------------------
    def get_k_svd(self):
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_driver: get_k_svd')
    
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
        self.j_svd_u, self.j_svd_w, self.j_svd_vt = numpy.linalg.svd(numpy.transpose(self.jacobian))
        
#  Define Inverse singular values.
        temp_work1 = numpy.where(self.j_svd_w > 0, 1.0/self.j_svd_w, 0)
#    U^T * e                                                                 (4)
        temp_work2 = numpy.dot(numpy.transpose(self.j_svd_u), self.e)
        
#  Perform pseudo-inversion for successive numbers of singular values retained.
        temp_work3 = numpy.zeros(len(self.model_workers))
        for i in range(0, len(self.j_svd_w)):
            temp_work3[i] = temp_work1[i]*temp_work2[i]
            self.delta_a[i + 1,:] = numpy.matmul(numpy.transpose(self.j_svd_vt), temp_work3)
    
#  Estimate the expected changes in g^2. Equation 22 in Hanson et. al.
#  doi: 10.1088/0029-5515/49/7/075031
        exp_dg2 = numpy.empty(len(self.j_svd_w) + 1, dtype=float)
        self.delta_a_len = numpy.empty(len(self.j_svd_w) + 1, dtype=float)
        exp_eff = numpy.empty(len(self.j_svd_w) + 1, dtype=float)
        for i in range(0, len(exp_dg2)):
            exp_dg2[i] = self.get_exp_dg2(self.delta_a[i,:])
            self.delta_a_len[i] = math.sqrt(numpy.dot(self.delta_a[i,:], self.delta_a[i,:]))
            
#  Equation 23 in Hanson et. al. doi: 10.1088/0029-5515/49/7/075031
            exp_eff[i] = abs(exp_dg2[i])/self.delta_a_len[i]
        
        
#  Although marginal efficiencies are not strictly defined for the Steepest
#  Descent (index 0) and just one singular value cases, for convenience define
#  them here.
        marg_exp_eff = numpy.empty(len(self.j_svd_w) + 1, dtype=float)
        marg_exp_eff[0:1] = exp_eff[0:1]
        for i in range(2, len(exp_dg2)):
            d_len = max(self.delta_a_len[i] - self.delta_a_len[i - 1], 1.0e-10)
            marg_exp_eff[i] = abs(exp_dg2[i]) - abs(exp_dg2[i - 1])/d_len

#  Check for cutoffs.
        largest_w = 1.0
        if self.j_svd_w[0] > 0:
            largest_w = self.j_svd_w[0]

#  Find the largest number of singular values to use.
        for i in range(len(self.j_svd_w), 1, -1):
            meets_cut =               (self.j_svd_w[i - 1]/largest_w >= self.cut_svd)
            meets_cut = meets_cut and (exp_eff[i]                    >= self.cut_eff)
            meets_cut = meets_cut and (marg_exp_eff[i]               >= self.cut_marg_eff)
            meets_cut = meets_cut and (self.delta_a_len[i]           >= self.cut_delta_a)
            meets_cut = meets_cut and (numpy.abs(exp_dg2[i])         >= self.cut_dg2)
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
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_driver: get_exp_dg2')
        
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
    
#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver lm_step method. Determines the normalized change in
#  parameters for a Levenberg-Marquardt step.
#
#-------------------------------------------------------------------------------
    def lm_step(self, step_size):
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_driver: lm_step')

        if step_size > 0.0 and self.delta_a_len[self.k_use] > step_size:
            ut_dot_e = numpy.matmul(self.e, self.j_svd_u)
        
#  Find the L-M parameter lambda that corresponds to a step length of step_size.
            _lambda = self.lm_get_lambda(step_size, ut_dot_e)

#  Find the step.
            return numpy.matmul(ut_dot_e[0:self.k_use]*self.j_svd_w[0:self.k_use]/(self.j_svd_w[0:self.k_use]**2.0 + _lambda), self.j_svd_vt[0:self.k_use])
        else:
            return self.delta_a[self.k_use]

#-------------------------------------------------------------------------------
#
#  QUASI-NEWTON Driver lm_get_lambda method. Finds a Levenberg-Marquardt
#  parameter lambda coresponding to the length of step_size.
#
#-------------------------------------------------------------------------------
    def lm_get_lambda(self, step_size, ut_dot_e):
        ScreenWriter.screen_output(self, 'verbose', 'quasi_newton_driver: lm_get_lambda')

#  Define a default value incase the root find fails. Note lambda is a python
#  builtin add an underscore to avoid this.
        _lambda = (self.j_svd_w[0]*self.delta_a_len[self.k_use]/step_size)**2.0
        
        f_sqrd = ut_dot_e[0:self.k_use]**2.0
        step_size_sqrd = step_size**2.0

#  Define a lambda function for the root finder.
        f = lambda x: numpy.sum(f_sqrd*(self.j_svd_w[0:self.k_use]/(self.j_svd_w[0:self.k_use]**2.0 + x))**2.0) - step_size_sqrd

#  Find the bracketing values of lambda. Look for a small value of lambda, to
#  give a positive function value. f(0) should be greater than zero since the
#  step size at k_use is larger than step_size.
        f_a = f(0)
        if f_a < 0.0:
            return _lambda

        lambda_b = self.j_svd_w[0]**2.0
        f_b = f(lambda_b)
        for i in range(0, 20):
            if f_b <= 0.0:
                break
            lambda_b = 4.0*lambda_b
            f_b = f(lambda_b)

        if f_b > 0.0:
            return _lambda

        if f_a*f_b > 0.0:
            return _lambda

#  Found to intervals that bracket the roots. Now find the root.
        return optimize.brentq(f, 0.0, lambda_b)
