#! /usr/bin/env python

from  component import Component
import os
import shutil
import re
import sys
import time
import generate_ftridyn_input
import analyze_ftridyn_simulations
import numpy as np
import netCDF4
import pickle
import math

print 'The generate_ftridyn_input path in ftridyn_worker_taskpool is: '
print os.path.abspath(generate_ftridyn_input.__file__)

class ftridynWorker(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)
        
	self.ftridyn_exe = self.FTRIDYN_EXE

    def init(self, timeStamp=0.0,**keywords):
        return
 
    def step(self, timeStamp=0.0,**keywords):
        
	print('ftridyn_worker: step (mpi version)')
        self.services.stage_plasma_state()
	#shutil.copyfile(self.FTMPI_IN,self.services.get_working_dir()+'/ftMPI.in')
        d = {}
        d['beam']=str(self.BEAM)
        d['target']=str(self.TARGET)
        d['Escale']=str(self.ESCALE)
        d['exe']=str(self.FTRIDYN_EXE)
        
        d['nE'] = int(self.NE)
	d['nA'] = int(self.NA)
	d['nR'] = int(self.NR)
	d['nEdist'] = int(self.NEDIST)
	d['nAdist'] = int(self.NADIST)
	d['nH']= int(self.NH)
        
	d['energy_start']=float(self.E_START)
	d['energy_end']=float(self.E_END)
	d['angle_start']=float(self.A_START)
	d['angle_end']=float(self.A_END)
	d['roughness_start']=float(self.R_START)
	d['roughness_end']=float(self.R_END)
	d['maxEdist']=float(self.MAXEDIST)
	f = open('ipsFTmpi.pkl', 'w')
        pickle.dump(d,f)
	f.close()

        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            'python',self.FTMPI,'--dictionary=ipsFTmpi.pkl',
                                            logfile='ftmpi.log') #,ppn=1)
        
	#monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('ftmpi_comp: step failed.')
        
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
