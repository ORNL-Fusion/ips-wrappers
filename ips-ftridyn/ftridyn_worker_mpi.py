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
import periodic
from collections import defaultdict

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
        d = defaultdict(list)
	if(int(self.GET_SPECIES)):
	    specZs=[];
	    specList = np.loadtxt('speciesList.txt', dtype='float', skiprows=1)
	    print specList, len(specList)
	    for i in range(len(specList)):
	        print "0",specZs
	        if ((specList[i,1] ==1.0) and (specList[i,3]>0.0)):
		    if(specList[i,1] == 1.0 and specList[i,2]==1.0):
		        d['beam'].append('H')
		        specZs.append(specList[i,1])
	                print "2",periodic.element(specList[i,1]).symbol
		    if(specList[i,1] == 1.0 and specList[i,2]==2.0):
		        d['beam'].append('D')
		        specZs.append(specList[i,1])
	                print "2",periodic.element(specList[i,1]).symbol
		    if(specList[i,1] == 1.0 and specList[i,2]==3.0):
		        d['beam'].append('T')
		        specZs.append(specList[i,1])
	                print "2",periodic.element(specList[i,1]).symbol
	        elif((specList[i,1] not in specZs) and (specList[i,3]>0.0)):
	                print periodic.element(specList[i,1]).symbol
		        specZs.append(specList[i,1])
		        d['beam'].append(periodic.element(specList[i,1]).symbol)
	                print "3",periodic.element(specList[i,1]).symbol
	    #shutil.copyfile(self.FTMPI_IN,self.services.get_working_dir()+'/ftMPI.in')
        d['beam'].append(str(self.BEAM))
        d['target']=str(self.TARGET)
        d['Escale']=str(self.ESCALE)
        d['exe']=str(self.FTRIDYN_EXE)
        
        d['nE'] = int(self.NE)
	d['nA'] = int(self.NA)
	d['nR'] = int(self.NR)
	d['nEdist'] = int(self.NEDIST)
	d['nEdist_ref'] = int(self.NEDIST_REF)
	d['nAdist'] = int(self.NADIST)
	d['nH']= int(self.NH)
        
	d['energy_start']=float(self.E_START)
	d['energy_end']=float(self.E_END)
	d['angle_start']=float(self.A_START)
	d['angle_end']=float(self.A_END)
	d['roughness_start']=float(self.R_START)
	d['roughness_end']=float(self.R_END)
	d['maxEdist']=float(self.MAXEDIST)
	d['maxEdist_ref']=float(self.MAXEDIST_REF)
	#f = open('ipsFTmpi.pkl', 'w')
        #pickle.dump(d,f)
	#f.close()
        print d
	ftMpiFile = open('ftMPI.in', 'w')
	for key in d:
	    if (key == 'beam'):
	        ftMpiFile.write("%s %s\n" %(key,' '.join(d[key])))
            else:		
	        ftMpiFile.write("%s %s\n" %(key,str(d[key])))
	ftMpiFile.close()	
        task_id = self.services.launch_task(self.NPROC,
                                            self.services.get_working_dir(),
                                            'python',self.FTMPI,task_ppn=self.TASK_PPN,#'--dictionary=ipsFTmpi.pkl',
                                            logfile='ftmpi.log') #,ppn=1)
        
	#monitor task until complete
        if (self.services.wait_task(task_id)):
            self.services.error('ftmpi_comp: step failed.')
        
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
