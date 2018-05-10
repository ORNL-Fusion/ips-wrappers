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
        #get name of FTridyn executable from config file
        #This is actually a shell script that calls FTridyn and pipes
        #the input to the executable
        self.ftridyn_exe = self.FTRIDYN_EXE

    def init(self, timeStamp=0.0,**keywords):
        #stage plasma state files for use on execution of FTridyn
        self.services.stage_plasma_state()
        #unpack arguments from init call in driver
        ffilename = keywords["ffilename"]
        beam = keywords["beam"]
        target = keywords["target"]
        nH = keywords["nH"]
        energy = keywords["eArg"]
        angle = keywords["aArg"]
        roughness = keywords["dArg"]

        for i in range(len(energy)):
            for j in range(len(angle)):
                for k in range(len(roughness)):
                    generate_ftridyn_input.beam_and_target(ffilename,beam,target,number_histories=nH,incident_energy=energy[i],incident_angle=angle[j],fractal_dimension=roughness[k])
                    pathString = "FTRIDYN_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
                    if not os.path.exists(pathString):
                        os.makedirs(pathString)
                    shutil.copyfile(ffilename+"0001.IN", pathString+"/"+ffilename+"0001.IN")
                    
                    surfname = str(1)+'p'+str(int((roughness[k]-0.99999)*1000.0))+'.surf'
                    shutil.copyfile(surfname, pathString+"/"+surfname)
                    if energy[i] < 0:
                        #copy energy distribution file
                        shutil.copyfile(self.INPUT_DIR+"/dist"+str(int(j))+".dat", pathString+"/"+ffilename+"0001.ED1") 
 
    def step(self, timeStamp=0.0,**keywords):
        start_time = time.time()
        print('ftridyn_worker: step (task pool version)')
        
        ffilename = keywords["ffilename"]
        beam = keywords["beam"]
        target = keywords["target"]
        nH = keywords["nH"] #number_histories=nH
        energy = keywords["eArg"]
        angle = keywords["aArg"]
        roughness = keywords["dArg"]
        file_namelist = []

        nEgrid = 100
        maxE = 100.0
        nAgrid = 50
        sputt = np.zeros(shape=(len(energy),len(angle)))
        refl = np.zeros(shape=(len(energy),len(angle)))
        eDistEgrid = np.linspace(0.0,maxE-maxE/nEgrid,nEgrid) 
        cosDistAgrid = np.linspace(0.0,90.0-90.0/nAgrid,nAgrid) 
        cosXDistribution = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosYDistribution = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosZDistribution = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosXDistributionRef = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosYDistributionRef = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        cosZDistributionRef = np.zeros(shape=(len(energy),len(angle),nAgrid)) 
        eDistribution = np.zeros(shape=(len(energy),len(angle),nEgrid)) 
        eDistributionRef = np.zeros(shape=(len(energy),len(angle),nEgrid)) 
        
        nFTruns = len(file_namelist)
        nFTrunsPerNode = int(self.FTMPI_PPN)
        nFTrunNodes = int(self.FTMPI_NODES)
        if(nFTruns <= nFTrunsPerNode*nFTrunNodes): 
            nFTrunsPerNode = int(math.ceil(1.0*nFTruns/nFTrunNodes))
            nFTpoolTasks = nFTrunNodes
        else:
            nFTpoolTasks = int(math.ceil(1.0*nFTruns/nFTrunsPerNode))
        
        cwd = self.services.get_working_dir()
        #pool = self.services.create_task_pool('pool')
        for i in range(len(energy)):
            for j in range(len(angle)):
                for k in range(len(roughness)):
                    pathString = "FTRIDYN_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
                    file_namelist.append(pathString)
                    #self.services.add_task('pool', 'task'+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k]), 1, cwd+"/"+pathString, self.FTRIDYN_EXE, ffilename+"0001.IN",logfile='task_'+pathString+'.log' )
        #ret_val = self.services.submit_tasks('pool')
        #print 'ret_val = ', ret_val
        #exit_status = self.services.get_finished_tasks('pool')
        #print exit_status
        #self.services.remove_task_pool('pool')
        nFTruns = len(file_namelist)
        nFTrunsPerNode = int(self.FTMPI_PPN)
        nFTrunNodes = int(self.FTMPI_NODES)
        if(nFTruns <= nFTrunsPerNode*nFTrunNodes): 
            nFTrunsPerNode = int(math.ceil(1.0*nFTruns/nFTrunNodes))
            nFTpoolTasks = nFTrunNodes
        else:
            nFTpoolTasks = int(math.ceil(1.0*nFTruns/nFTrunsPerNode))

        with open('ftridyn_file_namelist.pkl', 'wb') as f:
            pickle.dump(list(reversed(file_namelist)), f)

        pool = self.services.create_task_pool('pool')
        for i in range(nFTpoolTasks):
            self.services.add_task('pool', 'task'+str(i), nFTrunsPerNode, cwd, 'python',self.FTMPI_EXEC,str(i), str(nFTrunsPerNode),task_ppn= nFTrunsPerNode,logfile='task_pool'+str(i)+'.log' )
        ret_val = self.services.submit_tasks('pool')
        print 'ret_val = ', ret_val
        exit_status = self.services.get_finished_tasks('pool')
        print exit_status
        self.services.remove_task_pool('pool')

        exec_time = time.time()
        print("Execution of FTRIDYN Cases took --- %s seconds ---" % (exec_time - start_time))
        pool_py = self.services.create_task_pool('pool_py')
        for i in range(nFTpoolTasks):
            self.services.add_task('pool_py', 'task'+str(i), nFTrunsPerNode, cwd, 'python',self.FTMPI,str(i), str(nFTrunsPerNode),task_ppn= nFTrunsPerNode,logfile='task_'+str(i)+'.log' )
        ret_val = self.services.submit_tasks('pool_py')
        print 'ret_val = ', ret_val
        exit_status = self.services.get_finished_tasks('pool_py')
        print exit_status
        self.services.remove_task_pool('pool_py')

        #task_id = self.services.launch_task(self.FTMPI_TASKS,self.services.get_working_dir(),'python',self.FTMPI,task_ppn=self.FTMPI_PPN)
        ##monitor task until complete
        #if (self.services.wait_task(task_id)):
        #    self.services.error('ftmpicomp: step failed.')

        spyl_file = ffilename+'SPYL.DAT'
        #driver_out = self.services.get_config_param('EA_OUTPUT')
        #fid = open(driver_out,'a')
        #for i in range(len(file_namelist)):
        totalIndex=0
        for i in range(len(energy)):
            for j in range(len(angle)):
                pathString = file_namelist[totalIndex]
                yr = np.loadtxt(pathString+"/"+"YR.out", dtype='float')
                sputt[i,j] = yr[0]
                refl[i,j] = yr[1]
                nX = np.loadtxt(pathString+"/"+"nX.out", dtype='float')
                cosXDistribution[i,j,:] = nX
                nXref = np.loadtxt(pathString+"/"+"nXref.out", dtype='float')
                cosXDistributionRef[i,j,:] = nXref
                nenergy = np.loadtxt(pathString+"/"+"energy.out", dtype='float')
                eDistribution[i,j,:] = nenergy
                nenergyRef = np.loadtxt(pathString+"/"+"energyRef.out", dtype='float')
                eDistributionRef[i,j,:] = nenergyRef
                totalIndex = totalIndex+1
        #for i in range(len(energy)):
        #    for j in range(len(angle)):
        #        for k in range(len(roughness)):
        #            pathString = "FTRIDYN_"+str(energy[i]) + "_"+str(angle[j])+"_" + str(roughness[k])
        #            print('doing analysis on ' , pathString+"/"+spyl_file)
        #            WW=analyze_ftridyn_simulations.ftridyn_output_data(pathString+"/"+ffilename,'0001')
        #            thisSpyl = WW.calculate_total_sputtering_yield()
        #            print('sputtering yield', thisSpyl)
        #            sputt[i,j] = thisSpyl
        #            thisRefyl = WW.calculate_total_reflection_yield()
        #            print('reflection yield', thisRefyl)
        #            refl[i,j] = thisRefyl
        #            nP, nX, binsX,nY,binsY,nZ, binsZ = analyze_ftridyn_simulations.plot_sputtered_angular_distributions(pathString+"/"+ffilename,nAgrid)
        #            cosXDistribution[i,j,:] = nX
        #            cosYDistribution[i,j,:] = nY
        #            cosZDistribution[i,j,:] = nZ
        #            nP, nX, binsX,nY,binsY,nZ, binsZ = analyze_ftridyn_simulations.plot_reflected_angular_distributions(pathString+"/"+ffilename,nAgrid)
        #            cosXDistributionRef[i,j,:] = nX
        #            cosYDistributionRef[i,j,:] = nY
        #            cosZDistributionRef[i,j,:] = nZ
        #            nPenergy, nenergy, binsenergy = analyze_ftridyn_simulations.plot_sputtered_energy_distributions(pathString+"/"+ffilename,nEgrid)
        #            eDistribution[i,j,:] = nenergy
        #            analyze_ftridyn_simulations.plot_implantation_profile(pathString+"/"+ffilename)
        #            nPenergyRef, nenergyRef, binsenergyRef = analyze_ftridyn_simulations.plot_reflected_energy_distributions(pathString+"/"+ffilename,nEgrid)
        #            eDistributionRef[i,j,:] = nenergyRef

                    #fid.write(" ".join([str(energy[i]),str(angle[j]),str(roughness[k]),'  ',str(thisSpyl),'\n']))
        #fid.close()
        rootgrp = netCDF4.Dataset("ftridyn.nc", "w", format="NETCDF4")
        ne = rootgrp.createDimension("nE", len(energy))
        na = rootgrp.createDimension("nA", len(angle))
        nedistgrid = rootgrp.createDimension("nEdistBins", nEgrid)
        nadistgrid = rootgrp.createDimension("nAdistBins", nAgrid)
        spyld = rootgrp.createVariable("spyld","f8",("nE","nA"))
        rfyld = rootgrp.createVariable("rfyld","f8",("nE","nA"))
        ee = rootgrp.createVariable("E","f8",("nE"))
        aa = rootgrp.createVariable("A","f8",("nA"))
        cosxdist = rootgrp.createVariable("cosXDist","f8",("nE","nA","nAdistBins"))
        cosydist = rootgrp.createVariable("cosYDist","f8",("nE","nA","nAdistBins"))
        coszdist = rootgrp.createVariable("cosZDist","f8",("nE","nA","nAdistBins"))
        cosxdistref = rootgrp.createVariable("cosXDistRef","f8",("nE","nA","nAdistBins"))
        cosydistref = rootgrp.createVariable("cosYDistRef","f8",("nE","nA","nAdistBins"))
        coszdistref = rootgrp.createVariable("cosZDistRef","f8",("nE","nA","nAdistBins"))
        edist = rootgrp.createVariable("energyDist","f8",("nE","nA","nEdistBins"))
        edistref = rootgrp.createVariable("energyDistRef","f8",("nE","nA","nEdistBins"))
        edistegrid = rootgrp.createVariable("eDistEgrid","f8",("nEdistBins")) 
        cosdistagrid = rootgrp.createVariable("cosDistAgrid","f8",("nAdistBins")) 
        ee[:] = energy
        aa[:] = angle
        edistegrid[:] = eDistEgrid
        cosdistagrid[:] = cosDistAgrid
        spyld[:] = sputt
        rfyld[:] = refl
        cosxdist[:] = cosXDistribution
        cosydist[:] = cosYDistribution
        coszdist[:] = cosZDistribution
        cosxdistref[:] = cosXDistributionRef
        cosydistref[:] = cosYDistributionRef
        coszdistref[:] = cosZDistributionRef
        edist[:] = eDistribution
        edistref[:] = eDistributionRef
        rootgrp.close()
        post_time = time.time()
        print("Processing of FTRIDYN Cases took --- %s seconds ---" % (post_time - exec_time))
        #updates plasma state FTridyn output files
        self.services.update_plasma_state()
  
    def finalize(self, timeStamp=0.0):
        return
    
