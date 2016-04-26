# Run using 
# aprun -n 24 python solpsnode.py

import os, sys, shutil
from Namelist import Namelist
from mpi4py import MPI
import subprocess

baserun = 'baserun'
templaterun = 'templaterun'
runbasename = 'thread'
nml_boundary = 'b2.boundary.parameters'

comm = MPI.COMM_WORLD
myRank = comm.Get_rank()
nthreads = comm.Get_size()

if myRank == 0:
    print 'nThreads: ' + str(nthreads)

density_min =  1e18
density_step = 0.1e18
density_nsteps = 12 

for thread in range(0,density_nsteps):
    if thread%nthreads==myRank:
        #print "myRank: " + str(myRank) + " thread: " + str(thread)
        thisRun = runbasename + "_" + str(thread).zfill(2)            
        print "Staging the template run in " + thisRun        
        shutil.copytree(templaterun,thisRun,symlinks=True)
        os.chdir(thisRun)
        shutil.copy(nml_boundary,nml_boundary + "_original")
        obj = Namelist(nml_boundary)
        obj["BOUNDARY"]["CONPAR(0,1,1)"][1] = density_min + density_step * thread
        obj.write(nml_boundary)
        os.chdir('../')

for thread in range(0,density_nsteps):
    if thread%nthreads==myRank:
        thisRun = runbasename + "_" + str(thread).zfill(2)            
        launchSOLPScmd = "./launchSOLPS5.sh $(readlink -f baserun) $(readlink -f "+thisRun+")" 
        status = subprocess.call(launchSOLPScmd, shell=True)                    
 
