#! /usr/bin/python

#  Simple-minded script running cql3d for 20 step in restart mode
#  IPS_PHYS_ROOT is set by 
import sys   #things like sys.path, sys.argv[]
import os    #getcwd(),listdir(),chown,rename,mkdir,rmdir,system(command),etc.
import subprocess
import shutil


#os.system('cp cqlinput.2r.1 cqlinput')
#os.system('/home/bobh/cql3d/cql3d_cvs/code/cql3d_bipointered/cqlp>log.2r.1')

#os.system('cp cqlinput.2r.2 cqlinput')
#os.system('/home/bobh/cql3d/cql3d_cvs/code/cql3d_bipointered/cqlp>log.2r.2')


os.system('cp cqlinput.1r.1 cqlinput')
os.system('$HOME/cql3d_cswim_svn/trunk/xcql3d>log.2r.1')
os.system('cp ./iter_runaway.1r.1.nc distrfunc.nc')

os.system('cp cqlinput.1r.2 cqlinput')
os.system('$HOME/cql3d_cswim_svn/trunk/xcql3d>log.1r.2')

#Results on franklin:
#Time for 1 step per run, two steps is 0m41.3s
#Time for 2 steps with cql3d [no file write/read] is 0m40.0s
#  Very good result:  File overhead not a problem.
