#! /usr/bin/env python

#  Adaption of init_rf_ic_aorsa.py (BH070415)

import sys
import os
import subprocess
import getopt

IPS_ROOT=''
SIM_ROOT=''
COMPONENT_CLASS = 'fp/'
COMPONENT_SUBCLASS = ''
COMPONENT_NAME = 'cql3d/'
INPUT_FILES = ['cqlinput', 'grfont.dat', 'eqdskin']

def printUsageMessage():
    print 'Usage: %s --ipsroot=FULL_IPS_ROOT_PATH --simroot=FULL_PATH_TO_CURRENT_SIMULATION' % (sys.argv[0])


def main(argv=None):
    global IPS_ROOT, SIM_ROOT, COMPONENT_CLASS, COMPONENT_SUBCLASS, COMPONENT_NAME, INPUT_FILES
# Parse command line arguments
    if argv is None:
        argv = sys.argv
        try:
            opts, args = getopt.gnu_getopt(argv[1:],'', ["ipsroot=", "simroot="])
        except getopt.error, msg:
            printUsageMessage()
            return 1
    for arg,value in opts:
        print arg, value
        if (arg == '--ipsroot'):
            IPS_ROOT = value
        elif (arg == '--simroot'):
            SIM_ROOT = value
    if (IPS_ROOT == '' or SIM_ROOT == ''):
        printUsageMessage()
        return 1
#
#Assumptions:
#  1- Initial input files (files that do not change with each time step) are located
#     along with the component sources (e.g. in IPS_ROOT/components/fp/cql3d)
#  2- Input files are copied (or soft linked) to the work directory of
#     the current simulation run
#  3- The work directory for the current component in the current simulation run is
#      SIM_ROOT/work/COMPONENT_CLASS/COMPONENT_SUBCLASS/COMPONENT_NAME
#      (e.g. IPS_RUN_XYZ/work/fp/cql3d)  [BH: using COMPONENT_SUBCLASS='']
#
#BH   inputFiles_src = IPS_ROOT + '/components/' + COMPONENT_CLASS +'/' + COMPONENT_NAME
    inputFiles_src = IPS_ROOT + '/components/' + COMPONENT_CLASS + COMPONENT_NAME + '/src/'
    print 'inputFiles_src = ' + inputFiles_src
#
# Note: cqlinput has a static part for given problem AND a  dynamic part that is overwritten
#       as part of the call to prepare_cql3d_input

#
# Check existence and/or create working directory for the current run
    workdir = SIM_ROOT + '/work/'+ COMPONENT_CLASS + COMPONENT_SUBCLASS + COMPONENT_NAME
    try:
        os.chdir(workdir)
    except OSError, (errno, strerror):
        print 'Directory %s does not exist - will attempt creation' % (workdir)
        try:
            os.makedirs(workdir)
        except OSError, (errno, strerror):
            print 'Error creating directory %s : %s' % (workdir, strerror)
            return 1
    print 'INPUT_FILES = '
    INPUT_FILES
#    print os.getcwd()
# Copy static files into working directory
#    print  [inputFiles_src + file for file in INPUT_FILES]
    for f in [inputFiles_src + file for file in INPUT_FILES] :
        print f
        cmd = ['/bin/cp', f, workdir]
        retcode = subprocess.call(cmd)
        if (retcode != 0):
            print 'Error copying file %s to work directory %s' %(f, workdir)
            return 1

if __name__ == "__main__":
    sys.exit(main())
