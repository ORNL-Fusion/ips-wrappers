#! /usr/bin/env python

"""
monitor_comp_stand_alone.py 3-19-2012 DBB 

A script to drive the IPS monitor component.  Invoke this script with the path to a
plasma state file and it will generate a monitor file.  Of course there will only be
one time point in the monitor file.

This script needs to be in the directory with the monitor component script and also needs
the associated dummy IPS component script (named component.py) in the same directory.

There are several versions of the IPS monitor component out there.  I used to produce a
new stand alone component whenever the monitor component changed.  This version just
wraps the actual IPS component using the dummy version of the IPS component class mentioned
above.  So this script need not change.  EXCEPT that you have to be sure that the import
line below refers to the monitor_comp module you actually are using.  The IPS component
also generates a python pickle file "monitor_restart" which has no relevance to the stand
alone operation.

We get too soon old and too late schmart.
"""


import sys
import os
import shutil

from monitor_comp import *   # Make sure to import the right version of monitor_comp


#----------------------------------------------------------------------------------------------



if __name__ == '__main__':


    if debug: print(' sys.argv = ', sys.argv)

    if len(sys.argv) != 2:
        print(' sys.argv = ', sys.argv)
        print('Usage: this script takes one command line argument -> a Plasma State file name')
        sys.exit(1)

    stateFile_arg = sys.argv[1]

    # Generate the monitor file name <stateFile>.cdf -> <stateFile>_monitor.nc
    # Chop off extension of stateFile
    state_part = stateFile_arg[0: stateFile_arg.rfind('.')]
    monitor_fileName = state_part + '_monitor.nc'
    services = None
    config = None
    mf = monitor(services, config)
    mf.INCLUDE_ZONE_BASED_PROFILES = True

    print('monitor_file init')
    print(mf.init_monitor_file(stateFile_arg,  0.0))

    print('monitor_file step')
    print(mf.update_monitor_file(stateFile_arg,  1.0))

    shutil.move('monitor_file.nc', monitor_fileName)