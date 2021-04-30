
# basic_init_.py

A simple initializer whose purpose is to collect a complete set of initial state files and
stage them to the /work/state/ directory. It is simplified and adapted from
generic_ps_init.py, but eliminating reference to many features specific to plasma physics.
 The immediate application is to simple, example simulations which do not make use of the
SWIM Plasma State system. The SWIM Plasma State need not be used at all, but can be.

This script first touches all the files listed as STATE\_FILES in the config file.  If
config variable INIT\_MODE = TOUCH\_ONLY that is all that is done, then returns. If
INPUT\_FILES are listed in the [init] section of the config file these are then staged to
the working directory thereby overwriting the dummy files generated before.  It may be
convenient to initialize a state file from an input file of a different name, therefore a
service is provided to copy a file to another name before updating state.  This is
controlled by config parameters COPY\_FILES and COPIED\_FILES\_NEW\_NAMES

COPY\_FILES = list of files to be copied 
COPIED\_FILES\_NEW\_NAMES = list of new names for files, must match COPY_FILES

The term CURRENT\_STATE refers to a SWIM PLASMA\_STATE file if it is being used.  If a
CURRENT\_STATE is specified in the simulation config file, this script adds the Plasma
State variables to it -> run_id and time loop variables tinit, and tfinal.

If more work needs to be done before the individual components do their own init, one can
specify an INIT\_HELPER\_CODE (full path) in the config file which if present will be
executed here. If input files are needed for the helper code they must also be specified
in the [init] section of the config file.
