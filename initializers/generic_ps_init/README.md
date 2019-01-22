# generic_ps_init.py and generic_ps_init.f90

The Swiss army knife of Plasma State initializers.  It produces the intial CURRENT_STATE
and optionally the initial CURRENT_EQDSK.

This version combines several previous initializer routines and extends them.  There are
4 modes of initialization which must be specified by the config file variable INIT_MODE

INIT_MODE = touch_only
This mode only does a touch on all of the files listed as plasma state files so the
framework will have a complete set.  It does not actually put data in the plasma state file.

INIT_MODE = minimal
This is exactly the same as the previous minimal_state_init.py. It produces a CURRENT_STATE
that is empty except for some metadata:
time variables - ps%t0, ps%t1, ps%tinit, and ps%tfinal
simulation identifiers - ps%tokamak_id, ps%shot_number, ps%run_id.
ps%Global_label is set to run_id_tokamak_id_shot_number.
This data is set for all initialization modes, but for 'minimal' this is all the data
included.

INIT_MODE = existing_ps_file
This copies an existing input plasma state file and optionally an existing eqdsk file to
CURRENT_STATE and CURRENT_EQDSK.  If the config parameter GENERATE_EQDSK is set to 'True'
the CURRENT_EQDSK file is generated from equilibrium data in the INPUT_STATE_FILE.
The INPUT_STATE_FILE and INPUT_EQDSK_FILE must be specified in the config file.

INIT_MODE = mdescr
This initializes all machine description data from a plasma state machine description
file, e.g. <tokamak>.mdescr, as specified by config parameter MDESCR_FILE. In addition
if a shot configuration file config parameter, SCONFIG_FILE, is specified, the shot config
data is also loaded into CURRENT_STATE.  Machine description and shot configuration files
are namelist files that can be read and loaded using Plasma State subroutines ps_mdescr_read()
and ps_scongif_read().  Note:  machine description and shot configuration do not define
the MHD equilibrium or plasma profiles, so the equilibrium and profiles must be specified 
during further component initializations.

INIT_MODE = mixed
This combines existing_ps_file and mdescr modes.  This copies an existing input plasma state 
file and optionally an existing eqdsk file to CURRENT_STATE and CURRENT_EQDSK as in 
existing_ps_file mode.  But initializations from MDESCR_FILE and SCONFIG_FILE are also added.
This is a convenient way to get profiles and equilibrium from a plasma state but have
different initializations for source component from the MDESCR and SCONFIG files.

Caution is advised.  If the MDESCR_FILE or SCONFIG_FILE attempts to reallocate any of the
arrays already allocated in the CURRENT_STATE file a Plasma State error will occur.
Therefore one must provide in the simulation config file a list of which components are
to be initialized from mdescr or sconfig.  The config parameter is MDESCR_COMPONENTS.  For
example  "MDESCR_COMPONENTS = LH EC NBI"

IPS Components are

1  PLASMA    Thermal plasma parameters; fluid profile advance
2  EQ        MHD equilibrium
3  NBI       Neutral beam
4  IC        Ion cyclotron heating
5  LH        Lower Hybrid heating and current drive
6  EC        Electron cyclotron heating and current drive
7  RUNAWAY   Runaway electrons
8  FUS       Fusion product fast ions
9  RAD       Radiated Power & impurity transport
10  GAS      Neutral Gas sources & transport
11  LMHD     Linear MHD stability
12  RIPPLE   TF field ripple
13  ANOM     Anomalous transport

Except for possibly mode = existing_ps_file, all modes call on the fortran helper code
generic_ps_file_init.f90 to interact with the Plasma State. The fortran code is also used
in existing_ps_file mode to extract the CURRENT_EQDSK when GENERATE_EQDSK = true.