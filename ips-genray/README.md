
# New: rf_genray.py

Multi-frequency, config file programmable version: Batchelor (9/3/2018)
Supports both RFMODE = EC and LH.  Launcher power and launcher aiming programming is only 
implemented for ECH, so far.


Working notes:

9/3/2018 DBB
This version is adapted from the previous version genray_EC_p.py.  It retains the
programmability from the config file that was in genray_EC_p.py and the detection of zero
RF power that was in the 6/16/2014 version but that seemed to get dropped from the
genray_EC_p.py version.  The programming was implemented for ECH aiming angles, but no
such k-spectrum programming is implemented for Lower Hybrid, yet.
Switched NetCDF4.py from obsolete Scientific.IO.NetCDF.  Removed namelist editing routines
to /wrappers/utilities/simple_file_editing_functions.py.

11/25/2014 DBB
In many of the IPS components we have coding to allow the user to add a suffix to an input 
file name and then copy it to a generic file name, for example genray.in_myFile ->
genray.in.  However Bob Harvey has this functionality in the prepare_genray_input.f90
code.  He has "genraynml" as command line arg to prepare_genray_input, which then
writes to the generic name genray.in.  So to avoid confusion I have removed the coding
providing the suffix processing as used in other components.

Also found a bug in the exception handling coding for file copies, which I fixed 
(I think).

6/16/2014 DBB
About zero RF power:  We don't want to run an RF code when the RF power is zero.
The STEP function needs to produce a partial plasma state only containing the RF data
(i.e. ps_write_update_file not ps_store_plasma_state).  This really is done with
plasma state fortran code not with the python netcdf interface.
So I wrote a simple code called zero_RF_power to set
all RF source profiles in plasma state to zero and then write a partial plasma state.
For now this code lives in the GENRAY component directory and it gets built and
installed by the Makefile there.  The zero_RF_EC_power.f90 code is essentially the same
as the zero_RF_IC_power.f90, and the modifications to this python component are the
same as were needed for the rf_ic_toric_abr_mcmd.py component.


Programmable version: Batchelor (6/2/2013)
This version allows the aiming angles and power for multiple launchers to be set and
changed over time from within the simulation config file.  If this option is used, any 
programming of the launcher power from the EPA component is over-ridden.  This behavior 
is triggered by the presence a variable N_LAUNCHERS_PROGRAMMED in the [rf_genray] section 
of the simulation config file.  If N_LAUNCHERS_PROGRAMMED = 0, or is absent, the launcher
powers are taken from the current plasma state file.  
If > 0 the component looks for space delimited lists of time points for parameter changes 
and the associated parameters.  The config lists should be named: 
LAUNCHER1_TIMES,LAUNCHER1_alfast, LAUNCHER1_betast, LAUNCHER1_powtot_MW, 
LAUNCHER2_TIMES,LAUNCHER2_alfast, LAUNCHER2_betast, LAUNCHER2_powtot_MW,  etc

The parameter changes take effect the first time the simulation time at
the beginning of a time step (ps%t0) is >= LAUNCHERx_TIME then LAUNCHERx parameters are
changed.  LAUNCHERx_TIMES don't have to match simulation time step beginnings, but the 
changes won't take place until the beginning of the next time step.

Note: For ease of typing, enter powers in config file in MW, they get converted to Watts
in this component.  

The mechanism for changing the power is to modify the genray.in just before launching
the genray_EC step.  N_LAUNCHERS_PROGRAMMED must match the number of launchers found in
the input genray.in file.

EC version 0.0 1/10/2011
This version distinguishes between the different RF components that GENRAY can 
implement.  In this case the EC component.  So far the only place this affects
is in the writing the partial plasma state that the framework needs to merge.
#
Note: To merge plasma states the IPS framework expects the component to
produce a partial state file with the component-specific name:
<component>_<curr_state_filename>.  GENRAY works for EC, EC, and IC, i.e.
it can implement several different components.  Therefore process_genray_output_mcmd.f90
writes a partial state file with generic name 'RF_GENRAY_PARTIAL_STATE' and delegetes
to this python component script the copying of this to the proper 
component-specific update file name.  In this case RF_EC_<curr_state_filename>.

version 1.0 9/29/2010 (Batchelor)
This version adds checkpoint and restart functions and makes the exception 
handling more uniform.

version 0.0 3/1/2010 (Batchelor)
