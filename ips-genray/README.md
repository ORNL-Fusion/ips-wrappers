May 8, 2008  (BH)
===========
Compiled latest genray at
viz:/u/bharvey/genray_cvs/VIZ_genray_v7-2_080423/xgenray

See viz:/u/bharvey/genray_cvs/Coupling_to_PS
for README_onetwo-genray and some coding concerning
previous coupling of genray to ONETWO transport code.

GENRAY side of coupling is mainly in 
viz:/u/bharvey/genray_cvs/VIZ_genray_v7-2_080423/partner.f
We have, from genray.in namelist file:
!--------------------------------------------------------------------------
! partner  = 'disabled' to use input profiles from genray.dat or genray.in
!          = 'genray_profs_in.nc'  to use plasma profile data from the netCDF
!                       file: genray_profs_in.nc  created by transport code
!                       or otherwise;
!                  and to output power and current profiles on transport
!                      code radial grid to file: genray_profs_out.nc
!                      and genray_profs_out.txt (text file).
!          = 'genray_profs_in.txt'  to use plasma profile data from the text
!                      file: genray_profs_in.txt  created by transport code
!                      or otherwise;
!                  and to output power and current profiles on transport
!                      code radial grid to file: genray_profs_out.nc
!                      and genray_profs_out.txt (text file).
!------------------------------------------------------------------------


Start with partner='genray_profs_in.txt' for initial coupling
to plasma state.

prepare_genray_input.f90 will write genray_profs_in.txt/.nc using
PS data.  
process_genray_output.f90 will read genray_profs_out.txt/.nc [from
genray], and put the power deposition and current drive in the PS.

May 9, 2008
===========

For prepare_genray_input.f/process_genray_output.f start from
/p/swim1/berryla/ips/trunk/components/fp/rfmin/cql3d/prepare_fp_rfmin_cql3d_input.f and  process_fp_rfmin_cql3d_output.f
==>  prepare_genray_input.f/process_genray_output.f


viz:/u/bharvey/genray_cvs/Coupling_to_PS/onetwo_genray_related.txt
contains coding used  for coupling genray to onetwo [write of
           genray_profs_in.txt/genray_profs_in.nc
	   and reaad of genray_profs_out.txt/genray_profs_out.nc]

read_write_genray_input.f: From genray
partner.f: From genray, contains coding related to read
           of genray_profs_in.txt/genray_profs_in.nc
	   and write of genray_profs_out.txt/genray_profs_out.nc

See Smirnov email, 03/10/07 07:17 am, re
Latest genray with declaration of namelist variables


May 16, 2008
============
Have compiled version of prepare_genray_input.f90, on viz.
Accesses plasma state, and reads all genray namelist.
Working with Doug McCune to interpolate electrons and alla
ion species (abridged impurities, but complete) radial
density and temperature profiles to bin
boundaries, for passage to genray through genray_profs_in.txt and .nc
files.  For cases where electron and Zeff profiles are sufficient
(EC, EBW), don't need ion profiles.  For other cases (ICRF, HHFW, LH)
need ions.   Use abridged impurities (up to two impurities).
Will also have beam, RF minority, and fusion product nonthermal
distributions.   Allowing for extreme cases, we compile genray
for up to 11 ion species (12 species total, including electrons).


May 20, 2008
============

Adapting coding by McCune for gathering all ion species into a single
array of densities and temperatures, to transfer to genray_profs_in file.


May 23, 2008
============

The main complication in application of the Plasma State to genray
coupling is use of the abridged species set.  In general the plasma
state may provide data for a much larger list of ions--- perhaps
including multiple charge states of several impurities--- than is
desireable to be treated in genray.   The abridged set of ions reduces
the impurity species to .le.2 species, conserving Zeff and charge
neutrality.
See Abridged_Species_PS2.pdf, Abridged_PS2_ps_tha_summary.pdf,
McCune_frm_080516.eml, McCune_frm_080523.eml.

July 17, 2008
=============
Have working prepare_genray_input.f90.4 
and          process_genray_output.f90.3,
and
three sets of tests:
Cmod_LH/ ITER_EC_Budny/ NSTX_FW
[Have not worked yet on the D3D_ECCD case].

A new read_write_genray_input.f (old designated read_write_genray_input.f.0)
will be incorporated into prepare_genray_input.f90(.5).
The new read_write_genray_input.f simplifies (modularizes, i.e. separate
subroutines) the treatment of switching from mixed units (genray.dat 
conventions used in genray) to genray.in MKSA units.

The codes will also be updated to genray version="genray_v7-6_080714".
(The only change in the includes is the version parameter in param.i.)



Updating read_write_genray_input.f from genray distribution
===========================================================
-Move present read_write_genray_input.f to read_write_genray_input.f_xxxxxx
 where xxxxx is file date.
-Copy read_write_genray_input.f from distn to pwd.
-Compare new and old versions to check changes look OK. Don't forget to
 comment out indicated region.
-cd genray_includes, and diff with distn.   Update from distn.


August 31, 2008
===============
Updated prepare_genray_input.f90/process_genray_output.f90.
read_write_genray_input.f is from latest genray_v7-8_080801 (passes
regression tests, and including propagation outside LCFS).
[freqncy for EC case changed in PS from MHz to GHz.  Other small
changes].

Updated genray_includes/

Update tests/ITER_EC_Budny/test_1src

Other EC variables have been added to the PS Version_ID 2.011,
as indicated by runtime messages.   
Need an updated PS to apply these changes with, before modifying
prepare_genray_input.f90.

Other tests need re-checking.


Oct. 8, 2008
============

Bugs fixed in genray, per a_change.i in the distribution
in genray_v7-11_081007.tar.gz.

Putting this distribution in 
/p/swim1/bharvey/ips/trunk/components/rf/genray/genray
[i.e., ../genray], and adding it to SWIM svn repository.

Updating prepare_genray_input.f90 to frequencies in Hz,
for EC/LH/IC components for Plasma State PSv2.014,
per McCune Email, Oct. 9, 2008.


Feb. 11, 2009
=============
Update prepare_genray_input.f90, for recent changes in PS.
Compiled and loaded prepare_genray_input.f90/process_genray_output.f90
Added tests/ITER_EC_Budny/test_5src_NEW/ and 
tests/ITER_EC_Budny/test_5src_NEW_short/
Scripts in these directories run the
prepare_genray_input/xgenray_nowrite/process_genray_output sequences.
Full test in test_5src_NEW with 100 rays takes about 4.5 min on viz,
shorter test with 1 ray per source (5 sources) takes about 17 secs.

Feb. 13, 2009
=============
Added ipsmode='init' or 'step' option to the command line in
prepare_genray_input.   The 'init' option initializes some
rf related dimensions and variables in the plasma state.  There
is also the possibility for re-dimensioning there variables if
the radial dimension from the genray.in namelist file is different
from already set values in the namelist (as may occur if the 
PS file was obtained from another code).

April 9, 2008
=============
viz.pppl.gov
Updated prepare_genray_input.f90 and added genray_includes files
consistent with latest svn tagged genray version:
/p/swim1/bharvey/cswim/genray/tags/genray_v7-16_090406

Added environmental variable to facilitate use of tagged genray:
GENR_DISTN='/p/swim1/bharvey/cswim/genray/tags/genray_v7-16_090406'
Will also use this on franklin and jaguar.

Ran test of prepare/process genray given in
./tests/ITER_EC_Budny/test_5src_NEW_short/test_ITER_EC.sh_5src_NEW_short


Sept 27, 2009
=============
franklin.nersc.gov
Updated genray .f and genray_include files (added several also)
for genray version="genray_v8_1_090921", dynamic dimensioning version.






==========================================================================
==========================================================================

Seven-step program for updating my_ips_trunk/components/rf/genray/src
=====================================================================
after changes of genray affecting this area:
============================================
 
1. svn update genray svn distribution
2. If relevant changes, go to prep-proc genray (this area)
3. Execute update_prep-proc_genr_distn.sh, updating files in 
   this area and ./genray_includes
   (If have to update list of files in ./genray_includes, to the .sh file.)
4. make clean; make  gives update .o and prep-/proc-genray executables
5. Commit changes to svn
6. Go genray svn distrn;  update makefile_swim if necessary;
   make -f makefile_swim
7. commit any changes to makefile_swim to svn.




Sept. 29, 2009
==============
DBB may need additional command line inputs for 
prepare_genray_input.f90:  cur_state_file, cur_eqdsk.
process_genray_output.f90:  cur_state_file
These have been added (needs checking) as of ORNL
coding camping, in
prepare_genray_input.f90_new/process_genray_output.f90_new.
These changes are not yet committed.


Feb. 10, 2010
=============
Updated genray files, and remade mainline
prepare_genray_input.f90/process_genray_output.f90




Feb. 13, 2010
=============
Updated components/state directory, for latest PS.  Make looks OK.
Had to add ,ierr argument to ckerr subroutine.  Don't know how
it worked before?:  Maybe ierr was 0 by default in the subroutine?.
Reproduced genray.in and results thereof for
tests/ITER_EC_Budny/test_5src_NEW_short case.


Sept. 29-30, 2011
=================
Trying to update read_write_genray_input_prep.f,
but addition of rz-mesh to treat general toroidal density-temp
variation outside the LCFS has added several new subroutines.
Checking into reworking genray flow chart.


Oct. 31, 2011
=============
Smirnov removed calculations related to the new rz-mesh
from read_write_genray_input_prep.f, placing them
elsewhere in the code.  This enables 
read_write_genray_input_prep.f to contain only (almost)
data handling of the genray namelists, and cleans up
the storage to only (almost) namelist data only.

This is for 
version=genray_v10.0_111031


Nov. 11, 2011
=============
prepare_genray_input.f90 enhanced/debugged to handle
a TSC genrated PS file input, from Bonoli's C-MOD TSC/LSC case.

process_genray_output.f90 was modified according to
DBB process_genray_output_mcmd.f90 modifcation to output
partial PS file.
For time being, after make of process_genray_output,
copying it to process_genray_output_mcmd.
Put the executables in /global/homes/u/u650/my_ips_trunk/bin.

Also, checked that 
/global/homes/u/u650/my_ips_trunk/bin/rf_genray_LH.py
(and rf_genray_EC.py) are same as in 
/global/homes/u/u650/my_ips_trunk/components/rf/genray/src



Feb. 10, 2012
=============

prepare_genray_input.f90:
For ps_add_nml.eq.'disabled'
Modified code to NOT set frqncy in genray.in
from the PS, if the PS freq_lh (ec, ic, also) is 0.


Oct. 4, 2013
============
BH: Have recompiled prepare_genray_input.f90 (now in 
prepare_genray_input.f90.bak) using updated *.f files
from genray distribution, updated by inclusion of 
delim='apostrophe' clause in open statements.
However, it looks like a small change has been
introduced in the svn version of prepare_genray_input.f90
744,745c736
< c     eqdskin=trim(ps%eqdsk_file)
<       eqdskin='eqdsk'
---
>       eqdskin=trim(ps%eqdsk_file)

Probably need to recompile with the 
eqdskin='eqdsk' line???


Oct. 17, 2013
=============

Paul Bonoli (ptb) and I (BH) reviewed his changes of prepare_genray_input.f90
and process_genray_output.f90 (none in this file), rf_genray_LH.py,
to get IPS working for a time-dependent simulation of C-Mod using
TSC, genray, and cql3d.   BH took ptb prepare_genray_input.f90, updated
it with the delim='apostrophe' stanza, and recompiled giving 
prepare_genray_input/process_genray_output.

Note: genray also has been updated, changing read_write_genray_input_prep.f
and the include files in ./genray_includes.  genray is up-to-date
in swun svn.

BH svn committed changes to prep/proc genray, rf_genray_LH.py and 
the genray .f and .i files, and this file.
Also update and commit update_prep-proc_genr_distn.sh which get the 
necessary files from the genray svn trunk distribution, and puts them
in ./  and ./genray_includes.
Committed updated Makefile.
That is,
svn commit README_genray_prep-proc_devel Makefile prepare_genray_input.f90 \
read_write_genray_input_prep.f rf_genray_LH.py \
update_prep-proc_genr_distn.sh

Also, put prepare_genray_input/process_genray_output in 
~u650/my_ips_trunk/bin  and  checked there that have
ln -s process_genray_output process_genray_output_mcmd


Oct. 21, 2013
=============

Paul Bonoli and Francesca Poli had a problem with make (under
swim environment), since BH had mucked up the svn commit of the
./genray_includes directory.   This was fixed.

Also, updated update_prep-proc_genr_distn.sh.  This script
can be used to update the genray namelist related source files which
are used by prepare_genray_input.f90. 
The genray namelist related files should be updated whenever
using a new genray executable with changes in these files, to
keep the generated namelists in sync.


Oct. 23, 2013
=============

Changed an inappropriate delim='apostrophe' in prepare_genray_input.f90, 
since in this case, the code is simply transcribing lines, as character 
variables, from one file to the other (using subroutine transcibe).
The delim='apostrophe' clauses are kept in read_write_genray_input_prep.f,
used with prepare_genray_input.f90. Character variable lengths of namelist
variables eqdskin and netcdfnm were reduced from 512 to 256, in keeping
with an apparent ~300 character limit for (pgi) compiler namelist read.

rf_genray_EC_p.py has been updated in ~u650 genray/src area from
svn.  This is a "Programmable version: Batchelor (6/2/2013)".
This version allows the aiming angles and power for multiple launchers to be set and
changed over time from within the simulation config file.

This svn version is what is being tested by Poli, 2013-10-18.


Dec. 29, 2013
=============
! BH131229:  Added subroutine comma_check, which removes any ","
! BH131229:  which appears in columns 1:5 of the updated genray.in
! BH131229:  file.  Such commas appear when a namelist array name
! BH131229:  is of a specfic length giving a "," separated
! BH131229:  from two numbers by a line-feed. This appears as two
! BH131229:  separators (to pgi compiler, at least), giving an
! BH131229:  incorrect array length (as I understand the error).
!
Updated svn.


Dec. 15, 2014
=============
BH updated read_write_genray_input_prep.f and ./genray_includes
include files to genray version="genray_v10.8_141024"
Updated ./Makefile
svn updated ~u650/my_ips_trunk/swim.bashrc.hopper, and sourced it.
svn committed the changed files, including this file.
prepare_genray_input.f90 updated with minor spacing issue.
Remade prepare_genray_input and process_genray_output.