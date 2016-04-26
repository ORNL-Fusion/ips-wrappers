!
!
      program prepare_genray_input

!Invoked as:  prepare_genray_input [ipsmode rfmode isource_string
!                                   genraynml adj_read ps_add_nml]
!             Up to 6 command line arguments, as described below.
!
! BH131229:  Added subroutine comma_check, which removes any ","
! BH131229:  which appears in columns 1:5 of the updated genray.in
! BH131229:  file.  Such commas appear when a namelist array name
! BH131229:  is of a specfic length giving a "," separated
! BH131229:  from two numbers by a line-feed. This appears as two
! BH131229:  separators (to pgi compiler, at least), giving an
! BH131229:  incorrect array length (as I understand the error).
!
! BH131023:  Changed to delim='none', since writes here were not
! BH131023:  namelist writes per se, but rather simple character
! BH131023:  variable writes.  The namelist write is performedin 
! BH131023:  read_write_genray_input_prep.f. and which uses delim=
! BH131023:  'apostrophe'
! BH130927:  Added delim='apostrophe' to namelist writes.
!
! PB130829:  Add coding for LHRF where separate
!   grills do not necessarily correspond to separate sources (ngenray_src=1).
!   In this case we use the values of powers in powers(i) originally specified
!   in the genray.in template file to specifiy a fraction of the total power
!   for each LH grill.  BH extended to other rf types.
!
! BH+PB130328:  Adding coding to avoid divide by zero, when powers
!               of all sources are zero.
!
! BH111110:  Substantial re-work, to accomodate dynamic dimension
! BH111110:  genray modification, and other aspects of a TSC generated
! BH111110:  plasma state file.
!
! PS2, BH080512  (with Doug McCune)

!.......................................................................
!     For ipsmode='step':
!     Prepares new genray.in  namelist file from an input
!     "namelist template" file (gernraynml) specified on the command
!     line. The resulting genray.in file contains the "template" data
!     modified by Plasma State 2 (PS) plasma profiles and PS abridged
!     species list (if ions involved), plus PS source specifications
!     for a particular RF mode (rfmode is also command line input).
!     The PS state file is to be named cur_state.cdf.
!     A list of the PS source specification variables passed to 
!     genray.in will be maintained below.
!     Prepare_genray_input uses several include files from the genray
!     source, plus genray.in variables specified through genraynml,
!     and PS2 variables.
!
!     prepare_genray_input takes up to 6 optional command line
!     arguments:   1st: ipsmode, either 'init' or 'step' for coupling
!                       to SWIM Integrated Plasma Simulator (IPS).
!                       ipsmode='init', only initializes some
!                       component dimensions in the plasma state.
!                       default='step'  
!                  2nd: rfmode ('EC','LH','IC' are possible inputs)
!                       default='EC'
!                  3rd: isource_string=RF Source number.  Points to
!                       the first of one-or-more sources in the
!                       Plasma State for data passed to genray
!                       via the namelist file.
!                       default='1'
!                  4th: genraynml, full path specification of template
!                       genray nml file,   default='./genray.in' (in
!                       which case genray.in will be modified in place)
!                  5th: adj_read ('enabled', then read self-adjoint
!                       CD efficiency function from file, rather 
!                       than calculate it.)
!                       default='disabled'
!                  6th: ps_add_nml, indicates wherether some data,
!                       indicated below, should be added to the plasma
!                       state from genray namelist.
!                       Valid values: 'disabled','enabled','force'.
!                       E.G., Can be used to update a test PS.
!                       default='disabled' [the usual situation].
!
!     To set a given command line argument to non-default value, it
!     is necessary to set the preceeding values in the list.
!
!     rfmode indicates which PS component data will be used.
!     For example, rfmode='EC', electron profiles + zeff, plus
!     additional variables as specified below, will be inserted 
!     in the genray.in file.
!
!     isource_string keeps count of the source numbers for each
!     rfmode.   The aim with the associated process_genray_output.f90
!     code is that if isource_string is greater than 1, 
!     then powers and currents are accumulated in the PS,
!     as required, for given rfmode.
!     If the wave frequencies are equal for all sources for
!     given rfmode, usually all sources can be entered through
!     one genraynml file, and the power deposition and current
!     drive for all sources will be attributed in process_genray_output
!     to the isource_string source.  [Genray can only do one
!     frequency at a time.  Also, if treating several sources,
!     i.e, EC-like cones or LH/FW-like grills, for a given frequency,
!     it accumulates the power deposition and current drive
!     without recording which source is responsible.]
!     For sources with given rfmode but
!     different frequencies, isource_string is expected to have
!     the PS value of the first source in the set for that
!     frequency. (Usually it will be most convenient if sources
!     in the plasma state are grouped by frequency.)  For PS
!     storage which depends on isource_string, the genray data
!     for all sources specified in a given genraynml will be
!     stored in the locations associated with the given isource_string,
!     as if the data were for a single source.
!
!     If separate entries of power deposition and current drive
!     are desired in the PS for each source, then there needs to be
!     separate genraynml for each source, and separate calls for each
!     to prepare_genray_input, xgenray, <process_genray_output.
!     For each rfmode, the isource_string values would be advanced
!     by 1 for each source.
!
!     genraynml is a template of namelist input for the rfmode.
!     Most settings for genray will come from genraynml, setup
!     for a particular rfmode.
!
!     adj_read will pass data to genray indicating if the self-
!     adjoint current drive function should be recalculated
!     (a relatively CPU intensive CD calc, but only need be
!     recalculated for a new eqdsk).
!
!     There will be a check on whether plasma state RF variables
!     are defined.  Normally, they will be.
!     If not, and command line variable PS_add_nml.eq.'enabled',
!     then will add data to the PS from the genraynml data.
!     If RF variables are defined or not, but PS_add_nml.eq.'force',
!     the values from genraynml data will be inserted in the PS.
!     Normally, PS_add_nml will be 'disabled', and genray.in data
!     will be set from the PS data.
!
!
!.......................................................................



!.......................................................................
!     List of RF-related PS variables passed to genray.in, 
!     overwriting data from genraynml (in addition to plasma profiles
!     and the eqdsk file path), depending on rfmode:
!     (If PS_add_nml.ne.'disabled' these values may alternatively be
!     set in the PS using data from genraynml.)
!
!     EC   necrf_src=   number of ECRF sources
!          power_ec(necrf_src)= power on each ECRF source (W)
!          if PS freq_ec .ne. 0:
!          freq_ec(necrf_src)=frequency on each ECRF source (Hz)
!
!     LH   nlhrf_src= number of LHRF sources
!          power_lh(nlhrf_src)= power on each LHRF source (W)
!          if PS freq_lh .ne. 0:
!          freq_lh(nlhrf_src)= frequency on each LHRF source (Hz)
!
!     IC   nicrf_src= number of ICRF sources (includes FW, HHFW)
!          power_ic(nicrf_src)= power on each ICRF source
!          if PS freq_ic .ne. 0:
!          freq_ic(nicrf_src)= frequency on each ICRF source (Hz)
!
!     For PS_add_nml.eq.'enabled', only add appropriate freq_xx
!     to PS, if ps%freq_xx=0, i.e., undefined.
!
!
!          The PS number of sources should be .ge. the number
!          cumulatively specified in genraynml, for a given rfmode.
!          This will be checked.
!
!          There can be one or more calls to prepare_genray_input (and 
!          to genray) for a given rfmode.   If the frequencies are 
!          different, then as mentioned above the sources must be 
!          specified in different genraynml namelist files.
!          The command line variable isource_string specifies the
!          first PS source in the sequence of sources specified
!          in the genray.in file.  The above powers and frequencies
!          are transferred to the corresponding sources (beginning
!          from 1st in each genraynml), up to PS source number
!          (isource_string+ngenray_source), where ngenray_source
!          is either genray namelist variable ncone or ngrill
!          depending on whether input sources are ray cones 
!          (genraynml istart=1, as in EC cases) or grills (istart=2,3,
!          as in LH,FW,EBW) cases.
!          The maximum of (isource_string+ngenray_source) should
!          agree with the number of sources specified in the PS,
!          for given rfmode.
!
!          [The genraynml file will specify other aspects of the
!          source antennas, plus additional relevant input.
!          If additional PS variables not presently in the PS are 
!          to be passed to the genray namelist, then add these to the
!          above list and modify code below in an analogous manner.]
!
!.......................................................................
!     List of genray namelist variables modified by prepare_genray_input:
!     [This list should be maintained, if there are additions.]
!     
!     partner, idens(=1), ndens(=201), nbulk(PS abridged list plus 1,
!     if .gt.1), nonuniform_profile_mesh(='disabled'),
!     [cBH111108: temporarily, see below, nbulk=nspec_tha+1] 
!     den_scale()(=1.), temp_scale()(=1.), eqdskin, indexrho(=2),
!     i_calculate_or_read_adj_function, frqncy, powtot(), powers(),
!     (density,temperature-profiles, species=1,nbulk), (zeff-profiles),
!     i_calculate_or_read_adj_function
!.......................................................................
!
!
!
!.......................................................................
!
!     Some files are used from the genray distribution:
!     =================================================
!      -the subset of include files in ./genray_includes
!      -bcast.f
!      -read_write_genray_input_prep.f [This file is a modified version
!                                  of read_write_genray_input_prep.f, 
!                                  by commenting out the if clause which
!                                  includes 'call read_transport_prof'.]
!     Using these genray distn files is aimed at keeping genray and 
!     prepare_genray_input.f90 in sync.   When any of these files change 
!     in the genray distribution, the files should be updated in the
!     present code suite.
!
!     The subroutine prepare_genray_input in read_write_genray_input_prep.f
!     is not used for prepare_genray_input.f90, but is part of a test
!     code in the genray distribution, and provides a template for
!     reading and writing namelists in prepare_genray_input.f90.
!
!.......................................................................
!
!     Compile and load line:
!     ======================
!     Include files in subdir ./genray_includes:
!                             param.i adj_nml.i ETC
!
!viz.pppl.gov compile:
! ifort -c -I ./genray_includes read_write_genray_input_prep.f bcast.f
! ifort -fixed -80 -c -O -w -mp -fpe0 -heap-arrays 20 -I./genray_includes -I$NTCCHOME/mod -I$NETCDFHOME/include -o prepare_genray_input.o  prepare_genray_input.f90
!
!Check which modules are loaded with 'module list'.

!viz.pppl.gov load: (with module ntcc/intel.10.1_mkl.10)
! ifort -fpe0 -o prepare_genray_input    prepare_genray_input.o read_write_genray_input_prep.o bcast.o -L$NTCCHOME/lib -lplasma_state -lps_xplasma2 -lplasma_state_kernel -lxplasma2 -lgeqdsk_mds -lmdstransp -lvaxonly -lnscrunch -lfluxav -lr8bloat -lpspline -lezcdf -llsode -llsode_linpack -lcomput -lportlib -L$MDSPLUS_DIR/lib -L$GLOBUS_LOCATION/lib -lMdsLib -lglobus_xio_gcc64  -L$MKLHOME/lib/64 -lmkl_lapack -lmkl_ipf -lmkl_sequential -lmkl_core -lguide -lpthread -L$NETCDFHOME/lib -lnetcdf

!viz.pppl.gov load: (with module ntcc/intel.10.0_mkl.9)
! ifort -fpe0 -o prepare_genray_input    prepare_genray_input.o read_write_genray_input_prep.o  bcast.o -L$NTCCHOME/lib -lplasma_state -lps_xplasma2 -lplasma_state_kernel -lxplasma2 -lgeqdsk_mds -lmdstransp -lvaxonly -lnscrunch -lfluxav -lr8bloat -lpspline -lezcdf -llsode -llsode_linpack -lcomput -lportlib -L$MDSPLUS_DIR/lib -L$GLOBUS_LOCATION/lib -lMdsLib -lglobus_xio_gcc64  -L$MKLHOME/lib/64 -lmkl_lapack -lmkl_ipf -lpthread  -lguide -lpthread -L$NETCDFHOME/lib -lnetcdf

!.......................................................................


!.......................................................................
! The following USE gives access to the extensive plasma state data.
! (It is necessary to rename "output" in the plasma state to avoid
! conflict with genray namelist /output/.)

      USE plasma_state_mod, ps_output_local => output

!.......................................................................


      USE swim_global_data_mod, only :
     1 rspec, ispec,     ! kind specification for real (64-bit)
                         ! and integer (32-bit).
     1 swim_error,        ! error routine
     1 swim_filename_length   ! length of filename string
    
!      implicit none
      implicit integer (i-n), real*8 (a-h,o-z)

!.......................................................................
!     NAMELISTS AND THEIR STORAGE, from GENRAY FILES
!.......................................................................

      include 'param.i'
c----------------------------------------
c     type declaration and storage of all namelists variables
c-----------------------------------
      include 'cone_nml.i'
      include 'emissa_nml.i'
      include 'grill_nml.i'
      include 'ions_nml.i'
      include 'one_nml.i'
      include 'onetwo_nml.i'
      include 'output_nml.i'
      include 'scatnper_nml.i'    
      include 'six_nml.i'  !dens1(ndensa,nbulka),temp1(ndensa,nbulka)
                           !tpop1(ndensa,nbulka),vflow1(ndensa,nbulka)
                           !zeff1(ndensa)   
      include 'adj_nml.i'
c--------------------------------------------------------------
c     declaration of namelist variables which
c     are used in subroutine dinit_mr and
c     are not used in any common block
c--------------------------------------------------------------
      include 'dinit_nml.i' 
c--------------------------------------------------------------
      include 'rkutta.i' ! contains common /rkutta/ prmt(9)
c--------------------------------------------------------------
c     list of all namelists
c--------------------------------------------------------------
      include 'name.i'    ! list of all namelists
c--------------------------------------------------------------

      integer ndim,nray
      integer i,j,k

      save ndim  ! to give ndim to
                 ! subroutine write_all_namelists

      character*8 ipsmode,rfmode,isource_string,adj_read,ps_add_nml
cBH131022      character(swim_filename_length) genraynml
      character*256 genraynml
      
      integer isource,iargs
cBH131022      character(swim_filename_length) eqdsk_name
      character*256 eqdsk_name
      integer ngenray_source

      logical lfrq

c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus
      character filenc*128,ltitle*128,name*128
      integer dims_id(2),start(2),count(2),countm(2),char256dim
      integer nprofs
      integer nj_id,nspecgr_id
      integer nprofs_id,nbulk_id

c---------------------------------------------------------------
c     The _nc declarations are for organizing the data to be
c     inserted in the genray_profs_in.nc netcdf file.
c---------------------------------------------------------------

      integer, parameter :: nraa=201,nspec_maxa=12
      integer nspecgr
                    !nspecgr is the number of species passed to genray,
                    !including electrons, ions, possible hot plasma
                    !species.  Must be .le.nspec_maxa parameter above.
      integer nj,njp1
      REAL(KIND=rspec) dr_nc
      REAL(KIND=rspec) charge_nc(nspec_maxa),dmas_nc(nspec_maxa),
     +     r_nc(nraa),en_nc(nraa,nspec_maxa),temp_nc(nraa,nspec_maxa),
     +     zeff_nc(nraa)
      integer, dimension(0:nspec_maxa) :: incl_species=0 !indicator for
                                          !species to send to genray.
      common/pass_to_write_all_namelists/
     +     r_nc,en_nc,temp_nc,zeff_nc

c---------------------------------------------------------------
!  Unified plasma parameter arrays (on ps%nspec_alla list)
!    abridged thermal species + all fast species combined.

      REAL(KIND=rspec), dimension(:), allocatable :: rho_bdy_rezon
      REAL(KIND=rspec), dimension(:), allocatable :: dense,tempe,zeff
      REAL(KIND=rspec), dimension(:,:), allocatable :: densi,tempi
      REAL(KIND=rspec), dimension(:), allocatable :: wk_eperp,wk_epll,
     +                                               wk_dens
c---------------------------------------------------------------
      integer :: cclist(ps_ccount)  ! State partial copy controls.
                                    ! Ps_ccount is number of known
                                    ! components in the PS.
c-----------------------------------------------------------------------
!  Setup some allocatable variables from PS, which will be
!  checked if allocated:
      REAL(KIND=rspec),dimension(:),allocatable :: nbeami,nmini,nfusi

c-----------------------------------------------------------------------

      data start/1,1/

c---------------------------------------------------------------

c--------------------------------------------------------------
c     Read command line arguments
c--------------------------------------------------------------

c     Defaults:
      ipsmode='step'
      rfmode='EC'
      isource_string='1'
      genraynml='./genray.nc'
      adj_read='disabled'
      ps_add_nml='disabled'

c     F2003-syntax: command_argument_count()/get_command_argument(,)
      iargs=command_argument_count()
      write(*,*)'iargs=',iargs
      if (iargs.ne.6) then
         write(*,*)'prepare_genray_input usage: '
         write(*,*)'Up to six command line arguments, '
         write(*,*)
     +   'ipsmode rfmode genraynml isource_string adj_read ps_add_nml'
     +   //' (refer to code)'
      endif

      if (iargs.ge.1)   call get_command_argument(1,ipsmode)
      if (iargs.ge.2)   call get_command_argument(2,rfmode)
      if (iargs.ge.3)   call get_command_argument(3,isource_string)
      if (iargs.ge.4)   call get_command_argument(4,genraynml)
      if (iargs.ge.6)   call get_command_argument(5,adj_read)
      if (iargs.ge.6)   call get_command_argument(6,ps_add_nml)

      ipsmode=trim(ipsmode)
      rfmode=trim(rfmode)
      isource_string=trim(isource_string)
      genraynml=trim(genraynml)
      adj_read=trim(adj_read)
      ps_add_nml=trim(ps_add_nml)

      write(*,*)'prepare_genray_input command line arguments: ',
     +  ipsmode,'  ',rfmode,isource_string,trim(genraynml),'  ',
     +  adj_read,'  ',ps_add_nml

      read(isource_string, '(i4)') isource


c-----------------------------------------------------------------------
c     Copy genraynml to ./genray.in
c     Transcribe file, since not all fortrans have access
c     to shell commands (for which simple cp would work).
c-----------------------------------------------------------------------

      call transcribe(trim(genraynml),'./genray.in')


c-----------------------------------------------------------------------
c     the creation of default input data like in
c     genray.dat file format (mixed units, not all MKSA).
c     See genray.in_template file in the genray distribution for
c     description of the namelist variables in MKSA.
c-----------------------------------------------------------------------
      call default_in
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c     Input MKSA file genray.in set up above.  Following call will
c     change default data to genray.in file (MKSA) format.
c-----------------------------------------------------------------------

      call transform_input_data_to_MKSA

c-----------------------------------------------------------------------
c     Read all namelists from input genraynml MKSA file.
c     (genray_in_dat is o/p name of nml file: 'genray.in' here)
c-----------------------------------------------------------------------

      call read_all_namelists(genray_in_dat,ndim,nray)
      write(*,*)'after read_all_namelists'


!.......................................................................
!     Make available all variables in the plasma state module, specified
!     in swim_state_spec.dat (or corresponding .txt) in components/
!     state/src/plasma_state.  See also McCune's Plasma_State_Vnn
!     design doc at http://www.cswim.org/componentdescrip/.
!     These variable are available with the prepended "component
!     selector" for the current state, ps%, or for the previous state,
!     %psp, and a few others.
!.......................................................................

      call ps_get_plasma_state(ierr,filename='./cur_state.cdf')
      if(ierr .ne. 0) stop 'cannot get plasma state needed for profiles'


!.......................................................................
!     If ipsmode.eq.'init', then set approrpriate dimensions in the 
!     plasma state component section, dependent on rfmode
!.......................................................................

      if (ipsmode.eq.'init') then
         if (rfmode.eq.'EC') then
            
            nrho=NR  !NR is number of bin boundaries spec'd in genray.in

            write(*,*)'nrho frm nml, ps%nrho_ecrf frm PS =',
     +	               nrho,ps%nrho_ecrf

!           Necessary PS dimensions (user has to set his component dims):
            if(ps%nrho_ecrf.eq.nrho) then
               continue         ! dimension OK already
            else if(ps%nrho_ecrf.eq.0) then
               ps%nrho_ecrf=nrho
               call ps_alloc_plasma_state(ierr) !set these PS dims in PS
               write(*,*)'ps%nrho_ecrf set= ',nrho
c           Following if clause should not happen in ordinary IPS usage
c           (Can occur if using already initialized PS)
            else if(ps%nrho_ecrf.ne.nrho) then
               write(*,*) 
     +              ' * prepare_genray_input: reset EC profile size'
               write(*,*) ' * from ',ps%nrho_ecrf,' to ',nrho
               
!              copy all EXCEPT EC component profiles
               cclist = ps_cc_all
               cclist(ps_cnum_EC)=0

               call ps_copy_plasma_state(ps, aux, ierr, cclist = cclist)
               if(ierr.eq.0) then
                  ! OK, copy back to ps
                  call ps_copy_plasma_state(aux, ps, ierr)
                  
                  if(ierr.eq.0) then
                     ! set desired dimension
                     ps%nrho_ecrf=nrho
                     call ps_alloc_plasma_state(ierr) 
                                               !set these PS dims in PS
                  endif
               endif  !on ierr
            endif
            
            if (ierr.ne.0) then
               write(*,*)'prepare_genray_input:  EC'//
     +                     ' ps_alloc_plasma_state: ierr=',ierr
               stop
            endif

c           Initialize arrays, 
c           ps%rho_ecrf to be same as rho_bin as in genray.
c           A LITTLE AWKWARD, in that if the algrithm in genray
c           is changed from simple equi-spaced, then must transmit
c           the coding here also.
            hrho=1.d0/(NR-1)
            do i=1,NR
cgenray        rho_bin(i)=hrho*(i-1)
               ps%rho_ecrf(i)=hrho*(i-1)
            enddo

            ps%peech=0.0_rspec
            ps%curech=0.0_rspec

         elseif(rfmode.eq.'LH') then

            nrho=NR  !NR is number of bin boundaries spec'd in genray.in

            write(*,*)'nrho frm nml, ps%nrho_lhrf frm PS =',
     +	               nrho,ps%nrho_lhrf

!           Necessary PS dimensions (user has to set his component dims):
            if(ps%nrho_lhrf.eq.nrho) then
               continue         ! dimension OK already
            else if(ps%nrho_lhrf.eq.0) then
               ps%nrho_lhrf=nrho
               call ps_alloc_plasma_state(ierr) !set these PS dims in PS
               write(*,*)'ps%nrho_lhrf set= ',nrho
c           Following if clause should not be satisfied in ordinary IPS
c           usage (Can occur if using already initialized PS)
            else if(ps%nrho_lhrf.ne.nrho) then
               write(*,*) 
     +              ' * prepare_genray_input: reset LH profile size'
               write(*,*) ' * from ',ps%nrho_lhrf,' to ',nrho
               
!              copy all EXCEPT LH component profiles
               cclist = ps_cc_all
               cclist(ps_cnum_LH)=0

               call ps_copy_plasma_state(ps, aux, ierr, cclist = cclist)
               if(ierr.eq.0) then
                  ! OK, copy back to ps
                  call ps_copy_plasma_state(aux, ps, ierr)
                  
                  if(ierr.eq.0) then
                     ! set desired dimension
                     ps%nrho_lhrf=nrho
                     call ps_alloc_plasma_state(ierr)
                                               !set these PS dims in PS
                  endif
               endif
            endif
            
            if (ierr.ne.0) then
               write(*,*)'prepare_genray_input:  LH'//
     +                     ' ps_alloc_plasma_state: ierr=',ierr
               stop
            endif

c           Initialize arrays, 
c           ps%rho_lhrf to be same as rho_bin as in genray.
c           A LITTLE AWKWARD, in that if the algrithm in genray
c           is changed from simple equi-spaced, then must transmit
c           the coding here also.
            hrho=1.d0/(NR-1)
            do i=1,NR
cgenray        rho_bin(i)=hrho*(i-1)
               ps%rho_lhrf(i)=hrho*(i-1)
            enddo

            ps%pelh=0.0_rspec
            ps%pilh=0.0_rspec
            ps%curlh=0.0_rspec
	    ps%pelh_src=0.0_rspec
	    ps%pilh_src=0.0_rspec
	    curlh_src=0.0_rspec

         elseif(rfmode.eq.'IC') then

            nrho=NR  !NR is number of bin boundaries spec'd in genray.in

            write(*,*)'nrho frm nml, ps%nrho_icrf frm PS =',
     +	               nrho,ps%nrho_icrf

!           Necessary PS dimensions (user has to set his component dims):
            if(ps%nrho_icrf.eq.nrho) then
               continue         ! dimension OK already
            else if(ps%nrho_icrf.eq.0) then
               ps%nrho_icrf=nrho
               call ps_alloc_plasma_state(ierr) !set these PS dims in PS
               write(*,*)'ps%nrho_icrf set= ',nrho
c           Following if clause should not happen in ordinary IPS usage
c           (Can occur if using already initialized PS)
            else if(ps%nrho_icrf.ne.nrho) then
               write(*,*) 
     +              ' * prepare_genray_input: reset EC profile size'
               write(*,*) ' * from ',ps%nrho_icrf,' to ',nrho
               
!              copy all EXCEPT IC component profiles
               cclist = ps_cc_all
               cclist(ps_cnum_IC)=0

               call ps_copy_plasma_state(ps, aux, ierr, cclist = cclist)
               if(ierr.eq.0) then
                  ! OK, copy back to ps
                  call ps_copy_plasma_state(aux, ps, ierr)
                  
                  if(ierr.eq.0) then
                     ! set desired dimension
                     ps%nrho_icrf=nrho
                     call ps_alloc_plasma_state(ierr)
                                                !set these PS dims in PS
                  endif
               endif
            endif
            
            if (ierr.ne.0) then
               write(*,*)'prepare_genray_input:  EC'//
     +                     ' ps_alloc_plasma_state: ierr=',ierr
               stop
            endif

c           Initialize arrays, 
c           ps%rho_icrf to be same as rho_bin as in genray.
c           A LITTLE AWKWARD, in that if the algrithm in genray
c           is changed from simple equi-spaced, then must transmit
c           the coding here also.
            hrho=1.d0/(NR-1)
            do i=1,NR
cgenray        rho_bin(i)=hrho*(i-1)
               ps%rho_icrf(i)=hrho*(i-1)
            enddo

            ps%picrf_totals=0.0_rspec
            ps%picth=0.0_rspec
            ps%curich=0.0_rspec
            ps%cdicrf=0.0_rspec
            ps%curich=0.0_rspec

         endif  ! on rfmode

!     Stores the initialized plasma state--doesn't commit it. 
!     The previous one is still around

        call ps_store_plasma_state(ierr)
        if(ierr .ne. 0) then
          write(iout,*)
     +    'Cannot ps_store_plasma_state in prepare_genray_input:init'
        else
	  write(iout,*) "init'd cur_state.cdf stored"
        end if


        go to 999  ! go to end
                        
      endif  ! on ipsmode.eq.'init'

!.......................................................................
c     if only need electron profiles plus zeff (nbulk=1, from
c     genray namelist), then set nspecgr=1.  Else, obtain number
c     of ions for the plasma state abridged species list.
!.......................................................................
      if (nbulk.eq.1) then 
         nbulk=1       ! i.e., no change from genray.dat/.in input value
         nspecgr=1
      else             ! Get number of ion species from PS abridged list
!BH111108         nbulk=ps%nspec_alla+1
         nbulk=ps%nspec_tha+1
         nspecgr=nbulk
         if (nbulk.gt.nbulka) then
            write(*,*)
	    write(*,*)'nbulk from PS=',nbulk,' nbulka=',nbulka
            write(*,*)'STOP:  Recompile genray with larger nbulka?'
            write(*,*)
            stop
         endif
      endif

c---------------------------------------------------------------
c     Check/adjust plasma profile specs for passage from PS to
c     genray.in
c---------------------------------------------------------------

      if (idens.ne.1.or.ndens.ne.201) then
         write(*,*)'Resetting idens for uniform, 201-pt plasma profs'
         idens=1
         ndens=201
         nonuniform_profile_mesh='disabled'
      endif
!      BH111110:  genray mostly dynamically dimensioned:
!      nj=min(nra,nraa)  ! nra is from genray param.i
      nj=ndens
      njp1=nj+1
!BH111110      if (nj.ne.201) then
!BH111110         write(*,*)
!BH111110         write(*,*)'Prefer genray compiled with nra.ge.201'
!BH111110         STOP   'Clean this dimension coding up'
!BH111110         write(*,*)
!BH111110      endif


c     Check/adjust izeff
      if (izeff.ne.2) then
         write(*,*)'Resetting izeff=2 for all nbulk cases'
         izeff=2
      endif

c     Check/adjust partner
      if (partner.ne.'genray_profs_out.nc') then
         write(*,*)
     +        'Resetting partner=genray_profs_out.nc for PS coupling'
         partner='genray_profs_out.nc'
      endif

c     Make sure density and temp scaling turned off
      do i=1,nbulk
         den_scale(i)=1.d0
         temp_scale(i)=1.d0 
      enddo

c     Reset i_calculate_or_read_adj_function, per command line variable
      if (adj_read.eq.'enabled') then
         i_calculate_or_read_adj_function=1
      else
         i_calculate_or_read_adj_function=0
      endif

         


c     Reset genray.in eqdskin according to PS
c From Francesca Poli, 131018: we had to change the [eqdskin] line 
c to make sure genray reads the correct eqdsk file.
c I was getting an error, because GENRAY was looking for a file 
c that was not there.
c The ps_eqdsk variable had in memory the name of the default plasma
c state file, which starts with $INPUT_SUFFIX, rather than the correct
c plasma state eqdsk file that start with $SIM_NAME.
c      eqdskin=trim(ps%eqdsk_file)
      eqdskin='eqdsk'
c     Reset genray.in indexrho so radial coord proportional
c     to sqrt(tor flux)
      indexrho=2


c     Number or sources in genraynml:
cBH120829:  Will be treated as 1, for each genraynml, if cannot be 
cBH120829:  modify the PS, that is, if PS_add_nml.eq.'disabled'.
cBH120829:  In this case, if ncone or ngrill .gt.1, they are treated as
cBH120829:  a single source. Might need to re-examine this logic in 
cBH120829:  isource.gt.1 cases.
      if (PS_add_nml.eq.'disabled') then
         ngenray_source=1
      else
         if (istart.eq.1) then     !cone antenna
            ngenray_source=ncone
         elseif (istart.eq.2 .or. istart.eq.3) then !grill antenna
            ngenray_source=ngrill
         else
            write(*,*)'Problem with istart in genraynml!:   STOP'
            STOP
         endif
      endif
      

!.......................................................................
!     Usual operation: PS_add_nml will be 'disabled', and some genray.in
!     data will be set from the PS data.   The number of RF sources
!     in genraynml will be less than or equal to PS number of sources.
!
!     However, prepare_genray_input also treats situations where 
!     specification of RF sources may not be set in the PS (number
!     of sources = 0), or it is desired to reset the PS RF variables
!     based on genraynml data.
!     Check whether plasma state RF variables are defined.
!     If not, and command line variable PS_add_nml.eq.'enabled',
!     then will add data to the PS from the genraynml data.
!     If PS RF variables are defined, but are incompatible with
!     with genraynml data, then code will stop, else genraynml source
!     data will be inserted in the PS.
!     If RF variables are defined or not, but PS_add_nml.eq.'force',
!     the values from genraynml data will be inserted in the PS.
!.......................................................................


      if (PS_add_nml.eq.'disabled') then  !PS data inserted in genray.in

!.......................................................................
!     Check isource+(sources specified in genraynml) does not exceed
!     number of sources in the PS, and set genraynml powers and 
!     frequencies from the PS.
!.......................................................................



c     Usually rfmode='EC' will be associated with istart=1 (launch
c     cone of rays from outside the plasma), but if modeling EBW
c     waves, then usually use istart=2 or 3.
c     Usually istart=2 for LH and IC (FW or HHFW).
c     However, all istart possibilities are enabled for each rfmode.

        if (rfmode.eq.'EC') then
            
            if ((isource+ngenray_source-1).gt.ps%necrf_src) then
               write(*,*)'Incompatible source number specs'
	       write(*,*)'ps%necrf_src=',ps%necrf_src
               write(*,*)'(isource+ngenray_source-1)=',
     +                    (isource+ngenray_source-1)
               write(*,*) 'STOP'
               STOP
            endif
            
c     Set genray.in frqncy from PS
            lfrq=ps%freq_ec(1).ne.0   !test if ps freq set
            if (lfrq) frqncy=ps%freq_ec(1)

            if (istart.eq.1) then
               do is=1,ngenray_source
                  if(ncone .eq. 1) then       !PTB idea, as below for LH
                     powtot(is)=ps%power_ec(isource+(is-1))
                  endif
                  if (lfrq.and.frqncy.ne.ps%freq_ec(isource+(is-1)))then
                     write(*,*)'EC freqs must be = for each genray call'
                     write(*,*) 'STOP'
                     STOP
                  endif
                  if(ncone .gt. 1) then
                     sumgrillp = 0.0_rspec
                     do ig=1,ncone
                        sumgrillp = sumgrillp + powtot(ig)
                     enddo
                     if (sumgrillp.gt.0.0_rspec) then
                        do ig=1,ncone
                           powtot(ig) = (powtot(ig)/sumgrillp)*
     +                                  ps%power_ec(isource)
                        enddo
                     else
                        do ig=1,ncone
                           powtot(ig) = 0.0_rspec
                        enddo
                     endif
                  endif  !On ncone
               enddo  !On is
            elseif(istart.eq.2 .or. istart.eq.3) then
               do is=1,ngenray_source
                  if(ngrill .eq. 1) then
                     powers(is)=ps%power_ec(isource+(is-1))
                  endif
                  if (lfrq.and.frqncy.ne.ps%freq_ec(isource+(is-1)))then
                     write(*,*)'freqs must be same for each genray call'
                     write(*,*) 'STOP'
                     STOP
                  endif
                  if(ngrill .gt. 1) then
                     sumgrillp = 0.0_rspec
                     do ig=1,ngrill
                        sumgrillp = sumgrillp + powers(ig)
                     enddo
                     if (sumgrillp.gt.0.0_rspec) then
                        do ig=1,ngrill
                           powers(ig) = (powers(ig)/sumgrillp)*
     +                                   ps%power_ec(isource)
                        enddo
                     else
                        do ig=1,ngrill
                           powers(ig) = 0.0_rspec
                        enddo
                     endif
                  endif  !On ngrill
               enddo  !On is
            endif
            
         elseif (rfmode.eq.'LH') then
            if ((isource+ngenray_source-1).gt.ps%nlhrf_src) then
               write(*,*)'Incompatible source number specs'
               write(*,*) 'STOP'
               STOP
            endif
            
            lfrq=ps%freq_lh(1).ne.0   !test if ps freq set
            if (lfrq) frqncy=ps%freq_lh(1)
            
            if (istart.eq.1) then
               do is=1,ngenray_source
                  if(ncone .eq. 1) then
                     powtot(is)=ps%power_lh(isource+(is-1))
                  endif
                  if (lfrq.and.frqncy.ne.ps%freq_lh(isource+(is-1)))then
                     write(*,*)'freqs must be same for each genray call'
                     write(*,*) 'STOP'
                     STOP
                  endif
                  if(ncone .gt. 1) then
                     sumgrillp = 0.0_rspec
                     do ig=1,ncone
                        sumgrillp = sumgrillp + powtot(ig)
                     enddo
                     if (sumgrillp.gt.0.0_rspec) then
                        do ig=1,ncone
                           powtot(ig) = (powtot(ig)/sumgrillp)*
     +                                  ps%power_lh(isource)
                        enddo
                     else
                        do ig=1,ncone
                           powtot(ig) = 0.0_rspec
                        enddo
                     endif
                  endif  !On ncone
               enddo  !On is
            elseif(istart.eq.2 .or. istart.eq.3) then
               do is=1,ngenray_source
                  if(ngrill .eq. 1) then              ! PTB  08292012   
                     powers(is)=ps%power_lh(isource+(is-1))
                  endif                               ! PTB  08292012
                  if (lfrq.and.frqncy.ne.ps%freq_lh(isource+(is-1)))then
                     write(*,*)'freqs must be same for each genray call'
                     write(*,*) 'STOP'
                     STOP
                  endif
! PTB begins - 08292012
! PTB adds this to deal with the special situation for LHRF where separate
! grills do not necessarily correspond to separate sources (ngenray_src=1).
! In this case we use the values of powers in powers(i) originally specified
! in the genray.in template file to specifiy a fraction of the total power
! for each LH grill, that is presumed to be constant in time. If the total
! LH source power (ps%powerl_lh) in the PS changes in time, then the power
! in each grill is specified as the original fraction of the total source power.
! BH extended idea to other relevant sections.
                  if(ngrill .gt. 1) then
                     sumgrillp = 0.0_rspec
                     do ig=1,ngrill
                        sumgrillp = sumgrillp + powers(ig)
                     enddo
                     if (sumgrillp.gt.0.0_rspec) then
                        do ig=1,ngrill
                           powers(ig) = (powers(ig)/sumgrillp)*
     +                                  ps%power_lh(isource)
                        enddo
                     else
                        do ig=1,ngrill
                           powers(ig) = 0.0_rspec
                        enddo
                     endif
                  endif  !On ngrill
               enddo  !On is
            endif  !on istart
         elseif (rfmode.eq.'IC') then
            if ((isource+ngenray_source-1).gt.ps%nicrf_src) then
               write(*,*)'Incompatible source number specs'
               write(*,*) 'STOP'
               STOP
            endif
            
            lfrq=ps%freq_ic(1).ne.0   !test if ps freq set
            if (lfrq) frqncy=ps%freq_ic(1)
            
            if (istart.eq.1) then
               do is=1,ngenray_source
                  if(ncone .eq. 1) then
                     powtot(is)=ps%power_ic(isource+(is-1))
                  endif
                  if (lfrq.and.frqncy.ne.ps%freq_ic(isource+(is-1)))then
                     write(*,*)'freqs must be same for each genray call'
                     write(*,*) 'STOP'
                     STOP
                  endif
                  if(ncone .gt. 1) then
                     sumgrillp = 0.0_rspec
                     do ig=1,ncone
                        sumgrillp = sumgrillp + powtot(ig)
                     enddo
                     if (sumgrillp.gt.0.0_rspec) then
                        do ig=1,ncone
                           powtot(ig) = (powtot(ig)/sumgrillp)*
     +                                  ps%power_ic(isource)
                        enddo
                     else
                        do ig=1,ncone
                           powtot(ig) = 0.0_rspec
                        enddo
                     endif
                  endif  !On ncone
               enddo  !On is
            elseif(istart.eq.2 .or. istart.eq.3) then
               do is=1,ngenray_source
                  if(ngrill .eq. 1) then
                     powers(is)=ps%power_ic(isource+(is-1))
                  endif
                  if (lfrq.and.frqncy.ne.ps%freq_ic(isource+(is-1)))then
                     write(*,*)'freqs must be same for each genray call'
                     write(*,*) 'STOP'
                     STOP
                  endif
                  if(ngrill .gt. 1) then
                     sumgrillp = 0.0_rspec
                     do ig=1,ngrill
                        sumgrillp = sumgrillp + powers(ig)
                     enddo
                     if (sumgrillp.gt.0.0_rspec) then
                        do ig=1,ngrill
                           powers(ig) = (powers(ig)/sumgrillp)*
     +                                  ps%power_ic(isource)
                        enddo
                     else
                     	do ig=1,ngrill
                           powers(ig) = 0.0_rspec
                        enddo
                     endif
                  endif  !On ngrill
               enddo  !On is
            endif  !On istart
         else
            write(*,*)'Incorrect command line rfmode spec'
            write(*,*) 'STOP'
            STOP
         endif  !On rfmode


cMKSA issue:
         if (rfmode.eq.'EC') then
            frqncy=frqncy !In genray.in mode, genray expects Hz
                                      !PS provides Hz for EC
         elseif (rfmode.eq.'LH') then
            frqncy=frqncy !In genray.in mode, genray expects Hz
                                      !PS provides Hz for LH

         elseif (rfmode.eq.'IC') then
            frqncy=frqncy !In genray.in mode, genray expects Hz
                                      !PS provides Hz for IC

         endif
         
c        write(*,*)'frqncy(Hz), from PS',frqncy


         
      endif  !On PS_add_nml.eq.'disabled'
      


!.......................................................................
!     Usual operation: PS_add_nml will be 'disabled', and genray.in data
!     will be set from the PS data, as above.
!     However,
!     check whether plasma state RF variables are defined.
!     If not, and command line variable PS_add_nml.eq.'enabled',
!     then will add data to the PS from the genraynml data.
!     If RF variables are defined or not, but PS_add_nml.eq.'force',
!     the values from genraynml data will be inserted in the PS.
!.......................................................................

      if (PS_add_nml.eq.'enabled' .or. PS_add_nml.eq.'force') then
                                      !Some genraynml data added to PS

         i_check_rf=0           !Set =1 if PS and genraynml agree
                                !on number of sources.
         
         if (rfmode.eq.'EC') then
            nxxrf_src=ps%necrf_src
            if (ngenray_source.eq.nxxrf_src) i_check_rf=1
            
         elseif (rfmode.eq.'LH') then
            nxxrf_src=ps%nlhrf_src
            if (ngenray_source.eq.nxxrf_src) i_check_rf=1
            
         elseif (rfmode.eq.'IC') then
            nxxrf_src=ps%nicrf_src
            if (ngenray_source.eq.nxxrf_src) i_check_rf=1
            
         else
            write(*,*)'Incorrect command line rfmode spec'
            write(*,*) 'STOP'
            STOP
         endif

         write(*,*)'i_check_rf=',i_check_rf
         if (i_check_rf.eq.0 .and. PS_add_nml.eq.'enabled') then
            write(*,*)'prepare_genray_input: genraynml and PS disagree'
            write(*,*)'                      on number of rfsources, '
            write(*,*)'                      and ps_add_nml.eq.enabled'
            write(*,*)'                      ngenray_source=',ngenray_source
            write(*,*)'                      ps%nxxrf_src=',nxxrf_src
            write(*,*)'                      Need ps_add_nml.eq.force'
            write(*,*) 'STOP'
            STOP
         endif

c     For genray.in and PS agree on number of sources, but for
c     frequency not set, then put in in PS from genray.in.
         if (i_check_rf.eq.1  .and. PS_add_nml.eq.'enabled')  then
            if (rfmode.eq.'EC') then
               do is=1,ngenray_source
                   if(ps%freq_ec(isource+(is-1)).eq.0.0_rspec)
     +                ps%freq_ec(isource+(is-1))=frqncy
               enddo
            elseif (rfmode.eq.'LH') then
               do is=1,ngenray_source
                   if(ps%freq_lh(isource+(is-1)).eq.0.0_rspec)
     +                ps%freq_lh(isource+(is-1))=frqncy
               enddo
            elseif (rfmode.eq.'IC') then
               do is=1,ngenray_source
                   if(ps%freq_ic(isource+(is-1)).eq.0.0_rspec)
     +                ps%freq_ic(isource+(is-1))=frqncy
               enddo
            endif
         endif


c     For number of PS sources .ne. genraynml sources:
c     if  PS_add_nml.eq.'force',  modify PS to agree, and set PS power
c     freq_xx, number of sources, and source names.

         if (i_check_rf.eq.0 .and. PS_add_nml.eq.'force') then
            
            if (rfmode.eq.'EC') then
               if (ps%necrf_src.eq.0) then
                  ps%necrf_src=ngenray_source
                  call ps_alloc_plasma_state(ierr) !set PS dims for rf
               else
                  write(*,*) 
     +              'prepare_genray_input: resetting PS EC source count'
                  write(*,*) ' from ',ps%necrf_src,' to ',ngenray_source
                  
!     copy all EXCEPT EC component profiles
                  cclist = ps_cc_all
                  cclist(ps_cnum_EC)=0
                  
                  call ps_copy_plasma_state(ps,aux,ierr,cclist= cclist)
                  if(ierr.eq.0) then
!     OK, copy back to ps
                     call ps_copy_plasma_state(aux, ps, ierr)
                     
!                     if(ierr.eq.0) then
!     set desired dimension
                        ps%necrf_src=ngenray_source
                        call ps_alloc_plasma_state(ierr) !set dims in PS
!                     endif      !on ierr
                  endif         !on ierr
                  
               endif            !on ps%necrf_src

c     Set EC PS values from genraynml[NOTE: only 1 frqncy]

!     Reminder of data to be passed to PS:
!     EC   necrf_src=   number of ECRF sources
!          power_ec(necrf_src)= power on each ECRF source (W)
!          freq_ec(necrf_src)=frequency on each ECRF source (Hz)
!
!     LH   nlhrf_src= number of LHRF sources
!          power_lh(nlhrf_src)= power on each LHRF source (W)
!          freq_lh(nlhrf_src)= frequency on each LHRF source (Hz)
!
!     IC   nicrf_src= number of ICRF sources (includes FW, HHFW)
!          power_ic(nicrf_src)= power on each ICRF source
!          freq_ic(nicrf_src)= frequency on each ICRF source (Hz)
!
               do i=1,ngenray_source
                  ps%ecrf_src_name=achar(48+i) !achar(48)='0'
                  if (istart.eq.1) ps%power_ec(i)=powtot(i)
                  if (istart.gt.1) ps%power_ec(i)=powers(i)
 !MKSA issue
                  ps%freq_ec(i)=frqncy  !PS in Hz
               enddo
                  
               
            elseif (rfmode.eq.'LH') then
               ps%lhrf_src_name=achar(48+i) !achar(48)='0'
               if (ps%nlhrf_src.eq.0) then
                  ps%nlhrf_src=ngenray_source
                  call ps_alloc_plasma_state(ierr) !set PS dims for rf
               else
                  write(*,*) 
     +              'prepare_genray_input: resetting PS LH source count'
                  write(*,*) 'from ',ps%nlhrf_src,' to ',ngenray_source
                  
!     copy all EXCEPT EC component profiles
                  cclist = ps_cc_all
                  cclist(ps_cnum_LH)=0
                  
                  call ps_copy_plasma_state(ps,aux,ierr,cclist= cclist)
                  if(ierr.eq.0) then
!     OK, copy back to ps
                     call ps_copy_plasma_state(aux, ps, ierr)
                     
!                     if(ierr.eq.0) then
!     set desired dimension
                        ps%nlhrf_src=ngenray_source
                        call ps_alloc_plasma_state(ierr) !set dims in PS
!                     endif !on ierr
                  endif !on ierr

               endif !on ps%nlhrf_src

               do i=1,ngenray_source
                  ps%lhrf_src_name=achar(48+i)  !achar(48)='0'
                  if (istart.eq.1) ps%power_lh(i)=powtot(i)
                  if (istart.gt.1) ps%power_lh(i)=powers(i)
!MKSA issue
                  ps%freq_lh(i)=frqncy
               enddo
                  
            elseif (rfmode.eq.'IC') then
               if (ps%nicrf_src.eq.0) then
                  ps%nicrf_src=ngenray_source
                  call ps_alloc_plasma_state(ierr) !set PS dims
               else
                  write(*,*) 'prepare_genray_input: '//
     +                 'resetting PS IC source count'
                  write(*,*) 'from ',ps%nicrf_src,
     +                 ' to ',ngenray_source
                  
!     copy all EXCEPT EC component profiles
                  cclist = ps_cc_all
                  cclist(ps_cnum_IC)=0
                  
                  call ps_copy_plasma_state
     +                 (ps,aux,ierr,cclist=cclist)
                  if(ierr.eq.0) then
!     OK, copy back to ps
                     call ps_copy_plasma_state(aux, ps, ierr)
                     
!                     if(ierr.eq.0) then
!     set desired dimension
                        ps%nicrf_src=ngenray_source
                        call ps_alloc_plasma_state(ierr) !set dims in PS
!                     endif      ! on ierr
                  endif         ! on ierr
                  
               endif            ! on ps%nicrf_src

               do i=1,ngenray_source
                  ps%icrf_src_name=achar(48+i)  !achar(48)='0'
                  write(*,*)'icrf_src_name=',achar(48+i)
                  if (istart.eq.1) ps%power_ic(i)=powtot(i)
                  if (istart.gt.1) ps%power_ic(i)=powers(i)
!MKSA issue
                  ps%freq_ic(i)=frqncy
               enddo
               
            endif     ! on rfmode
            
         endif        ! on i_check_rf.eq.0 .and. PS_add_nml.eq.'force'
         
      endif           ! on (PS_add_nml.eq.'enabled' .or. PS_add_nml.eq.'force
                     

!.......................................................................
!     In general, as many genray namelist parameters as possible will be
!     specified only through the genray.in/.dat file, usually set up
!     for a particular rf mode and machine.  Depending on these
!     settings, additional data will be obtained from the plasma state.
!.......................................................................



!.......................................................................
c     Interpolate PS profiles to user (i.e., genray) r_nc-grid:
c     Use rezoning subroutine from PS, onto 
c     rho_bdy_rezone(...) which gives PS values at the r_nc()
c     bin boundaries (except values centered in last half bin at
c     plasma center and edge are taken as r_nc=0. and 1. values).
c
c     That is, size(rho_bdy_rezon) = size(r_nc) + 1
c     This gives r_nc bins with r_nc+1 bin boundaries.
c     Plasma profile values at the centers of the bins are taken 
c     to be values at the user r_nc-grid bin boundaries.
!.......................................................................

c     Equispaced, bin-boundary grid for profiles to genray
      dr_nc=1.0_rspec/(nj-1)
      do l=1,nj
         r_nc(l)=(l-1)*dr_nc
      enddo

      allocate(rho_bdy_rezon(njp1))
      rho_bdy_rezon(1)=r_nc(1)
      rho_bdy_rezon(njp1)=r_nc(nj)
      rho_bdy_rezon(2:nj) =
     +     .5_rspec*(r_nc(1:nj-1)+r_nc(2:nj))

!.......................................................................
c     First treat electrons and zeff:
!.......................................................................
      
      allocate(dense(nj),tempe(nj),zeff(nj))

      call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(0), dense, ierr,
     +     zonesmoo=.TRUE.)
      write(*,*)'after ps_rho_rezone: ierr=',ierr
      call ckerr('ps_rho_rezone (U1)',ierr)

      call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(0), tempe, ierr,
     +     zonesmoo=.TRUE.)
      call ckerr('ps_rho_rezone (U2)',ierr)

      call ps_rho_rezone(rho_bdy_rezon, ps%id_Zeff, zeff, ierr,
     +     zonesmoo=.TRUE.)
      call ckerr('ps_rho_rezone (U2)',ierr)


c     Checking interpolation
      write(*,*)
      write(*,*)'Checking interpolation'
      write(*,*)'   l  r_nc    rho_bdy_rezon   dense    tempe    zeff'
      do l=1,nj
      write(*,1001)l,r_nc(l),rho_bdy_rezon(l),dense(l),tempe(l),zeff(l)
      enddo
 1001 format(i4,5(1pe12.4))


!.......................................................................
!     Fill in _nc electron variables for genray.in file
!.......................................................................

      one=1._rspec
      charge_nc(1)=one
      dmas_nc(1)=one
      do j=1,nj
         en_nc(j,1)=dense(j)
         temp_nc(j,1)=tempe(j)
         zeff_nc(j)=zeff(j)
      enddo
      incl_species(0)=1   !Using PS index convention here.
      
      
      
!.......................................................................
!     If ion species are to be passed to genray, 
!     fill in arrays for electrons and ions, and write out
!     to genray.in file.
!.......................................................................

      if (nbulk.gt.1) then
         
         allocate(densi(nj,ps%nspec_alla),tempi(nj,ps%nspec_alla))
         densi(:,:)=0.0_rspec
         tempi(:,:)=0.0_rspec
         
!     get the abridged thermal species data first...
         
         call ps_tha_rezone(rho_bdy_rezon, ierr, zonesmoo=.TRUE.,
     +       ns = densi(:,1:ps%nspec_tha), Ts = tempi(:,1:ps%nspec_tha))
         call ckerr('ps_tha_rezone',ierr)

         incl_species(1:ps%nspec_tha)=1
         
!     get the fast species; in general these will be on different grids.
         
         allocate(wk_dens(nj),wk_eperp(nj),wk_epll(nj))
         
!     Need to initialize ierr,ierr1,ierr2,ierr2. for the case the
!     ps_rho_rezone  is not to be executed, per PB experience:
         ier=0
         ier1=0
         ier2=0
         ier3=0

         do ii=ps%nspec_tha + 1, ps%nspec_alla
            jj = ps%alla_index(ii) ! index in "all species" list

            write(*,*)'ii,jj,ps%nspec_tha,ps%nspec_alla',
     +                 ii,jj,ps%nspec_tha,ps%nspec_alla

            write(*,*)'ps%alla_type(0:ps%nspec_alla),ps_beam_ion',
     +                 ps%alla_type(0:ps%nspec_alla),ps_beam_ion

!     jj=-1 if ps%alla_type(ii).eq.ps_impurity, but this 
!     should not occur as here we loop over fast ions only.
!     Still, I test this...
            if(jj.le.0) then
               write(*,*) '? unexpected alla_index(ii) value: ',jj
            else
!     OK, jj=alla_index(ii) points to the element index of whatever
!     fast species list element corresponds to element ii of the alla list
               wk_dens(:)=0.0_rspec
               wk_eperp(:)=0.0_rspec
               wk_epll(:)=0.0_rspec

               if(ps%alla_type(ii).eq.ps_beam_ion) then
 
                if (allocated(nbeami)) then
                incl_species(ii)=1
                call ps_rho_rezone(rho_bdy_rezon, ps%id_nbeami(jj),
     +                 wk_dens, ier1, zonesmoo=.TRUE.)
                call ckerr('ps_rho_rezone: ier1',ier1)
                call ps_rho_rezone(rho_bdy_rezon, ps%id_eperp_beami(jj),
     +                 wk_eperp, ier2, zonesmoo=.TRUE.)
                call ps_rho_rezone(rho_bdy_rezon, ps%id_epll_beami(jj),
     +                 wk_epll, ier3, zonesmoo=.TRUE.)
                endif
              
                  
               else if(ps%alla_type(ii).eq.ps_rf_minority) then
                if (allocated(nmini)) then
                incl_species(ii)=1
                call ps_rho_rezone(rho_bdy_rezon, ps%id_nmini(jj),
     +                 wk_dens, ier1, zonesmoo=.TRUE.)
                call ps_rho_rezone(rho_bdy_rezon, ps%id_eperp_mini(jj),
     +                 wk_eperp, ier2, zonesmoo=.TRUE.)
                call ps_rho_rezone(rho_bdy_rezon, ps%id_epll_mini(jj),
     +                 wk_epll, ier3, zonesmoo=.TRUE.)
                endif
                  
               else if(ps%alla_type(ii).eq.ps_fusion_ion) then
                if (allocated(nfusi)) then
                incl_species(ii)=1
                call ps_rho_rezone(rho_bdy_rezon, ps%id_nfusi(jj),
     +                 wk_dens, ier1, zonesmoo=.TRUE.)
                call ps_rho_rezone(rho_bdy_rezon, ps%id_eperp_fusi(jj),
     +                 wk_eperp, ier2, zonesmoo=.TRUE.)
                call ps_rho_rezone(rho_bdy_rezon, ps%id_epll_fusi(jj),
     +                 wk_epll, ier3, zonesmoo=.TRUE.)
                endif
                  
               endif  !  On ps%alla_type(ii)
               
               ierr = ier1 + ier2 + ier3
               call ckerr('ps_rho_rezone (fast specie)',ierr)
               
               densi(:,ii) = wk_dens

!     McCune rational for following is that wk_perp and wk_epll
!     are average perp and parallel energies.  (For isotropic
!     distns, wk_eperp will be 2.*wkpell.)

               tempi(:,ii) = (2.0_rspec/3.0_rspec)*(wk_eperp + wk_epll)
               
            endif  !  On jj
         enddo  !  On ii

         write(*,*)'incl_species(0:ps%nspec_all)=',
     +              incl_species(0:ps%nspec_all)

!     Determine reduced species list if any densities are ~zero inside
!     LCFS, or temperatures zero, including LCFS.  
!     Indicator is stored in incl_species(0:ps%nspec_alla).

         eps_dens=1.e-12_rspec
         eps_temp=1.e-12_rspec
         write(*,*)'eps_dens,eps_temp=',eps_dens,eps_temp

!     Check electron density/temperature
         do j=1,nj-1
            if(abs(dense(j)).lt.eps_den) incl_species(0)=0
            if(abs(tempe(j)).lt.eps_temp) incl_species(0)=0
         enddo
         if(abs(tempe(nj)).lt.eps_temp) incl_species(0)=0

         if (incl_species(0).eq.0) then
            write(*,*)
            write(*,*)'Electron density or temperature ~zero'
            write(*,*)"Can't ray trace. Will stop"
            write(*,*)
            write(*,*) 'STOP'
            STOP
         endif

!     Check ion densities/temperatures
         do ii=1,ps%nspec_alla
            do j=1,nj-1
               if(abs(densi(j,ii)).lt.eps_dens) incl_species(ii)=0 
               if(abs(tempi(j,ii)).lt.eps_temp) incl_species(ii)=0
            enddo
            if(abs(tempi(nj,ii)).lt.eps_temp) incl_species(ii)=0
         enddo

         write(*,*)'incl_species(0:ps%nspec_alla)=',
     +              incl_species(0:ps%nspec_alla)

         jj=1  ! genray index for electrons
         do ii=1,ps%nspec_alla
            if (incl_species(ii).ne.0) then
               jj=jj+1
               charge_nc(jj)=abs(ps%q_alla(ii)/ps_xe)
               dmas_nc(jj)=ps%m_alla(ii)/ps_me
               do j=1,nj
                  en_nc(j,jj)=densi(j,ii)
                  temp_nc(j,jj)=tempi(j,ii)
               enddo ! On j
            endif
         enddo ! On ii
         nspecgr=jj
         nbulk=nspecgr
         ndens=nj

         write(*,*)'ps%nspec_alla=',ps%nspec_alla
         write(*,*)'ps%nspec_tha=',ps%nspec_tha
         write(*,*)'nspecgr, nbulk=',nspecgr,nbulk
         write(*,*)'charge_nc=',charge_nc(1:nspecgr)
         write(*,*)'dmas_nc=',dmas_nc(1:nspecgr)


         ll=1
         iout=5
c         write(iout,'(A,1pe12.5,i3)') '  density      Temperature'//
         write(*,'(A,1pe12.5,i3)') '  density      Temperature'//
     +        '   atm_no   charge   Species @ r_nc(ll),ll= ',r_nc(ll),ll

!BH111108         do ii=1,ps%nspec_alla
         do ii=1,nbulk
!             Write(iout,'(4(1x,1pe12.5),3x,a)') densi(1,ii+1),
             Write(*,'(4(1x,1pe12.5),3x,a)') en_nc(ll,ii),
     +       temp_nc(ll,ii), dmas_nc(ii), charge_nc(ii)
         enddo  ! On ii

      endif ! On nbulk.gt.1

      do k=1,nbulk
         charge(k)=charge_nc(k)
         dmas(k)=dmas_nc(k)
      enddo


      do k=1,nbulk
         do j=1,ndens
            dens1(j,k)=en_nc(j,k)
            temp1(j,k)=temp_nc(j,k)
            tpop1(j,k)=1.d0
            vflow1(j,k)=0.d0
         enddo
      enddo
      do j=1,ndens
         zeff1(j)=zeff_nc(j)
      enddo


c---------------------------------------------------------------
c     write all namelist to  genray.in file, and check/adjust for
c     the comma problem (see subroutine comma_check)
c----------------------------------------------------------------

      write(*,*)'before write_all_namelists'

      call write_all_namelists(ndim)
      call comma_check('genray.in')

      write(*,*)'after write_all_namelists and comma_check'

c---------------------------------------------------------------
c     if have modified PS, then store it.
c---------------------------------------------------------------
      if (PS_add_nml.ne.'disabled') then
c         write(iout,*)
         write(*,*)
     +        'prepare_genray_input: --storing genraynml data '//
     +        'in current PS'
!     Stores the modified plasma state--doesn't commit it. 
!     The previous one is still around
         
         call ps_store_plasma_state(ierr)
         if(ierr .ne. 0) then
c            write(iout,*)
            write(*,*)
     +           'Cannot ps_store_plasma_state in prepare_genray_input'
         endif
         
      endif                     ! on PS_add_nml

c---------------------------------------------------------------
c     end
c---------------------------------------------------------------
      write(*,*)
      write(*,*)'prepare_genray_input: Normal end'
      write(*,*)

 999  continue

      end  ! prepare_genray_input.f90

c----------------------------------------------------------------
c----------------------------------------------------------------



      SUBROUTINE ckerr(sbrtn,ierr)
      character*(*), intent(in) :: sbrtn
      write(*,*) 'ckerr: ierr=',ierr
      IF(ierr.NE.0) then
         write(6,*) ' ?plasma_state_test: error in call: '//trim(sbrtn)
         stop
      ENDIF
      END SUBROUTINE ckerr



      subroutine transcribe(genraynml,genraynml1)
c
c     Transcribes file genraynml to genraynml1, except do nothing
c     if genraynml1 is same file as genraynml.
c     This 'copy' is in lieu of access to operating system commands.
c
      character(len=*), intent(in) :: genraynml,genraynml1
      character(len=256) :: genraynml2
      parameter(long_enough=100000)
      character(len=long_enough) :: line1          !Automatic array
      logical logic1

      write(*,*)
      write(*,*)  'transcribe genraynml to '//trim(genraynml1)
      write(*,*)

      max_length=0

      write(*,*)'transcribe: genraynml =',genraynml
      inquire(file=genraynml,iostat=kiostat,opened=logic1,
     +     number=inumber)
      write(*,*)'transcribe: inquire on ',genraynml,', opened=',logic1,
     1   'iostat=',kiostat,'unit=',inumber
      
c     Adjust for case genraynml1 does/does not include "./"
      genraynml2=trim(genraynml1)
      if (genraynml2(1:2) .eq. './') then
         genraynml2(1:2)='  '
         genraynml2=trim(genraynml2)
      endif
      if (trim(genraynml).eq.genraynml1
     1    .or. (trim(genraynml).eq.genraynml2) ) then
         write(*,*)'transcribe: genraynml1 is same file as genraynml'
         return
      endif

c     Open genraynml and transcribe it
      open(unit=20,file=genraynml,delim='none',status="old",
     +     iostat=kiostat)
      if (kiostat.ne.0) write(*,*)'transcribe: kiostat=',kiostat
      if (kiostat.ne.0) STOP 'transcribe: prblm with genraynml'

      open(unit=21,file=genraynml1,delim='none',status='replace',
     +     iostat=kiostat)
      if (kiostat.ne.0) write(*,*)'transcribe: kiostat=',kiostat
      if (kiostat.ne.0) STOP 'transcribe: prblm with opening genraynml1'

 3    read(unit=20,fmt='(a)',end=4) line1
      len_line1=len_trim(line1)
      if (len_line1.gt.max_length) max_length=len_line1
cBH131022      write(unit=21,fmt=*) trim(line1)
      write(unit=21,fmt='(a)') trim(line1)
      if (len_line1.ge.(long_enough-1)) then
         write(*,*)'transcribe:len_line1,long_enough',
     +             len_line1,long_enough
         STOP 'Adjust long_enough'
      endif
      go to 3
 4    continue

      write(*,*)'transcribe:  max line length =',max_length
      close(20)
      close(21)
      write(*,*)'transcribe:  complete'
      write(*,*)
 
      return
      end subroutine transcribe



      subroutine comma_check(genraynml)
c
c     This is a fix of a "comma-problem" occuring for namelist
c     writes in the intel compiler, swim setup.  Under special
c     circumstances of namelist real*8 array variable name length, 
c     the line wrap causes an out of place comma to be generated,
c     which is read as an extra delimiter when subsequently
c     reading the namelist (at least with pgi compiled genray).
c
c     genraynml (usually a genray.in-like file) is read, checked
c     for a comma before other non-blank characters, in columns 1:5.
c     If a comma is found, it is replaced by a blank. genraynml is
c     then written to file genray.in.

      character(len=*), intent(in) :: genraynml
      character(len=256) :: genraynml_tmp
      parameter(long_enough=100000)
      character(len=long_enough) :: line1          !Automatic array
      logical logic1


      write(*,*)'comma_check: genraynml =',genraynml
      inquire(file=genraynml,iostat=kiostat,opened=logic1,
     +     number=inumber)
      write(*,*)'comma_check: inquire on ',genraynml,', opened=',logic1,
     1   'iostat=',kiostat,'unit=',inumber
      if (kiostat.ne.0) STOP 'comma_check: prblm with genraynml'

c     Open genraynml and transcribe it the temporary file (can\'t
c     change genraynml in place)
      genraynml_tmp=trim(genraynml)//'.tmp'
      call transcribe(genraynml,genraynml_tmp)

      max_length=0
      n_extra_commas=0

c     Open genraynml_tmp
      open(unit=20,file=genraynml_tmp,delim='none',status="old",
     +     iostat=kiostat)
      if (kiostat.ne.0) write(*,*)'comma_check: kiostat=',kiostat
      if (kiostat.ne.0) STOP 'comma_check: prblm with genraynml_tmp'

      open(unit=21,file=genraynml,delim='none',status='replace',
     +     iostat=kiostat)
      if (kiostat.ne.0) write(*,*)'transcribe: kiostat=',kiostat
      if (kiostat.ne.0) STOP 'transcribe: prblm with opening genraynml'

 3    read(unit=20,fmt='(a)',end=4) line1
      len_line1=len_trim(line1)
      if (len_line1.gt.max_length) max_length=len_line1

c     Check for extra comma
      do i=1,5
	if (line1(i:i).ne." ") then
	  if (line1(i:i).eq.",") then
	     line1(i:i)=" "
	     n_extra_commas=n_extra_commas+1
	  endif
	  exit  !Only allow for 1 comma, exit do-loop if other than comma.
	endif
      enddo	

      write(unit=21,fmt='(a)') trim(line1)
      if (len_line1.ge.(long_enough-1)) then
         write(*,*)'transcribe:len_line1,long_enough',
     +             len_line1,long_enough
         STOP 'Adjust long_enough'
      endif
      go to 3
 4    continue

      if (n_extra_commas.ne.0) then
	write(*,*)'comma_check: n_extra_commas=',n_extra_commas
      endif

      close(20)
      close(21)
 
      return
      end subroutine comma_check
