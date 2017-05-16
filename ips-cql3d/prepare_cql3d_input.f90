!
!
      program prepare_cql3d_input

!Invoked as:  prepare_cql3d_input [ips_mode cql3d_mode cql3d_output
!             cql3d_nml restart nsteps_str deltat_str ps_add_nml, enorm (optional) ]
!             Up to 9 command line arguments, as described below.
!             Output is cqlinput namelist file, adjusted with PS data.
!
!             Probably better to code some of these inputs as optional
!             arguments.
!
! BH120705: Adjustments for prepare_cql3d_input.

!.......................................................................
!     For ips_mode='step':
!     Prepares new cqlinput  namelist file from an input
!     namelist "template" file (cql3d_nml) specified on the command
!     line. The resulting cqlinput file contains the "template" data
!     modified by Plasma State 2 (PS) plasma profiles and PS abridged
!     species list (if ions involved), plus PS source specifications
!     for particular cql3d_mode and cql3d_output command line inputs.
!     A list of the PS source specification variables passed to 
!     cqlinput will be maintained below.
!     Prepare_cql3d_input uses several include files from the cql3d
!     source, plus cqlinput variables specified through cql3d_nml,
!     and PS2 variables.
!
!     See Plasma_State_V2.003.pdf (May 2008) re Abridged Species.
!
!     prepare_cql3d_input takes up to 8 optional command line
!     arguments:   1st: ips_mode, either 'init' or 'step' for coupling
!                       to SWIM Integrated Plasma Simulator (IPS).
!                       ips_mode='init', only initializes some
!                       component dimensions in the plasma state.
!                       default='step'
!                  2nd: cql3d_mode: 'el-only' 'el+ions' 'ions+el'
!                       are possible inputs.  
!                       In the 'el-only' case, only electron  and Zeff
!                       PS plasma profiles and Zeff are used. Two additional
!                       ions are setup for cql3d in accord with the input
!                       template cqlinput file.
!                       In 'el+ions' case, electron and ion densities from
!                       the PS abridged list are set up.  Electrons are
!                       a cql3d general species. Zeff not used.
!                       In 'ions-abr' case, PS electron +abridged ion list
!                       profiles are input to cqlinput.  The first ion
!                       will be a cql3d general species.  Future work
!                       can modify this to two or more ion general species.
!                       default='el-only'
!                  3rd: cql3d_output, gives types of output power 
!                       depostion and currents
!                       (One of 'LH', 'EC' 'IC' 'NBI' 'RW' 
!                       'LH+RW' 'NBI+IC').
!                       (Plasma state component RW refers to RunaWay elec.)
!                       Present implementation is 'LH', 'EC', 'IC', 'NBI', 
!                       'RW', 'LH+RW'.  (Future: dimension cql3d_output
!                       and add number of elements to treat simultaneously.)
!                  4th: cql3d_nml, full path specification of template
!                       cql3d nml file,   default='./cqlinput' (in
!                       which case cqlinput will be modified in place)
!                  5th: restart, 'disabled' indicates initial step, starting
!                       from internal Maxwl; 'enabled', restarts from previously
!                       generated distrfunc.nc netcdf file, which is generally
!                       the netcdf output files from a previous run containing
!                       grids and the general distribution function.
!                  6th: nsteps_str, number of time steps by cql3d (default='10')
!                  7th: deltat_str, cql3d time step in secs (default='0.001')
!                  8th: ps_add_nml, indicates whether some data,
!                       indicated below, should be added to the plasma
!                       state from cqlinput namelist.
!                       Valid values: 'disabled','enabled','force'.
!                       E.G., Can be used to update a test PS.
!                       default='disabled' [the usual situation].
!                  9th: enorm (optional) (added by DBB 5-16-2017) Set in call from 
!                       fp_cql3d_general.py which gets it from IPS config
!                       parameter.
!
!     CONSIDER: adding enorm and jx arguments, to adjust if changed
!     from cqlinput.
!
!     To set a given command line argument to non-default value, it
!     is necessary to set the preceeding values in the list.
!     (Consider making some arguments optional.)
!
!     In general, for RF cases, genray (or with minor adjustments
!     an alternative RF code) will be called BEFORE cql3d.
!     This will set PS RF variables.
!     However, there will be a check on whether plasma state
!     RF variables are defined.
!     If not, and command line variable PS_add_nml.eq.'enabled',
!     then will add data to the PS from the genraynml data.
!     If RF variables are defined or not, but PS_add_nml.eq.'force',
!     the values from genraynml data will be inserted in the PS.
!     Normally, PS_add_nml will be 'disabled', and genray.in data
!     will be set from the PS data.
!
!     cql3d requires several auxiliary input files.  Full paths can
!     be given in cqlinput.  Usually, it will be easiest to simply
!     stage them in the work directory, and give them the default
!     names used in cql3d (e.g., genray.nc, du0u0_input).
!     The eqdsk location is set according to the ps%eqdsk_file.
!
!     This code can be easily adjusted for coupling to aorsa
!     or toric.  Perhaps cql3d_mode='ions+el' will be sufficient.
!     Alternatively, add cql3d_mode='ions-aor' and/or 'ions-tor'.
!
!.......................................................................



!.......................................................................
!     List of RF-related PS variables passed to cqlinput, 
!     overwriting data from cql3d_nml (in addition to plasma profiles
!     and the eqdsk file path), depending on cql3d_mode/cql3d_output:
!     (If PS_add_nml.ne.'disabled' these values may alternatively be
!     set in the PS using data from cql3d_nml.)
!
!
!     For PS_add_nml.eq.'enabled', only add appropriate freq_xx
!     to PS, if ps%freq_xx=0, i.e., undefined.
!
!
!          The PS number of sources should be .ge. the number
!          cumulatively specified in genraynml, for a given rfmode.
!          This will be checked.
!
!
!.......................................................................
!     List of cql3d namelist variables modified by prepare_cql3d_input
!     in addition to plasma profile data:
!     [This list should be maintained, if there are additions.]
!     
!     Various scaling of profiles set for unit conversion, as necessary.
!     pwrscale(1)=1.0  RF power
!     bptor(1)=ps%power_nbi(1)   NBI power passed from ps to cqlinput
!     
!     
!.......................................................................
!
!
!.......................................................................
!
!     Some files are used from the cql3d distribution:
!     =================================================
!      -the subset of include files in ./cql3d_includes
!      -aindfpa.f,aindflt.f,eqindflt.f,urfindfl.f,frinitl.f
!
!     Using these cql3d distn files is aimed at keeping cql3d and 
!     prepare_cql3d_input.f90 in sync.   When any of these files change 
!     in the cql3d distribution, the files should be updated in the
!     present code suite.
!
!.......................................................................
!
!     Compile and load line:
!     ======================
!     Include files in subdir ./cql3d_includes:  param.h, etc.
!                             Above .f files.
!viz.pppl.gov compile:
!     ifort --fix -o prepare_input prepare_cql3d_input.f aindfpa.f aindflt.f eqindflt.f urfindfl.f frinitl.f
! ifort -fixed -80 -c -O -w -mp -fpe0 -heap-arrays 20 -I./cql3d_includes -I$NTCCHOME/mod -I$NETCDFHOME/include -o prepare_cql3d_input.o  prepare_cql3d_input.f90
!
!Check which modules are loaded with 'module list'.
!
!.......................................................................


!.......................................................................
! The following USE gives access to the extensive plasma state data.
! (For genray, it is necessary to rename "output" in the plasma state to
! avoid conflict with genray namelist /output/. For cql3d, it is benign.)

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

      include 'param.h'
c----------------------------------------
c     type declaration and storage of all namelists variables
c-----------------------------------
      include 'name.h'
      include 'name_decl.h'
      include 'frname.h'
      include 'frname_decl.h'


      character*8 ips_mode,cql3d_mode,cql3d_output,restart,ps_add_nml
      character(swim_filename_length) cql3d_nml
      character(len = 128) :: deltat_str, nsteps_str, enorm_str = ' '
      character(len = 128) :: nlrestrt_str
      REAL(KIND=rspec) :: deltat
      LOGICAL :: enorm_str_present

      integer nj,njp1,nrho
      integer nsteps
      REAL(KIND=rspec) :: dryain

!!$      REAL(KIND=rspec) charge_nc(nspec_maxa),dmas_nc(nspec_maxa),
!!$     +     r_nc(nraa),en_nc(nraa,nspec_maxa),temp_nc(nraa,nspec_maxa),
!!$     +     zeff_nc(nraa)
      integer, dimension(:), allocatable :: incl_species !indicator for
                                          !species to send to genray.

      integer iargs

c---------------------------------------------------------------
!  Unified plasma parameter arrays (in ps%nspec_alla list)
!  abridged thermal species + all fast species combined.

      REAL(KIND=rspec), dimension(:), allocatable :: rho_bdy_rezon
      REAL(KIND=rspec), dimension(:), allocatable :: dense,tempe,zeff
      REAL(KIND=rspec), dimension(:,:), allocatable :: densi,tempi
      REAL(KIND=rspec), dimension(:), allocatable :: wk_eperp,wk_epll,
     +                                               wk_dens

      REAL(KIND=rspec), dimension(:), allocatable :: rho_user
c---------------------------------------------------------------
      integer :: cclist(ps_ccount)  ! State partial copy controls.
                                    ! Ps_ccount is number of known
                                    ! components in the PS.
c-----------------------------------------------------------------------
!  Setup some allocatable variables from PS, which will be
!  checked if allocated:
      REAL(KIND=rspec),dimension(:),allocatable :: nbeami,nmini,nfusi

      REAL(KIND=rspec) :: pi,twopi,me,mp,xe,zero,one

c---------------------------------------------------------------

c     From plasma_state_spec.dat:
      pi=3.1415926535897931_rspec    ! pi
      twopi=6.2831853071795862_rspec !2*pi
      me  = 9.1094e-31_rspec         ! electron mass, kg
      mp  = 1.6726e-27_rspec         ! proton mass, kg
      xe  = 1.6022e-19_rspec         ! Joule/ptcl-eV (also, electron charge)
      zero     = 0.0_rspec           ! Zero
      one      = 1.0_rspec           ! one
      

c--------------------------------------------------------------
c     Read command line arguments
c--------------------------------------------------------------

c     Defaults:
      ips_mode='step'
      cql3d_mode='el-only'
      cql3d_output='LH'
      cql3d_nml='./cqlinput'
      restart='disabled'
      nsteps_str='10'
      deltat_str='0.1e-3'  !secs
      ps_add_nml='disabled'
      enorm_str_present = .FALSE.

c     F2003-syntax: command_argument_count()/get_command_argument(,)
      iargs=command_argument_count()
      write(*,*)'iargs=',iargs
      if (iargs.ne.8) then
         write(*,*)'prepare_cql3d_input usage: '
         write(*,*)'Up to eight command line arguments, '
         write(*,*)
     +   'ips_mode cql3d_mode cql3d_output cql3d_nml restart nsteps_str '
     +   //'deltat_str ps_add_nml (refer to code)'
      endif

      if (iargs.ge.1)   call get_command_argument(1,ips_mode)
      if (iargs.ge.2)   call get_command_argument(2,cql3d_mode)
      if (iargs.ge.3)   call get_command_argument(3,cql3d_output)
      if (iargs.ge.4)   call get_command_argument(4,cql3d_nml)
      if (iargs.ge.5)   call get_command_argument(5,restart)
      if (iargs.ge.6)   call get_command_argument(6,nsteps_str)
      if (iargs.ge.7)   call get_command_argument(7,deltat_str)
      if (iargs.ge.8)   call get_command_argument(8,ps_add_nml)
      if (iargs.ge.9) then
        call get_command_argument(9,enorm_str)
        enorm_str_present = .TRUE.
        enorm_str=trim(enorm_str)
      end if

      ips_mode=trim(ips_mode)
      cql3d_mode=trim(cql3d_mode)
      cql3d_output=trim(cql3d_output)
      cql3d_nml=trim(cql3d_nml)
      restart=trim(restart)
      nsteps_str=trim(nsteps_str)
      deltat_str=trim(deltat_str)
      ps_add_nml=trim(ps_add_nml)

      write(*,*)'prepare_cql3d_input command line arguments: ',
     +  ips_mode,'  ',cql3d_mode,'  ',cql3d_output,'  ',cql3d_nml,
     +  '  ',restart,'  ',nsteps_str,'  ',deltat_str,'  ',ps_add_nml,
     + ' ', enorm_str


c-----------------------------------------------------------------------
c     Copy cql3d_nml to ./cqlinput
c     Transcribe file, since not all fortrans have access
c     to shell commands (for which simple cp would work).
c-----------------------------------------------------------------------

      call transcribe(cql3d_nml)


c-----------------------------------------------------------------------
c     Create default input data
c-----------------------------------------------------------------------

      call aindfpa
      call aindflt
      call eqindflt
      call urfindfl
      call frinitl

c-----------------------------------------------------------------------
c     Read all namelists [storage setup in name_decl.h,frname_decl.h].
c-----------------------------------------------------------------------

      open(unit=2, file = 'cqlinput', delim='apostrophe',
     1     status = 'old', form = 'formatted', iostat = ierr)
      if(ierr .ne. 0) stop 'cannot open cqlinput cql3d namelist file'

      read(2,setup0,iostat=istat)
      if (istat.ne.0) then
         write(*,*)
         write(*,*)
     1   'No setup0?: Check have new setup0/setup namelist structure.'
         write(*,*)
         STOP 1
      endif
      read(2,setup)     !Note:  Now only one namelist setup section
      read(2,trsetup)
      read(2,sousetup)
      read(2,eqsetup)
      read(2,rfsetup)
      read(2,frsetup)
      close(2)

      write(*,*)'after read_all_namelists'
      write(*,*) 'fmass (2) = ', fmass(2) ! ptb:
      write(*,*) 'nstop = ', nstop ! ptb:
!.......................................................................
!     Make available all variables in the plasma state module, specified
!     in swim_state_spec.dat (or corresponding .txt) in components/
!     state/src/plasma_state.  See also McCune's Plasma_State_Vnn
!     design doc at http://www.cswim.org/componentdescrip/.
!     Also, swim_state_spec.dat (or corresponding *.txt) in
!     components/state/src/plamsa_state.
!     These variables are available with the prepended "component selector"
!     for the current state, ps%, or for the previous state,  %psp,
!     and a few others.
!.......................................................................

      call ps_get_plasma_state(ierr,filename='./cur_state.cdf')
      if(ierr .ne. 0) stop 'cannot get plasma state needed for profiles'


!.......................................................................
!     If ips_mode.eq.'init', then set approrpriate dimensions in the 
!     plasma state component section, dependent on rfmode
!.......................................................................

      if (ips_mode.eq.'init') then
         if (cql3d_output.eq.'EC') then
            
            nrho=lrz+1  !lrz is number of fp bin centers spec'd in cqlinput

            write(*,*)'nrho=lrz+1 frm nml, ps%nrho_ecrf frm PS =',
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
     +             ' * prepare_cql3d_input: reset PS EC profile size'
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
                     call ps_alloc_plasma_state(ierr) !set these PS dims in PS
                  endif
               endif  !on ierr
            endif
            
            if (ierr.ne.0) then
               write(*,*)'prepare_cql3d_input:  EC'//
     +                     ' ps_alloc_plasma_state: ierr=',ierr
               stop
            endif

c           Set ec calc code name
            ps%ec_code_info='cql3d'

c           Initialize arrays
c           BH120731:  Consider use of ps_rho_rezone (as below).
            call zone_check(nrho,ps%rho_ecrf,rya(1))  !rya dims 0:...

            ps%peech=0.0_rspec    !dimension 1:nrho
            ps%curech=0.0_rspec   !dimension 1:nrho

         elseif(cql3d_output.eq.'LH' .or. cql3d_output.eq.'LH+RW' ) then

            nrho=lrz+1  !lrz is number of bin centers spec'd in cqlinput

            write(*,*)'nrho frm nml, ps%nrho_lhrf frm PS =',
     +	               nrho,ps%nrho_lhrf

!           Necessary PS dimensions (user has to set his component dims):
            if(ps%nrho_lhrf.eq.nrho) then
               continue         ! dimension OK already
            else if(ps%nrho_lhrf.eq.0) then
               ps%nrho_lhrf=nrho
               call ps_alloc_plasma_state(ierr) !set these PS dims in PS
               write(*,*)'ps%nrho_lhrf set= ',nrho
c           Following if clause should not be satisfied in ordinary IPS usage
c           (Can occur if using already initialized PS)
! ptb don't execute this part of the else if endif 
!            else if(ps%nrho_lhrf.ne.nrho) then
!               write(*,*) 
!     +             ' * prepare_cql3d_input: reset PS LH profile size'
!               write(*,*) ' * from ',ps%nrho_lhrf,' to ',nrho
!               
!              copy all EXCEPT LH component profiles
!               cclist = ps_cc_all
!               cclist(ps_cnum_LH)=0
!
!               call ps_copy_plasma_state(ps, aux, ierr, cclist = cclist)
!               if(ierr.eq.0) then
!                  ! OK, copy back to ps
!                  call ps_copy_plasma_state(aux, ps, ierr)
!                  
!                  if(ierr.eq.0) then
!                     ! set desired dimension
!                     ps%nrho_lhrf=nrho
!                     call ps_alloc_plasma_state(ierr) !set these PS dims in PS
!                  endif
!               endif
            endif
            
            if (ierr.ne.0) then
               write(*,*)'prepare_cql3d_input:  LH'//
     +                     ' ps_alloc_plasma_state: ierr=',ierr
               stop
            endif

c           Set lh calc code name
            ps%lh_code_info='cql3d'

c           Initialize arrays, 
            write(*,*) 'Before rezoning: rho_lhrf =', ps%rho_lhrf(:)
            write(*,*) 'CQL3D grid is rya =', rya(1:nrho-1)
! ptb            call zone_check(nrho,ps%rho_lhrf,rya(1))  !rya dims 0:...
! ptb added this hack to get the preprare_cql3d_input to run through
! will need to do something better later on
! Won't need to recalculate ps%rho_lhrf if we do not allow it to be changed
!            do itr = 2,nrho-1
!               ps%rho_lhrf(itr) = 0.5* (rya(itr)+rya(itr-1))
!            enddo
!            ps%rho_lhrf(1) = 0.0_rspec
!            ps%rho_lhrf(nrho) = 1.0_rspec
! ptb - end of hack
!            write(*,*) 'After rezoning: rho_lhrf =', ps%rho_lhrf(1:nrho)

            ps%pelh=0.0_rspec    !dimension 1:nrho
            ps%pilh=0.0_rspec
            ps%curlh=0.0_rspec
	    ps%pelh_src=0.0_rspec
	    ps%pilh_src=0.0_rspec
	    curlh_src=0.0_rspec

         elseif(cql3d_output.eq.'IC' .or. cql3d_output.eq.'NBI+IC') then

            nrho=lrz+1  !lrz is number of bin centers spec'd in cqlinput


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
     +              ' * prepare_cql3d_input: reset PS IC profile size'
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
                     call ps_alloc_plasma_state(ierr) !set these PS dims in PS
                  endif
               endif
            endif
            
            if (ierr.ne.0) then
               write(*,*)'prepare_cql3d_input:  IC'//
     +                     ' ps_alloc_plasma_state: ierr=',ierr
               stop
            endif

c           Set ec calc code name
            ps%ic_code_info='cql3d'

c           Initialize arrays, 
            call zone_check(nrho,ps%rho_icrf,rya(1))  !rya dims 0:...

            ps%picrf_totals=0.0_rspec
            ps%picth=0.0_rspec
            ps%curich=0.0_rspec
            ps%cdicrf=0.0_rspec
            ps%curich=0.0_rspec

         elseif(cql3d_output.eq.'NBI'.or. cql3d_output.eq.'NBI+IC') then

            nrho=lrz+1  !lrz is number of bin centers spec'd in cqlinput


            write(*,*)'nrho frm nml, ps%nrho_icrf frm PS =',
     +	               nrho,ps%nrho_icrf

            nrho=lrz+1  !lrz is number of bin centers spec'd in cqlinput


            write(*,*)'nrho frm nml, ps%nrho_nbi frm PS =',
     +	               nrho,ps%nrho_nbi

!           Necessary PS dimensions (user has to set his component dims):
            if(ps%nrho_nbi.eq.nrho) then
               continue         ! dimension OK already
            else if(ps%nrho_nbi.eq.0) then
               ps%nrho_nbi=nrho
               call ps_alloc_plasma_state(ierr) !set these PS dims in PS
               write(*,*)'ps%nrho_nbi set= ',nrho
c           Following if clause should not happen in ordinary IPS usage
c           (Can occur if using already initialized PS)
            else if(ps%nrho_nbi.ne.nrho) then
               write(*,*) 
     +              ' * prepare_cql3d_input: reset PS NBI profile size'
               write(*,*) ' * from ',ps%nrho_nbi,' to ',nrho
               
!              copy all EXCEPT NBI component profiles
               cclist = ps_cc_all
               cclist(ps_cnum_NBI)=0

               call ps_copy_plasma_state(ps, aux, ierr, cclist = cclist)
               if(ierr.eq.0) then
                  ! OK, copy back to ps
                  call ps_copy_plasma_state(aux, ps, ierr)
                  
                  if(ierr.eq.0) then
                     ! set desired dimension
                     ps%nrho_nbi=nrho
                     call ps_alloc_plasma_state(ierr) !set these PS dims in PS
                  endif
               endif
            endif
            
            if (ierr.ne.0) then
               write(*,*)'prepare_cql3d_input:  NBI'//
     +                     ' ps_alloc_plasma_state: ierr=',ierr
               stop
            endif


c           Set nbi calc code name
            ps%nbi_code_info='cql3d'

c           Initialize arrays, 
            call zone_check(nrho,ps%rho_nbi,rya(1))  !rya starts at 0:...

            ps%curbeam=0.0_rspec    !dimension 1:nrho
            ps%eperp_beami=0.0_rspec
	    ps%epll_beami=0.0_rspec
	    ps%pbe=0.0_rspec    !Power to electrons
	    ps%pbi=0.0_rspec    !Power to ions
	    ps%pbth=0.0_rspec   !Power of beam thermalization
c           More to come, as get more specific with beam


         elseif(cql3d_output.eq.'RW'.or. cql3d_output.eq.'LH+RW' ) then

            nrho=lrz+1  !lrz is number of bin centers spec'd in cqlinput


            write(*,*)'nrho frm nml, ps%nrho_rw frm PS =',
     +	               nrho,ps%nrho_rw

c           BUNCH MORE CODING HERE [BH120724]
!           Necessary PS dimensions (user has to set his component dims):
            if(ps%nrho_rw.eq.nrho) then
               continue         ! dimension OK already
            else if(ps%nrho_rw.eq.0) then
               ps%nrho_lhrf=nrho
               call ps_alloc_plasma_state(ierr) !set these PS dims in PS
               write(*,*)'ps%nrho_nrho set= ',nrho
c           Following if clause should not be satisfied in ordinary IPS usage
c           (Can occur if using already initialized PS)
            else if(ps%nrho_rw.ne.nrho) then
               write(*,*) 
     +           ' * prepare_cql3d_input: reset PS RunaWay profile size'
               write(*,*) ' * from ',ps%nrho_rw,' to ',nrho
               
!              copy all EXCEPT RW component profiles
               cclist = ps_cc_all
               cclist(ps_cnum_RW)=0

               call ps_copy_plasma_state(ps, aux, ierr, cclist = cclist)
               if(ierr.eq.0) then
                  ! OK, copy back to ps
                  call ps_copy_plasma_state(aux, ps, ierr)
                  
                  if(ierr.eq.0) then
                     ! set desired dimension
                     ps%nrho_rw=nrho
                     call ps_alloc_plasma_state(ierr) !set these PS dims in PS
                  endif
               endif
            endif
            
            if (ierr.ne.0) then
               write(*,*)'prepare_cql3d_input:  RunaWay'//
     +                     ' ps_alloc_plasma_state: ierr=',ierr
               stop
            endif

c           Set runaway calc code name
            ps%runaway_code_info='cql3d'

c           Initialize arrays, 
            call zone_check(nrho,ps%rho_rw,rya(1))  !rya starts at 0:...

            ps%cur_rw=0.0_rspec    !dimension 1:nrho
            ps%dcur_rw_dvloop=0.0_rspec
            ps%eperp_rw=0.0_rspec
	    ps%epll_rw=0.0_rspec
	    ps%nrw=0.0_rspec

         endif  ! on cql3d_output

!     Stores the initialized plasma state, and write file (ps)
!     --but doesn't commit it.
!     The previous one (psp)is still around

        call ps_store_plasma_state(ierr)
        if(ierr .ne. 0) then
          write(iout,*)
     +    'Cannot ps_store_plasma_state in prepare_cql3d_input:init'
        else
	  write(iout,*) "init'd cur_state.cdf stored"
        end if


        go to 999  ! go to end
                        
      endif  ! on ips_mode.eq.'init'

!.......................................................................
!     Set cqlinput for initial step of restart
!.......................................................................
      write(*,*) 'Start initial step' 
    
      if (restart .eq. 'enabled') then
         nlrestrt='ncdfdist'
         nlwritf='ncdfdist'
      else
         nlrestrt='disabled'
         nlwritf='ncdfdist'
      endif

!.......................................................................
!     Set cqlinput nsteps, dtr and enorm
!.......................................................................

      write(*,*) 'nsteps_str = ', nsteps_str
      write(*,*) 'deltat_str = ', deltat_str

      read(nsteps_str, '(i10)') nsteps
! ptb:      read(nsteps_string, '(i10)') nsteps
      nstop=nsteps     !Reset these variables in namelist
      nplot=nsteps
      nplt3d=nsteps
      read(deltat_str, '(e14.7)') deltat
! ptb:      read(deltat_string, '(e15.7)') deltat
      dtr=deltat
      write(*,*) 'deltat from deltat_string = ', deltat  ! ptb:
      write(*,*) 'nsteps from nsteps_string = ', nsteps  ! ptb:

	 if (enorm_str_present .and. trim(enorm_str) /= 'None') then
        write(*,*) 'enorm_str = ', enorm_str
		read(enorm_str,*) enorm
	 end if

!.......................................................................
c     If only need electron profiles plus zeff (ngen=1 and 
c     iprozeff=spline from cql3d namelist), then specify
c     an ion species in cqlinput, and cql3d will add a further
c     high-Z impurity, if necesary.
c     (The cqlinput must be setup for electron gen species.
c      Some checking is done here.)  
c     Else, obtain number
c     of ions for the plasma state abridged species list.
!     See Plasma_State_V2.003.pdf (May 2008) re Abridged Species.
!.......................................................................

      write(*,*) 'Step at cql3d_mode'       

      if (cql3d_mode.eq.'el-only') then

      if (ngen.eq.1 .and. iprozeff.eq.'spline' .and. izeff.eq.'backgrnd' 
     +   .and. bnumb(1).eq.-1.d0) then 
         ! i.e.,  cqlinput set up for single electron gen species
         continue
      else
         write(*,*)'STOP: cqlinput not setup for electron gen species'
         write(*,*)'CHECK: ngen,iprozeff,izeff,bnumb(1)'
         stop
      endif

      elseif (cql3d_mode.eq.'el-ions' .or. cql3d_mode.eq.'ions-el') then

         ! Get number of ion species from PS abridged list
         ! BH111108         ntotal=ps%nspec_alla+1
!        cqlinput ngen=1 is assumed  for time being.  It is
!        expected the single gen species will be consistent
!        with the PS.  [Probably an ion in this case, but e might work.]
         nmax=ps%nspec_alla+1   !PS thermal species + electrons

c        Check/adjust izeff
         if (izeff.ne.'ion') then
            write(*,*)'Resetting izeff=ion for all ion spec-d cases'
            izeff='ion'
         endif
      endif

c---------------------------------------------------------------
c     Check/adjust plasma profile specs for passage from PS to
c     cqlinput (set cql3d namelist njene=101).
c---------------------------------------------------------------

      if (njene.ne.101) then
         write(*,*)'Resetting njene for uniform, 101-pt plasma profs'
         njene=101
         nj=njene    !local variable
         njp1=nj+1   !local variable
         if (njene.gt.njenea) then
            write(*,*)'Stop:  Need njenea bigger in param.h'
            stop
         endif
      endif


c     Make sure cql3d namelist density and temp scaling turned off,
c     except introduce any required unit changes here.
      enescal=1.d-6
      tescal=1.d0  !cqlinput is in keV.
      tiscal=1.d0
      zeffscal=1.d0
! This must be wrong because the value of elecfld in the cql3d output file is 
! too small by just this factor ! Also I notice that the default value of
! elecscal in aindflt.f is 1.0. So leave elecscal at 1.0.
! ptb:      elecscal=1./(2.*pi*ps%r_axis*1.d2)  !v/cm from loop volts
      vphiscal=ps%r_axis*1.d2  !Change to cm/sec from omega (rad/sec)
      ennscal=1.d-6   !Change to /cm**3


c     Reset cqlinput eqdskin according to PS
c From Francesca Poli, 131018: we had to change the [eqdskin] line 
c to make sure genray reads the correct eqdsk file.
c I was getting an error, because GENRAY was looking for a file 
c that was not there.
c The ps_eqdsk variable had in memory the name of the default plasma
c state file, which starts with $INPUT_SUFFIX, rather than the correct
c plasma state eqdsk file that start with $SIM_NAME.
c      eqdskin=trim(ps%eqdsk_file)
      eqdskin='eqdsk'
c     Reset cqlinput radcoord to ensure radial coord proportional
c     to sqrt(tor flux)
      radcoord='sqtorflx'

      

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


      if (PS_add_nml.eq.'disabled') then  !PS data inserted in cqlinput
                                          !and not other way around

!.......................................................................
!     Assume RF powers have been prevously set in a call to genray.
!     Alternatively, could set power in genray=1 MW for each
!     source, and renormalize here.   But, we will comment out
!     this possibility.
!     However, we use internal cql3d-freya code for NBI, and thus
!     set that power here.
!.......................................................................

            
c     Below, set some cqlinput data from the PS.
c     Add to this section, as required by details of the modeling.
c
c     Assume ray data in genray.nc, and that the genray
c     power is set to 1 MW.  Then adjust with cql3d nml pwrscale(1).
c     Check that no more than 1 nbi source (Else need to adjust code.)

         if (cql3d_output.eq.'EC') then
            
!             pwrscale(1)=ps%power_ec(1)*1.d-6  !Convert to Watts
             pwrscale(1)=1.d0

! ptb         elseif (rfmode.eq.'LH') then
         elseif (cql3d_output.eq.'LH') then
            
!              pwrscale(1)=ps%power_lh(1)*1.d-6
              pwrscale(1)=1.d0

! ptb         elseif (rfmode.eq.'LH+RW') then
         elseif (cql3d_output.eq.'LH+RW') then
            
!              pwrscale(1)=ps%power_lh(1)*1.d-6
              pwrscale(1)=1.d0

! ptb         elseif (rfmode.eq.'IC') then
         elseif (cql3d_output.eq.'IC') then
            
!              pwrscale(1)=ps%power_ic(1)*1.d-6
              pwrscale(1)=1.d0

! ptb         elseif (rfmode.eq.'NBI') then
         elseif (cql3d_output.eq.'NBI') then
            
              if (ps%nbeam.ne.1) then
                 write(*,*)'STOP:  Need to adjust prepare_cql3d_input'
                 write(*,*)'       and PS for nbeam.gt.1'
                 stop
              endif
              bptor(1)=ps%power_nbi(1)

! ptb         elseif (rfmode.eq.'NBI+IC') then
         elseif (cql3d_output.eq.'NBI+IC') then
 
              if (ps%nbeam.ne.1) then
                 write(*,*)'STOP:  Need to adjust prepare_cql3d_input'
                 write(*,*)'       for nbeam.gt.1'
                 stop
              endif
              bptor(1)=ps%power_nbi(1)
!              pwrscale(1)=ps%power_ic(1)*1.d-6
              pwrscale(1)=1.d0

         else
            write(*,*)'Incorrect command line cql3d_output'
            STOP
         endif

         
      endif                     ! on PS_add_nml.eq.'disabled'
      


!!$!.......................................................................
!!$!     Usual operation: PS_add_nml will be 'disabled', and cqlinput data
!!$!     will be set from the PS data, as above.
!!$!     However,
!!$!     check whether plasma state RF variables are defined.
!!$!     If not, and command line variable PS_add_nml.eq.'enabled',
!!$!     then will add data to the PS from the genraynml data.
!!$!     If RF variables are defined or not, but PS_add_nml.eq.'force',
!!$!     the values from genraynml data will be inserted in the PS.
!!$!.......................................................................
!!$
!!$      if (PS_add_nml.eq.'enabled' .or. PS_add_nml.eq.'force') then
!!$                                      !Some genraynml data added to PS
!!$
!!$         i_check_rf=0           !Set =1 if PS and genraynml agree
!!$                                !on number of sources.
!!$         
!!$         if (rfmode.eq.'EC') then
!!$            nxxrf_src=ps%necrf_src
!!$            if (ngenray_source.eq.nxxrf_src) i_check_rf=1
!!$            
!!$         elseif (rfmode.eq.'LH') then
!!$            nxxrf_src=ps%nlhrf_src
!!$            if (ngenray_source.eq.nxxrf_src) i_check_rf=1
!!$            
!!$         elseif (rfmode.eq.'IC') then
!!$            nxxrf_src=ps%nicrf_src
!!$            if (ngenray_source.eq.nxxrf_src) i_check_rf=1
!!$            
!!$         else
!!$            write(*,*)'Incorrect command line rfmode spec'
!!$            STOP
!!$         endif
!!$
!!$         write(*,*)'i_check_rf=',i_check_rf
!!$         if (i_check_rf.eq.0 .and. PS_add_nml.eq.'enabled') then
!!$            write(*,*)'prepare_genray_input: genraynml and PS disagree'
!!$            write(*,*)'                      on number of rfsources, '
!!$            write(*,*)'                      and ps_add_nml.eq.enabled'
!!$            write(*,*)'                      ngenray_source=',ngenray_source
!!$            write(*,*)'                      ps%nxxrf_src=',nxxrf_src
!!$            write(*,*)'                      Need ps_add_nml.eq.force'
!!$            STOP
!!$         endif
!!$
!!$c     For genray.in and PS agree on number of sources, but for
!!$c     frequency not set, then put in in PS from genray.in.
!!$         if (i_check_rf.eq.1  .and. PS_add_nml.eq.'enabled')  then
!!$            if (rfmode.eq.'EC') then
!!$               do is=1,ngenray_source
!!$                   if(ps%freq_ec(isource+(is-1)).eq.0.0_rspec)
!!$     +                ps%freq_ec(isource+(is-1))=frqncy
!!$               enddo
!!$            elseif (rfmode.eq.'LH') then
!!$               do is=1,ngenray_source
!!$                   if(ps%freq_lh(isource+(is-1)).eq.0.0_rspec)
!!$     +                ps%freq_lh(isource+(is-1))=frqncy
!!$               enddo
!!$            elseif (rfmode.eq.'IC') then
!!$               do is=1,ngenray_source
!!$                   if(ps%freq_ic(isource+(is-1)).eq.0.0_rspec)
!!$     +                ps%freq_ic(isource+(is-1))=frqncy
!!$               enddo
!!$            endif
!!$         endif
!!$
!!$
!!$c     For number of PS sources .ne. genraynml sources:
!!$c     if  PS_add_nml.eq.'force',  modify PS to agree, and set PS power
!!$c     freq_xx, number of sources, and source names.
!!$
!!$         if (i_check_rf.eq.0 .and. PS_add_nml.eq.'force') then
!!$            
!!$            if (rfmode.eq.'EC') then
!!$               if (ps%necrf_src.eq.0) then
!!$                  ps%necrf_src=ngenray_source
!!$                  call ps_alloc_plasma_state(ierr) !set PS dims for rf
!!$               else
!!$                  write(*,*) 
!!$     +              'prepare_genray_input: resetting PS EC source count'
!!$                  write(*,*) ' from ',ps%necrf_src,' to ',ngenray_source
!!$                  
!!$!     copy all EXCEPT EC component profiles
!!$                  cclist = ps_cc_all
!!$                  cclist(ps_cnum_EC)=0
!!$                  
!!$                  call ps_copy_plasma_state(ps,aux,ierr,cclist= cclist)
!!$                  if(ierr.eq.0) then
!!$!     OK, copy back to ps
!!$                     call ps_copy_plasma_state(aux, ps, ierr)
!!$                     
!!$!                     if(ierr.eq.0) then
!!$!     set desired dimension
!!$                        ps%necrf_src=ngenray_source
!!$                        call ps_alloc_plasma_state(ierr) !set dims in PS
!!$!                     endif      !on ierr
!!$                  endif         !on ierr
!!$                  
!!$               endif            !on ps%necrf_src
!!$
!!$c     Set EC PS values from genraynml[NOTE: only 1 frqncy]
!!$
!!$!     Reminder of data to be passed to PS:
!!$!     EC   necrf_src=   number of ECRF sources
!!$!          power_ec(necrf_src)= power on each ECRF source (W)
!!$!          freq_ec(necrf_src)=frequency on each ECRF source (Hz)
!!$!
!!$!     LH   nlhrf_src= number of LHRF sources
!!$!          power_lh(nlhrf_src)= power on each LHRF source (W)
!!$!          freq_lh(nlhrf_src)= frequency on each LHRF source (Hz)
!!$!
!!$!     IC   nicrf_src= number of ICRF sources (includes FW, HHFW)
!!$!          power_ic(nicrf_src)= power on each ICRF source
!!$!          freq_ic(nicrf_src)= frequency on each ICRF source (Hz)
!!$!
!!$               do i=1,ngenray_source
!!$                  ps%ecrf_src_name=achar(48+i) !achar(48)='0'
!!$                  if (istart.eq.1) ps%power_ec(i)=powtot(i)
!!$                  if (istart.gt.1) ps%power_ec(i)=powers(i)
!!$ !MKSA issue
!!$                  ps%freq_ec(i)=frqncy  !PS in Hz
!!$               enddo
!!$                  
!!$               
!!$            elseif (rfmode.eq.'LH') then
!!$               ps%lhrf_src_name=achar(48+i) !achar(48)='0'
!!$               if (ps%nlhrf_src.eq.0) then
!!$                  ps%nlhrf_src=ngenray_source
!!$                  call ps_alloc_plasma_state(ierr) !set PS dims for rf
!!$               else
!!$                  write(*,*) 
!!$     +              'prepare_genray_input: resetting PS LH source count'
!!$                  write(*,*) 'from ',ps%nlhrf_src,' to ',ngenray_source
!!$                  
!!$!     copy all EXCEPT EC component profiles
!!$                  cclist = ps_cc_all
!!$                  cclist(ps_cnum_LH)=0
!!$                  
!!$                  call ps_copy_plasma_state(ps,aux,ierr,cclist= cclist)
!!$                  if(ierr.eq.0) then
!!$!     OK, copy back to ps
!!$                     call ps_copy_plasma_state(aux, ps, ierr)
!!$                     
!!$!                     if(ierr.eq.0) then
!!$!     set desired dimension
!!$                        ps%nlhrf_src=ngenray_source
!!$                        call ps_alloc_plasma_state(ierr) !set dims in PS
!!$!                     endif !on ierr
!!$                  endif !on ierr
!!$
!!$               endif !on ps%nlhrf_src
!!$
!!$               do i=1,ngenray_source
!!$                  ps%lhrf_src_name=achar(48+i)  !achar(48)='0'
!!$                  if (istart.eq.1) ps%power_lh(i)=powtot(i)
!!$                  if (istart.gt.1) ps%power_lh(i)=powers(i)
!!$!MKSA issue
!!$                  ps%freq_lh(i)=frqncy
!!$               enddo
!!$                  
!!$            elseif (rfmode.eq.'IC') then
!!$               if (ps%nicrf_src.eq.0) then
!!$                  ps%nicrf_src=ngenray_source
!!$                  call ps_alloc_plasma_state(ierr) !set PS dims
!!$               else
!!$                  write(*,*) 'prepare_genray_input: '//
!!$     +                 'resetting PS IC source count'
!!$                  write(*,*) 'from ',ps%nicrf_src,
!!$     +                 ' to ',ngenray_source
!!$                  
!!$!     copy all EXCEPT EC component profiles
!!$                  cclist = ps_cc_all
!!$                  cclist(ps_cnum_IC)=0
!!$                  
!!$                  call ps_copy_plasma_state
!!$     +                 (ps,aux,ierr,cclist=cclist)
!!$                  if(ierr.eq.0) then
!!$!     OK, copy back to ps
!!$                     call ps_copy_plasma_state(aux, ps, ierr)
!!$                     
!!$!                     if(ierr.eq.0) then
!!$!     set desired dimension
!!$                        ps%nicrf_src=ngenray_source
!!$                        call ps_alloc_plasma_state(ierr) !set dims in PS
!!$!                     endif      ! on ierr
!!$                  endif         ! on ierr
!!$                  
!!$               endif            ! on ps%nicrf_src
!!$
!!$               do i=1,ngenray_source
!!$                  ps%icrf_src_name=achar(48+i)  !achar(48)='0'
!!$                  write(*,*)'icrf_src_name=',achar(48+i)
!!$                  if (istart.eq.1) ps%power_ic(i)=powtot(i)
!!$                  if (istart.gt.1) ps%power_ic(i)=powers(i)
!!$!MKSA issue
!!$                  ps%freq_ic(i)=frqncy
!!$               enddo
!!$               
!!$            endif     ! on rfmode
!!$            
!!$         endif        ! on i_check_rf.eq.0 .and. PS_add_nml.eq.'force'
!!$         
!!$      endif           ! on (PS_add_nml.eq.'enabled' .or. PS_add_nml.eq.'force
!!$                     
!!$
!!$!.......................................................................
!!$!     In general, as many genray namelist parameters as possible will be
!!$!     specified only through the genray.in/.dat file, usually set up
!!$!     for a particular rf mode and machine.  Depending on these
!!$!     settings, additional data will be obtained from the plasma state.
!!$!.......................................................................
!!$


!.......................................................................
c     Interpolate PS profiles to user (i.e., cqlinput) ryain-grid:
c     Use ps_rho_rezone rezonning subroutine from PS, onto 
c     rho_bdy_rezone(...) which gives PS values at the ryain()
c     bin boundaries (except values centered in last half bin at
c     plasma center and edge are taken as ryain()=0. and 1. values).
c
c     That is, size(rho_bdy_rezon) = size(ryain) + 1
c     This gives size(ryain) bins with size(ryain)+1 bin boundaries.
c     Plasma profile values at the centers of the bins are taken 
c     to be values at the user ryain-grid bin boundaries.
c
c     Alternatively, could use ps_interp_1d [cf Bonoli, toric prep code]
!.......................................................................

c     Equispaced, bin-boundary grid for profiles to cql3d

      dryain=1.0_rspec/(nj-1)
      do l=1,nj
         ryain(l)=(l-1)*dryain
      enddo

      allocate(rho_bdy_rezon(njp1))
      rho_bdy_rezon(1)=ryain(1)
      rho_bdy_rezon(njp1)=ryain(nj)
      rho_bdy_rezon(2:nj) =
     +     .5_rspec*(ryain(1:nj-1)+ryain(2:nj))

!.......................................................................
!     cql3d_mode: 'el-only' 'el+ions' 'ions+el'
c     First treat electrons-only and zeff case:
c     See Plasma_State2_V0.pdf, fr. 18, for ps_rho_rezone.
!.......................................................................

      if (cql3d_mode.eq.'el-only') then
c     Choose ion species to be highest central dens of abrdgd ions (SA)
      iden_sa=0
!ptb      den_sa_max=0d0
      den_sa_max=0.0_rspec
      do i=1,ps%nspec_tha
! ptb hacking
! It turns out that care must be taken when using ps%sa_index(j) since 
! ps%sa_index(j)=-1 indicates that element of (j) is derived from a composite
! of ions in the "S" list using Zeff and quasinetrality.
!         ii=ps%sa_index(i)
         if (ps%sa_index(i) .eq. -1) then
            ii = i
         else
            ii = ps%sa_index(i)
         endif
!end of ptb hacking
         if (ps%ns(1,ii).gt.den_sa_max) then
            den_sa_max=ps%ns(1,ii)
            iden_sa=ii
         endif
      enddo
      write(*,*)'Main thermal ion species, name_sa, iden_sa ptr in S=',
     +          trim(ps%sa_name(iden_sa)), iden_sa

c     Set cql3d species parameters:

      ngen=1 !Can be increased in future, as needed.
      nmax=2 !Will be increased by 1 in cql3d in accord with zeff.

      bnumb(1)=-1.d0
      bnumb(2)=ps%q_sa(iden_sa)/xe
      bnumb(3)=-1.d0
! ptb: clq3d namelist description says that fmass should be in cgs - g. 
! ptb: But me and ps%m_SA are in kG. So they should be multiplied by 10^3
! ptb: to convert to g (not multiplied by 10^-3 !!!).
! ptb:      fmass(1)=me*1.d-3
! ptb:      fmass(2)=ps%m_sa(iden_sa)*1.d-3
      fmass(1)=me*1.d+3
      fmass(2)=ps%m_sa(iden_sa)*1.d+3
      fmass(3)=fmass(1)
      kspeci(1,1)='e'
      kspeci(2,1)='general'
      kspeci(1,2)=trim(ps%sa_name(iden_sa))
      kspeci(2,2)='maxwell'
      kspeci(1,3)='e'
      kspeci(2,3)='maxwell'


c     Check set spline profiles are expected for electrons, ions,
c     and zeff.
      if (iprone.ne.'spline' .or. iprote.ne.'spline' .or.
     +    iproti.ne.'spline' .or. iprozeff.ne.'spline') then
          write(*,*)'STOP: Need iprone/iprote/iproti/iprozeff = spline'
          stop
      endif
      if (iproelec.eq.'spline') then
         efswtch='method1'  !To ensure simple spline input
      else
         write(*,*)'iproelec.ne.spline:  Set tor E-fld parabolically'
         write(*,*)'(Can set elecfld(0:1)=0.)'
         iproelec='parabola'
         efswtch='method1'
      endif

      ntotal=ngen+nmax
      kelecg=1
      kelecm=3

!!$      Don't need to do this, since simply setting kelecg/kelecm values
!!$      At any rate, these are not namelist variables.
!!$      do k=1,ngen
!!$        if (fmass(k) .lt. 1.e-27) then
!!$          kelecg=k
!!$          goto 2001
!!$        endif
!!$      continue
!!$ 2001 continue
!!$      kelecm=0
!!$
!!$      do  k=ngen+1,ntotal
!!$        if (fmass(k) .lt. 1.e-27) then
!!$          kelecm=k
!!$          goto 2003
!!$        endif
!!$      enddo
!!$ 2003 continue
!!$
!!$c     What is this?
!!$      do k=1,ntotal
!!$         if (k.eq.kelecg .or. k.eq.kelecm) then
!!$            incl_species(0)=1   !Using PS index convention here.
!!$         endif
!!$      enddo

!.......................................................................
!     Fill in  electron radial profiles for cqlinput file
!.......................................................................
     
      do k=1,ntotal
         if (k.eq.kelecg .or. k.eq.kelecm) then
!            call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(0), enein(1,k),
! ptb            call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(0), enein(:,k),
! ptb     +                         ierr, zonesmoo=.TRUE.)
            call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(0), enein(1:nj,k),
     +                         ierr, zonesmoo=.TRUE.)
            write(*,*)'after ps_rho_rezone: ierr=',ierr
            call ckerr('ps_rho_rezone (U1)',ierr)
         endif
      enddo

! ptb      call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(0), tein(:), ierr,
! ptb     +     zonesmoo=.TRUE.)
      call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(0), tein(1:nj), ierr,
     +     zonesmoo=.TRUE.)
      call ckerr('ps_rho_rezone (U2)',ierr)

c     iden_sa species (above) is main ion species:
! ptb      call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(iden_sa), tiin(:), ierr,
! ptb     +     zonesmoo=.TRUE.)
      call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(iden_sa), tiin(1:nj), ierr,
     +     zonesmoo=.TRUE.)
      call ckerr('ps_rho_rezone (U2)',ierr)

! ptb      call ps_rho_rezone(rho_bdy_rezon, ps%id_Zeff, zeffin(:), ierr,
! ptb     +     zonesmoo=.TRUE.)
      call ps_rho_rezone(rho_bdy_rezon, ps%id_Zeff, zeffin(1:nj), ierr,
     +     zonesmoo=.TRUE.)
      call ckerr('ps_rho_rezone (U2)',ierr)

      do l=1, ps%nrho-1
         if (zeffin(l).lt.1.0) then
            write(*,*)'Problem: PS gives zeffin.lt.1: l, zeffin(l)',
     +                l, zeffin(l)
            stop
         endif
      enddo

! ptb: added this coding to form elecfld (V/cm) from data in the PS in terms 
! ptb: of the loop voltage in the PS (ps%vsur) - note sign convention:
! ptb: Positive (CCW) electric field pushes electrons in the negative (CW) 
! ptb  direction, so that the ohmic current is negative - which is the case 
! ptb  in C-Mod.
!      elecfld(0) = +ps%vsur / (2.*pi*ps%r_axis*1.0d+02)
!      elecfld(1) = +ps%vsur / (2.*pi*ps%r_axis*1.0d+02)
!      mpwrelec = 1.0_rspec
!      npwrelec = 1.0_rspec
! Realized that we need the actual electric field as a function radius from TSC
! and the Plasma State File
      iproelec = 'spline'
      call ps_rho_rezone(rho_bdy_rezon, ps%id_V_loop, elecin(1:nj), ierr,
     +     zonesmoo=.TRUE.)
      call ckerr('ps_rho_rezone (U2)',ierr)
      elecin(1:nj) = elecin(1:nj) / (2.*pi*ps%r_axis*1.0d+02) 
      write(*,*) 'elecin(i)', elecin(1:nj)
c     Checking interpolation
      write(*,*)
      write(*,*)'Checking interpolation'
      write(*,*)'  l  ryain  rho_bdy_rezon   enein   tein  tiin  zeffin'
      do l=1,nj
      write(*,1001)l,ryain(l),rho_bdy_rezon(l),enein(l,1),tein(l),
     +             tiin(l),zeffin(l)
      enddo
 1001 format(i4,6(1pe11.3))

      endif  !On cql3d_mode.eq.'el-only'


!.......................................................................
!     For additional interfacing of to cql3d, the alla species
!     (electrons, abridged thermal ions [fully stripped light ions
!     and two impuritied], and the nonthermal species [NBI, Minority
!     ion species for RF, and fusion products]) are brought in.
!     These can then be rearranged as required.
!
!     Use Doug McCune coding from his plasma_state_test.f90, to get
!     density/temperatures, and augment with rotation rate.
!     The interpolation is on to the above rho_bdy_rezon, appropriate
!     for cql3d namelist input.
!.......................................................................

      if (cql3d_mode .ne. 'el-only') then

  !  May 2008: demonstrate  (DMC, in relation to BH genray work).
  !  creation of single arrays of all the densities and temperatures
  !  for the nspec_alla list:
  !       abridged thermal ion species + all fast species
  !
  !  rezoned to a boundary oriented grid; use rho_bdy_rezone(...) for
  !  this rezoning.
  !

!      nrho_user=201
!      nrhop1=nrho_user+1
      nrho_user=nj
      nrhop1=nj+1

!  allocate(rho_user(nrho_user),rho_bdy_rezon(nrhop1))
!
!  do ii=1,nrho_user
!     rho_user(ii) = (ii-1)*1.0_rspec/(nrho_user-1)
!  enddo
      allocate(rho_user(nrho_user))
      rho_user(:)=ryain(:)

      !  rho_bdy_rezon(...) has a half width zone at the center and edge, and a
      !  full width zone around each interior boundary in rho_user(...).
      !  size(rho_bdy_rezon) = size(rho_user) + 1

      rho_bdy_rezon(1)=rho_user(1)
      rho_bdy_rezon(nrhop1)=rho_user(nrho_user)
      rho_bdy_rezon(2:nrho_user) =
     +      0.5_rspec*(rho_user(1:nrho_user-1)+rho_user(2:nrho_user))

      allocate(dense(nrho_user),tempe(nrho_user))
      allocate(densi(nrho_user,ps%nspec_alla),tempi(nrho_user,ps%nspec_alla))
      
      call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(0), dense, ierr, zonesmoo=.TRUE.)
      call ckerr('ps_rho_rezone (U1)')
      
      call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(0), tempe, ierr, zonesmoo=.TRUE.)
      call ckerr('ps_rho_rezone (U2)')
      
      !  get the abridged thermal species list first...
      
      call ps_tha_rezone(rho_bdy_rezon, ierr, zonesmoo=.TRUE.,
     +     ns = densi(:,1:ps%nspec_tha), Ts = tempi(:,1:ps%nspec_tha))
      call ckerr('ps_tha_rezone')
      
      !  get the fast species; in general these will be on different grids.
      
      allocate(wk_dens(nrho_user),wk_eperp(nrho_user),wk_epll(nrho_user))
      
      do ii=ps%nspec_tha + 1, ps%nspec_alla
         jj = ps%alla_index(ii)   ! index in "all species" list
         !       jj=-1 if ps%alla_type(ii).eq.ps_impurity, but this should not
         !       occur as here we loop over fast ions only.  Still, I test this...
         if(jj.le.0) then
            write(iout,*) ' ? unexpected alla_index(ii) value: ',jj
         else
            !  OK, jj=alla_index(ii) points to the element index of whatever
            !  fast species list element corresponds to element ii of the alla list
            
            if(ps%alla_type(ii).eq.ps_beam_ion) then
               call ps_rho_rezone(rho_bdy_rezon, ps%id_nbeami(jj), wk_dens,
     +              ier1, zonesmoo=.TRUE.)
               call ps_rho_rezone(rho_bdy_rezon, ps%id_eperp_beami(jj), wk_eperp,
     +              ier2, zonesmoo=.TRUE.)
               call ps_rho_rezone(rho_bdy_rezon, ps%id_epll_beami(jj), wk_epll,
     +              ier3, zonesmoo=.TRUE.)
               
            else if(ps%alla_type(ii).eq.ps_rf_minority) then
               call ps_rho_rezone(rho_bdy_rezon, ps%id_nmini(jj), wk_dens,
     +              ier1, zonesmoo=.TRUE.)
               call ps_rho_rezone(rho_bdy_rezon, ps%id_eperp_mini(jj), wk_eperp,
     +              ier2, zonesmoo=.TRUE.)
               call ps_rho_rezone(rho_bdy_rezon, ps%id_epll_mini(jj), wk_epll,
     +              ier3, zonesmoo=.TRUE.)
               
            else if(ps%alla_type(ii).eq.ps_fusion_ion) then
               call ps_rho_rezone(rho_bdy_rezon, ps%id_nfusi(jj), wk_dens,
     +              ier1, zonesmoo=.TRUE.)
               call ps_rho_rezone(rho_bdy_rezon, ps%id_eperp_fusi(jj), wk_eperp,
     +              ier2, zonesmoo=.TRUE.)
               call ps_rho_rezone(rho_bdy_rezon, ps%id_epll_fusi(jj), wk_epll,
     +              ier3, zonesmoo=.TRUE.)
               
            endif
            
            ierr = ier1 + ier2 + ier3
            call ckerr('ps_rho_rezone (fast specie)')
            
            densi(:,ii) = wk_dens
            tempi(:,ii) = (2.0_rspec/3.0_rspec)*(wk_eperp + wk_epll)
            
         endif
      enddo  !On ii=ps%nspec_tha + 1, ps%nspec_alla
      
      write(iout,'(A,1pe12.5)')
     +    '  density      Temperature   Species @ rho_user(7) = ',
     +    rho_user(7)
      write(iout,'(2(1x,1pe12.5),3x,a)') dense(7),tempe(7),
     +     trim(ps%alla_name(0))
      do ii=1,ps%nspec_alla
         write(iout,'(2(1x,1pe12.5),3x,a)') densi(7,ii),tempi(7,ii),
     +        trim(ps%alla_name(ii))
      enddo

      call ps_rho_rezone(rho_bdy_rezon, ps%id_Zeff, zeffin(:), ierr,
     +     zonesmoo=.TRUE.)
      call ckerr('ps_rho_rezone (U2)',ierr)

!     Check ion densities/temperatures
         do ii=1,ps%nspec_alla
            do j=1,nj-1
               if(abs(densi(j,ii)).lt.eps_dens) incl_species(ii)=0 
!!$               if(abs(tempi(j,ii)).lt.eps_temp) incl_species(ii)=0
            enddo
!!$            if(abs(tempi(nj,ii)).lt.eps_temp) incl_species(ii)=0
         enddo

!!$         write(*,*)'incl_species(0:ps%nspec_alla)=',
!!$     +              incl_species(0:ps%nspec_alla)
         write(*,*)'incl_species(0:ps%nspec_tha)=',
     +              incl_species(0:ps%nspec_tha)


      endif  !cql3d_mode.ne.'el-only'



      if (cql3d_mode .eq. 'el+ions') then
      
!.......................................................................
!     Electrons are assumed to be the general cql3d (first) species.
!     Will put the Maxwln electrons at the end of the species list.
!     Ion Maxwln species are species 1:nspec_tha (abridged species list).
!     CQL3D does not treat separate thermal and fast ion species
!     such as from NBI or fusion;  the fast+thermal ions are one
!     non-Maxwln species.
!     In processing cql3d output:
!     it would be possible, to separate the "thermal" and
!     "nonthermal" species, say by fitting a Maxwln to the bulk,
!     or, to divide at some appropriate velocity ~2-3 vth.
!
!     For the time being, will only set up electrons (first general
!     species) and thermal ions, and thermal equivalents of the the
!     fast ion species.   This will provide an alternative
!     treatment of above 'el-only' case, and in cases where the 
!     wave particle interaction is with electrons-only, and no
!     nonthermal ions, it should give same results.
!
!     If there is a nonthermal species, this can affect Zeff;
!     we thus include them in to the following.
!
!     CQL3D spline fits use a single "background" ion temperature.
!     We use a PS density weighted average of the nspec_tha ions,
!     obtained with plasma state module routine ps_Ti_rezone.
!     CQL3D can be readily adapted to multiple ion temperature input,
!     if there is need.  For example, it might be useful to provide
!     collisional power going to each of the ion species.
!
!     Ref: Plasma_State_V2.003.pdf (May 2008) re Abridged Species.
!.......................................................................


c     Set cql3d species parameters:

      ngen=1 !Can be increased in future, as needed
      nmax=ps%nspec_alla !ions (abridged and fast) plus maxwl electrons

      bnumb(1)=-1.d0
      kspeci(1,1)='e'
      kspeci(2,1)='general'
! ptb: clq3d namelist description says that fmass should be in cgs - g. 
! ptb: But me and ps%m_SA are in kG. So they should be multiplied by 10^3
! ptb: to convert to g (not multiplied by 10^-3 !!!).
      fmass(1)=me*1.d+3
! ptb:      fmass(1)=me*1.d-3
      do i=2,ps%nspec_alla
         bnumb(i)=ps%q_alla(i)/xe
! ptb         fmass(i)=ps%m_alla(i)*1.d-3
         fmass(i)=ps%m_alla(i)*1.d+3
         kspeci(1,2)=trim(ps%alla_name(i))
         kspeci(2,2)='maxwell'
      enddo
      ntotal=ngen+nmax

      !Add in electron maxwl species after maxwl ions
      bnumb(ntotal)=-1.d0
      fmass(ntotal)=fmass(1)
      kspeci(1,ntotal)='e'
      kspeci(2,ntotal)='maxwell'

      iprone='spline'
      iprote='spline'
      iproti='spline'  !Will ion density weighted value, from PS
      iprozeff='disabled'  !Will be calculated from given ions
      iproelec='spline'
      efswtch='method1'
      iprovphi='spline'  !Will use ion density weighted value, from PS

c     Set spline profiles expected for electrons, ions, efld, vphi

      
      
      
!BH120807:  Replace this section with DMC coding, as above
!BH120807:  We may need some following commented coding to eliminate
!BH120807:  unallocated fast species, so keep it for time being.

!!$       allocate(densi(nj,0:ps%nspec_alla+1),tempi(nj)) !densi incls elec
!!$       densi(:,:)=0.0_rspec
!!$       tempi(:,:)=0.0_rspec
!!$
!!$!     get electrons
!!$
!!$       call ps_rho_rezone(rho_bdy_rezon, ps%id_ns(0), densi(:,0),
!!$     +                         ierr, zonesmoo=.TRUE.)
!!$       write(*,*)'after ps_rho_rezone: ierr=',ierr
!!$       call ckerr('ps_rho_rezone (U1)',ierr)
!!$
!!$       call ps_rho_rezone(rho_bdy_rezon, ps%id_Ts(0), tein(:),
!!$     +                         ierr, zonesmoo=.TRUE.)
!!$       call ckerr('ps_rho_rezone (U2)',ierr)
!!$
!!$       enein(:,ntotal)=enein(:,1)
!!$            
!!$         
!!$!     get the alla ion species (abridged thermal + fast ions)
!!$         
!!$         call ps_tha_rezone(rho_bdy_rezon, ierr, zonesmoo=.TRUE.,
!!$     +       ns = densi(1:nj,1:ps%nspec_alla))
!!$         call ckerr('ps_tha_rezone',ierr)
!!$
!!$!     get ni-wtd average temp and tor rotation for thermal ion species
!!$
!!$         call ps_ti_rezone(rho_bdy_rezon, ierr, zonesmoo=.TRUE.,
!!$     +                     Ti=tempi(:), omegat=vphi(:))
!!$
!!$         incl_species(1:ps%nspec_tha)=1
         
!BH120802:  Skip nonthermal species for now.  Need to decide how
!BH120802:  how to best interface cql3d with them.
!!$
!!$!     get the fast species; in general these will be on different grids.
!!$         
!!$         allocate(wk_dens(nj),wk_eperp(nj),wk_epll(nj))
!!$         
!!$         do ii=ps%nspec_tha + 1, ps%nspec_alla
!!$            jj = ps%alla_index(ii) ! index in "all species" list
!!$.
!!$            write(*,*)'ii,jj,ps%nspec_tha,ps%nspec_alla',
!!$     +                 ii,jj,ps%nspec_tha,ps%nspec_alla
!!$
!!$            write(*,*)'ps%alla_type(0:ps%nspec_alla),ps_beam_ion',
!!$     +                 ps%alla_type(0:ps%nspec_alla),ps_beam_ion
!!$
!!$!     jj=-1 if ps%alla_type(ii).eq.ps_impurity, but this 
!!$!     should not occur as here we loop over fast ions only.
!!$!     Still, I test this...
!!$            if(jj.le.0) then
!!$               write(*,*) '? unexpected alla_index(ii) value: ',jj
!!$            else
!!$!     OK, jj=alla_index(ii) points to the element index of whatever
!!$!     fast species list element corresponds to element ii of the alla list
!!$               wk_dens(:)=0.0_rspec
!!$               wk_eperp(:)=0.0_rspec
!!$               wk_epll(:)=0.0_rspec
!!$
!!$               if(ps%alla_type(ii).eq.ps_beam_ion) then
!!$ 
!!$                if (allocated(nbeami)) then
!!$                incl_species(ii)=1
!!$                call ps_rho_rezone(rho_bdy_rezon, ps%id_nbeami(jj),
!!$     +                 wk_dens, ier1, zonesmoo=.TRUE.)
!!$                call ckerr('ps_rho_rezone: ier1',ier1)
!!$                call ps_rho_rezone(rho_bdy_rezon, ps%id_eperp_beami(jj),
!!$     +                 wk_eperp, ier2, zonesmoo=.TRUE.)
!!$                call ps_rho_rezone(rho_bdy_rezon, ps%id_epll_beami(jj),
!!$     +                 wk_epll, ier3, zonesmoo=.TRUE.)
!!$                endif
!!$              
!!$                  
!!$               else if(ps%alla_type(ii).eq.ps_rf_minority) then
!!$                if (allocated(nmini)) then
!!$                incl_species(ii)=1
!!$                call ps_rho_rezone(rho_bdy_rezon, ps%id_nmini(jj),
!!$     +                 wk_dens, ier1, zonesmoo=.TRUE.)
!!$                call ps_rho_rezone(rho_bdy_rezon, ps%id_eperp_mini(jj),
!!$     +                 wk_eperp, ier2, zonesmoo=.TRUE.)
!!$                call ps_rho_rezone(rho_bdy_rezon, ps%id_epll_mini(jj),
!!$     +                 wk_epll, ier3, zonesmoo=.TRUE.)
!!$                endif
!!$                  
!!$               else if(ps%alla_type(ii).eq.ps_fusion_ion) then
!!$                if (allocated(nfusi)) then
!!$                incl_species(ii)=1
!!$                call ps_rho_rezone(rho_bdy_rezon, ps%id_nfusi(jj),
!!$     +                 wk_dens, ier1, zonesmoo=.TRUE.)
!!$                call ps_rho_rezone(rho_bdy_rezon, ps%id_eperp_fusi(jj),
!!$     +                 wk_eperp, ier2, zonesmoo=.TRUE.)
!!$                call ps_rho_rezone(rho_bdy_rezon, ps%id_epll_fusi(jj),
!!$     +                 wk_epll, ier3, zonesmoo=.TRUE.)
!!$                endif
!!$                  
!!$               endif  !  On ps%alla_type(ii)
!!$               
!!$               ierr = ier1 + ier2 + ier3
!!$               call ckerr('ps_rho_rezone (fast specie)',ierr)
!!$               
!!$               densi(:,ii) = wk_dens
!!$
!!$!     McCune rational for following is that wk_perp and wk_epll
!!$!     are average perp and parallel energies.  (For isotropic
!!$!     distns, wk_eperp will be 2.*wkpell.)
!!$
!!$               tempi(:,ii) = (2.0_rspec/3.0_rspec)*(wk_eperp + wk_epll)
!!$               
!!$            endif  !  On jj
!!$         enddo  !  On ii
!!$
!!$         write(*,*)'incl_species(0:ps%nspec_all)=',
!!$     +              incl_species(0:ps%nspec_all)
!!$
!!$!     Determine reduced species list if any densities are ~zero inside
!!$!     LCFS, or temperatures zero, including LCFS.  
!!$!     Indicator is stored in incl_species(0:ps%nspec_alla).
!!$
!!$         eps_dens=1.e-6_rspec
!!$         eps_temp=1.e-6_rspec
!!$         write(*,*)'eps_dens,eps_temp=',eps_dens,eps_temp
!!$
!!$!     Check electron density/temperature
!!$         do j=1,nj-1
!!$            if(abs(densi(j,0)).lt.eps_den) incl_species(0)=0
!!$            if(abs(tempe(j)).lt.eps_temp) incl_species(0)=0
!!$         enddo
!!$         if(abs(tempe(nj)).lt.eps_temp) incl_species(0)=0
!!$
!!$         if (incl_species(0).eq.0) then
!!$            write(*,*)
!!$            write(*,*)'Electron density or temperature ~zero'
!!$            write(*,*)"Can't FP. Will stop"
!!$            write(*,*)
!!$            STOP
!!$         endif


!!$         jj=1  ! genray index for electrons
!!$         do ii=1,ps%nspec_alla
!!$         do ii=1,ps%nspec_tha
!!$            if (incl_species(ii).ne.0) then
!!$               jj=jj+1
!!$               charge_nc(jj)=abs(ps%q_alla(ii)/xe)
!!$               dmas_nc(jj)=ps%m_alla(ii)/me
!!$               do j=1,nj
!!$                  en_nc(j,jj)=densi(j,ii)
!!$                  temp_nc(j,jj)=tempi(j,ii)
!!$               enddo ! On j
!!$            endif
!!$         enddo ! On ii
!!$         nspecgr=jj
!!$         nbulk=nspecgr
!!$         ndens=nj

         write(*,*)'ps%nspec_alla=',ps%nspec_alla
         write(*,*)'ps%nspec_tha=',ps%nspec_tha


!!$         ll=1
!!$         iout=5
!!$c         write(iout,'(A,1pe12.5,i3)') '  density      Temperature'//
!!$         write(*,'(A,1pe12.5,i3)') '  density      Temperature'//
!!$     +        '   atm_no   charge   Species @ r_nc(ll),ll= ',r_nc(ll),ll
!!$
!!$!BH111108         do ii=1,ps%nspec_alla
!!$         do ii=1,nbulk
!!$!             Write(iout,'(4(1x,1pe12.5),3x,a)') densi(1,ii+1),
!!$             Write(*,'(4(1x,1pe12.5),3x,a)') en_nc(ll,ii),
!!$     +       temp_nc(ll,ii), dmas_nc(ii), charge_nc(ii)
!!$         enddo  ! On ii
!!$
!!$      endif ! On nbulk.gt.1
!!$
!!$      do k=1,nbulk
!!$         charge(k)=charge_nc(k)
!!$         dmas(k)=dmas_nc(k)
!!$      enddo
!!$
!!$
!!$      do k=1,nbulk
!!$         do j=1,ndens
!!$            dens1(j,k)=en_nc(j,k)
!!$            temp1(j,k)=temp_nc(j,k)
!!$            tpop1(j,k)=1.d0
!!$            vflow1(j,k)=0.d0
!!$         enddo
!!$      enddo
!!$      do j=1,ndens
!!$         zeff1(j)=zeff_nc(j)
!!$      enddo

         do i=1,ps%nspec_alla
            enein(:,i)=densi(:,i)
         enddo
         tein(:)=tempe(:)
         tiin(:)=tempi(:,1)
c        There is no input for FI temperatures into cql3d.
c        This is not difficult to change in cql3d, if/when needed.
c        (For the general cql3d FP species, bulk and FI are 
c        a single species).

      endif  !On cql3d_mode .eq. 'el+ions'


c---------------------------------------------------------------
c     Add in coding for 'ions+el' here.
c     Rearrange species for above 'el+ions' case.
c     May be sufficient for 'ions-aor' or 'ions-tor'
c----------------------------------------------------------------

      if (cql3d_mode .eq. 'ions+el') then
         write(*,*)
         write(*,*)'STOP: Need to add coding for this case'
         write(*,*)
         stop
      endif

c---------------------------------------------------------------
c     write all namelist to  cqlinput file
c----------------------------------------------------------------

      write(*,*)'before write_all_namelists'

      open(unit=64, file='cqlinput_new', delim = 'apostrophe',
     1         status = 'unknown', form = 'formatted')

      write(64, nml = setup0)
      write(64, nml = setup)
      write(64, nml = trsetup)
      write(64, nml = sousetup)
      write(64, nml = eqsetup)
      write(64, nml = rfsetup)
      write(64, nml = frsetup)

      close(64)
      write(*,*)'after write all namelists'

c     if have modified PS, then store it.
      if (PS_add_nml.ne.'disabled') then
         write(*,*)
     +        'prepare_cql3d_input: --storing cqlinput data '//
     +        'in current PS'
!        Stores the modified plasma state--doesn't commit it. 
!        The previous one is still around
         
         call ps_store_plasma_state(ierr)
         if(ierr .ne. 0) then
            write(*,*)
     +           'Cannot ps_store_plasma_state in prepare_cql3d_input'
         endif
         
      endif                     ! on PS_add_nml

      write(*,*)'prepare_genray_input: Normal end'
      write(*,*)

 999  continue

      end program prepare_cql3d_input


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



      subroutine transcribe(genraynml)
c
cBH131016:  Would be less confusing if genraynml ==> cql3d_nml
cBH131016L  everywhere in this subroutine.
c
c     Transcribes file genraynml to ./genray.in, except do nothing
c     if genray.in already exits.
c     Transcribes file genraynml to ./cqlinput, except do nothing
c     if cqlinput already exits.
c     This 'copy' is in lieu of access to operating system commands.
c
      character(len=*), intent(in) :: genraynml
      parameter(long_enough=100000)
      character(len=long_enough) :: line1          !Automatic array
      logical logic1

      write(*,*)
! ptb      write(*,*)  'transcribe genraynml to ./genray.in'
      write(*,*)  'transcribe genraynml to ./cqlinput'
      write(*,*)

      max_length=0

      write(*,*)'transcribe: genraynml =',genraynml
      inquire(file=genraynml,iostat=kiostat,opened=logic1,
     +     number=inumber)
      write(*,*)'transcribe: inquire on ',genraynml,', opened=',logic1,
     1   'iostat=',kiostat,'unit=',inumber
      
! ptb      if (trim(genraynml).eq.'genray.in' .or.
! ptb     +    trim(genraynml).eq.'./genray.in') then
! ptb         write(*,*)'transcribe: genray.in already exists'
! ptb         return
! ptb      endif

      if (trim(genraynml).eq.'cqlinput' .or.
     +    trim(genraynml).eq.'./cqlinput') then
         write(*,*)'transcribe: cqlinput already exists'
         return
      endif

c     Open genraynml and transcribe it
cBH131030      open(unit=20,file=genraynml,delim='apostrophe',
cHB131030: Use 'none' here, for formated read/write
      open(unit=20,file=genraynml,delim='none',
     +     status="old",iostat=kiostat)
      if (kiostat.ne.0) write(*,*)'transcribe: kiostat=',kiostat
      if (kiostat.ne.0) STOP 'transcribe: prblm with genraynml'

! ptb      open(unit=21,file='genray.in',status='replace',iostat=kiostat)
! ptb      if (kiostat.ne.0) write(*,*)'transcribe: kiostat=',kiostat
! ptb      if (kiostat.ne.0) STOP 'transcribe: prblm with opening genray.in'

!PTB 11-14-2013      open(unit=21,file='cqlinput',delim='apostrophe',
      open(unit=21,file='cqlinput',delim='none',
     +     status='replace',iostat=kiostat)
      if (kiostat.ne.0) write(*,*)'transcribe: kiostat=',kiostat
      if (kiostat.ne.0) STOP 'transcribe: prblm with opening cqlinput'

 3    read(unit=20,fmt='(a)',end=4) line1
      len_line1=len_trim(line1)
      if (len_line1.gt.max_length) max_length=len_line1
      write(unit=21,fmt=*) trim(line1)
      if (len_line1.ge.(long_enough-1)) then
         write(*,*)'ain_transcribe,long_enough',len_line1,long_enough
         STOP 'Adjust long_enough'
      endif
      go to 3
 4    continue

      write(*,*)'ain_transcribe:  max line length =',max_length
      close(20)
      close(21)
 
      return
      end subroutine transcribe



          
      subroutine zone_check(nrho, x_out, x_in)
    
!     unravel zone points to boundary points
!     nrho-1 = dimension of x_in bin-centered data
!     nrho   = dimension of x_out bin-boundary data

      use swim_global_data_mod, only :
     1 rspec, ispec,     ! int: kind specification for real and integer
!            & swim_string_length,  ! length of strings for names, files, etc.
     1 swim_error         ! error routine


      integer, intent(in) :: nrho    ! # of pts covering [0:1]
      REAL(KIND=rspec), intent(in) :: x_in(nrho-1)  ! zone values
      REAL(KIND=rspec), intent(out) :: x_out(nrho)  ! boundary values
      
      integer :: irho
      
      x_out(nrho) = 0.5 * (3.0 * x_in(nrho -1) - x_in(nrho-2))
      do irho = nrho - 1, 1, -1
         x_out(irho) = 2.0 * x_in(irho) - x_out(irho + 1)
      end do

      end subroutine zone_check

      
