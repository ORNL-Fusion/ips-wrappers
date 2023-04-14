!toric preprocessor heavily based on prepare_aorsa_input.f90
!JCW 2007, 2008
!Last mod 25 Jan 2008 - added optional command line arg for state file

!--------------------------------------------------------------------------
!
! Multi-toroidal mode version 11/18/2012  (Batchelor)
!
!	This version allows for concurrently running multiple TORIC instances 
!	with multiple toroidal modes and possibly for multiple power sources. The
!	The source/mode data comes in from the plasma state.  The the user
!	should put the source and mode data into plasma state machine description and
!	shot configuration namelist files, which are read by do_toric_init.f90 and
!	put into the initial plasma state during the init phase. The multi-mode
!	python component rf_ic_toric_mm.py generates separate work sub-directoies
!	for each instance.  This code generates a separate torica.inp file for each source/
!	mode and writes it to the appropriate subdirectory.  It also generates other files
!	that are identical for all instances.  The python component puts links to these in
!	each sub-directory.
!
!	1) ICRF relevant machine description data now comes in through the plasma state 
!		 The ICRF relevant items in the machine description are:
!        icrf_src_name  = number & name of ICRF sources
!        ant_model = antenna model filenames (1 per antenna source)
!        nrz_antgeo = number of (R,Z) points, antenna geometries
!        max_nrz_antgeo = Maximum size of variable length enumeration: NRZ_ANTGEO
!        R_antgeo = antenna geo: R pts 
!        Z_antgeo = antenna geo: Z pts 
!        dx_fshield = distance, antenna to Faraday shield
!
!	2) ICRF relevant shot configuration now comes in throught he plasma state. 
!		 The ICRF relevant items in the shot configuration section are:
!        RFMIN = ICRF minority species
!        freq_ic = frequency on each ICRF source
!        n_straps = number of straps in the antenna
!        max_n_straps = Maximum size of variable length enumeration: N_STRAPS
!        R_ant = major radius of antenna
!        Z_ant = height of antenna
!        Z_mid_ant = center of antenna relative to TF coil midplane
!        num_nphi_vac = number of non-zero n_phi values
!        max_num_nphi_vac = Maximum size of variable length enumeration: NUM_NPHI_VAC
!        Real_ant_coef = real part of Fourier Coef
!        Imag_ant_coef = imag part of Fourier Coef
!        num_nphi = number of non-zero n_phi values
!        max_num_nphi = Maximum size of variable length enumeration: NUM_NPHI
!        nphi = n_phi wave spectrum from antenna
!        wt_nphi = vacuum spectrum n_phi weight
!
!	3)  Very important note: The TORIC variables FREQCY and NPHI can also appear in the
!		machine.inp file.  For multi-mode operation this doesn't make sense because they
!		are now arrays.  So in multi-mode use these should be removed from the machine.inp
!		file.  However, for backward compatibility for use with single modes, I am
!		testing for non-zero values (i.e. non-default values) coming in from machine.inp.
!		If present they will over-ride the values in plasma state and go into the 
!		torica.inp file.  This is bad practice though.  The values should be put into the
!		mdescr.nml and sconfig.nml files that are read by do_toric_init_mm.f90.
!
!	3)	The code takes one command-line argument, the name of the current plasma state
!		to be initialized.  The names of the namelist files are assumed to be
!		mdescr_namelist.dat and sconfig_namelist.dat.
!
!--------------------------------------------------------------------------


      program prepare_toric_input_mm
      USE plasma_state_mod
!--------------------------------------------------------------------------

    use swim_global_data_mod, only : &
            & rspec, ispec, &                ! int: kind specification for real and integer
            & swim_error                     ! error routine
!            & swim_string_length, &         ! length of strings for names, files, etc.

      implicit none

      integer :: ierr, i, isp, iwarn, irho ! PTB
      logical :: lex
      integer, parameter :: swim_string_length = 256  !for compatibility LAB

      character(len =swim_string_length) :: cur_state_file, cur_geq_file, &
      				subdirectory_path, path


      real(rspec), dimension(:),allocatable :: tmp_prof
      real(rspec), dimension(:),allocatable :: vol_int
! PTB begins -
      real(rspec), dimension(:,:),allocatable :: ns_tha
      real(rspec), dimension(:,:),allocatable :: ts_tha
      real(rspec), dimension(:,:),allocatable :: v_pars
      real(rspec), dimension(:),allocatable :: aa_prof
      real(rspec), dimension(:),allocatable :: bb_prof
      real(rspec), dimension(:),allocatable :: x_orig
      real(rspec), dimension(:),allocatable :: x_toric
      real(rspec) :: tol_zero = 1.0E-12_rspec, dVol, dVol_int, Q_ps, Q_int

! PTB ends



!------Namelist inputs-------------
! see the rf component write-up in the svn repository for more description of
! the inputs

!  Mode of TORIC utilization "equil" to preprocess equilibrium, "toric" to run solver
      character(10):: toricmode='toric'

! Dimensions of the problem
      integer :: nvrb=3       ! Generally three vector components
! Poloidal resolution
      integer :: ntt = 64
      integer :: nmod =31     ! nmod will be set as nmod=ntt/2-1, 
                              ! unless user enters a non-zero positive value
                              ! less than ntt/2-1 to conserve memory at 
                              ! the cost of efficiency.
! Radial resolution: Number of radial elements in the plasma
      integer :: nelm = 240
!N.B. max memory usage for double precision is approx 16*3*(6*Nmod)^2*Nelm/Nprocs
! on a machine with 2GB of ram, the practical in core limit in Nm=255, Nelm=240, but
! this run will take several hours, so most likely, memory will not be a problem since
! many processors will be used to make run times reasonable

!Don't change these three unless you really know what you're doing
! Number of points in vacuum
      integer :: nptvac = -1
! Number of poloidal modes in vacuum
      integer :: mxmvac = 15
! Number of gaussian quadrature points in radial finite elements
      integer :: ngaus = 3

! Namelist inputs for choice of the model, rarely need to be changed
      integer ::   isol=1     ! to solve or not to solve
      integer ::   mastch=2   ! solver choice (=2 uses the new 3D parallel solver)
      integer ::   iwdisk = 1 ! output results for plotting in file
      integer ::   iflr = 1   ! flr effects
      integer ::   iqlmin =0  ! minority tail profile read or not
      integer ::   iqtor = 1  ! toroidal broadening in Z functions
      integer ::   ibweld=1   ! eld damping of ibw's on
      integer ::   icosig=0   ! extra damping for IBW when (kperp rho_i)^2>2
      integer ::   ibpol=1    ! poloidal field on
      integer ::   imdedg=2   ! JPW pollution scheme for vacuum layer modelling
      integer ::   iezvac=1   ! Ez suppressed in vacuum
      integer ::   icoll=0    ! No collisions
      integer ::   iclres=0   ! adds collisions around isolated ion-ion 
                              ! resonances. Use with care!
      integer ::   iregax=1   ! turn on regularization at the magnetic axis 
                              ! (changes equil near axis a bit)
      integer ::   bscale=12  ! blocksize scaling factor for parallel runs only
      integer ::   pcblock=4  ! number of processors used for the sub-blocks
      logical ::   use_incore=.true. 
                              ! default is to use in-core memory
                              ! but be careful to have enough processors

! Namelist inputs for control of the output, most are off to avoid too much data dumped
      integer ::   iout=0, ipltht=0, idlout=1
! NetCDF output for plots on by default.  torica.sol file has solutions but not necessarily
! many of the parameters
      integer ::   io_ncdf = 1

! Wave and antenna parameters (default values for now, later 
!   use ps%ant_model file for machine state)
      integer     :: nphi=0
      real(rspec) :: freqcy=0.0_rspec

!units are in cm for lengths
      real(rspec) :: antlen=24.0_rspec, antlc=1.0_rspec, &
         theant=0.0_rspec

! Parameters for Ehst-Karney current drive estimate and collisions
      real(rspec) ::  zeff=1._rspec, enhcol=1
      real(rspec) ::  dnures=1.0_rspec !width of layer for iclres (cm)
      real(rspec) ::  tnures=1.0_rspec !strength of layer for iclres (A.U.)
! profiling
      logical :: timing_on=.true.
!
! MHD equilibrium and profiles data
!
      integer :: idprof=1, igsmhd=1 !use equil_file
      character(24) :: equil_file='equigs.data', profnt_file='equidt.data', &
                       qlmin_file='', gfile
! Path to the file containing the namelists
      character(24) :: inputpath = './'
! Path to the buffer file for the solver using virtual memory
      character(50) :: scratchpath = '/scratch/scratchdirs/u1565/tmp'

! Namelist for numerical magnetic equilibrium
! PTB - changed the deafult fitting for the MHD equilibrium form Chebyschev 
! polynomials to splines as the they seem to work better, especially near
! the separatrix.
!     integer :: nmhd, intchb=7
      integer :: nmhd, intchb=0
! Namelist inputs for positioning the plasma with respect to the wall
      real(rspec), dimension(:) :: dist_plawall(4)=        &
     &                         (/10._rspec,10._rspec,10._rspec,10._rspec/)
      real(rspec) ::  dist_plafars=0._rspec, dist_plaant=0.5_rspec
      real(rspec) ::  so_thickness=0.5_rspec

! Namelist inputs for profile specification
!  Density and temperature profiles - namelist and file inputs
!
      integer, parameter :: nspmx = 16
! NOTE: A maximum of 15 ion species allowed
! The place,nspec+1, is reserved for the electrons in mod_direl.
! mainsp is used to impose charge neutrality
      integer ::   nspec,  mainsp, iudsym=0

      integer, dimension(:) ::                                       &
     &             atm(nspmx)=0._rspec,       azi(nspmx)=0._rspec

      real(rspec), dimension(:) ::                                   &
     &             aconc(nspmx)=0._rspec,     tempic(nspmx)=0._rspec,& 
     &             tisepr(nspmx)=0._rspec,    glti(nspmx)=0._rspec,  &
     &             pptii(nspmx)=0._rspec,     pptie(nspmx)=0._rspec

! Namelist controls for the RF minority population

      real(rspec) :: q_rfmin = 1.0_rspec, qatom_rfmin = 1.0_rspec, &
     &               m_rfmin = 1.0_rspec, fracmin = 0.04_rspec
      integer :: isThermal = 1
      character(24) :: rfmin_name = 'H_min', kdens_rfmin = 'fraction'

      integer :: i2mex_remap = 1, i2mex_mx = 2, i2mex_npsi = 0, &
     &  i2mex_kb = 0, i2mex_direct = 1, inputFormat = 3, &
     &  ic1, inc, imom1, i2mex_last_norm_surface_is, ntheta
          ! if direct representation BR, BZ, BPhi
          ! = BR, BZ, BPhi(R, Z) should be allowed
          ! EFIT/geqdsk input file format


!TORIC namelist blocks (some variables are not written for simplification)
!Most of these variables are described in man_toric and 
!initalized in t4_mod_public.F

!originally in t4_aamain.F
!specifies task for toric
      namelist /toric_mode/ toricmode

!originally in t4_torica.F
!specifies general wave parameters, some numerical parameters
      namelist /toricainp/ &
     &   nvrb,   nmod,   ntt,    nelm,   nptvac, mxmvac, &
     &   freqcy, nphi,   antlen, antlc,  theant, &
     &   iflr,   ibpol,  iqtor,  icoll,  enhcol, &
     &   imdedg, iezvac, ibweld, icosig, iregax, &
     &   isol,   mastch,         iout,   idlout, io_ncdf, &
     &   iwdisk, ipltht, zeff,   iclres, dnures, tnures, &
     &   timing_on, scratchpath, bscale, use_incore, pcblock

!originally in t4_mod_toi2mex.F
!specifies numerical equilibrium (EFIT usually) settings
      namelist /gfreadinp/ &
     &   gfile, inputFormat, nmhd, ntheta, &
     &   i2mex_remap, i2mex_mx, i2mex_npsi, i2mex_kb, i2mex_direct, &
     &   i2mex_last_norm_surface_is, ic1, inc, imom1

!originally in t4_mod_equil.F
!specifies plasma equilibrium settings and profiles
      namelist/equidata/ igsmhd, intchb, idprof, &
     &   nspec, iudsym, mainsp, atm,  azi, so_thickness,  &
     &   dist_plafars,   dist_plaant,    dist_plawall,  &
     &   inputpath,      equil_file,     profnt_file

! This namelist group was added to specify parameters controlling nonthermal 
! ion populations such as the ICRF minority in the machine.inp file.
! This namelist group is NOT written to the torica.inp file. It only
! serves to transfer information between the machine.inp file and the
! Fortran wrappers do_toric_init and prepare_toric_input. 

      namelist /nonthermals/ &
     &   fracmin, q_rfmin, qatom_rfmin, m_rfmin, rfmin_name, &
     &   kdens_rfmin, isThermal 

!other variables for output of profiles
    integer :: nprodt, nproeq, kdiff_idens=1, kdiff_itemp=1
    
!Not used yet
    real(rspec) :: prfin
!I/O units
    integer :: inp_unit, out_unit, iarg
!Other magic numbers
    real(rspec), parameter :: MHz=1.e6_rspec, cubic_cm=1.e-6_rspec
    integer, parameter :: isrc=1 !isrc == 1..nicrf_src in plasma state
!
!   PLASMA STATE DATA that will be given to toric via a namelist
!
!--------------------------------------------------------------------------
!JCW we leave this here for future use, but no S_ variables are being used right now
!6 Mar 2007   
    !-----------------------------------
    ! Time at beginning and end of time step
    !-----------------------------------
    REAL (KIND = rspec) ::     &
      & S_t0,                  &   ! time at beginning of step [msec]
      & S_t1                       ! time at end of step [msec]
   
    !-----------------------------------
    ! Basic Geometry
    !-----------------------------------
    REAL (KIND = rspec) ::     &
        S_r_axis,              & ! major radius of magnetic axis [m]
        S_z_axis,              & ! Z of magnetic axis [m]
        S_r0_mach,             & ! Z of machine center [m]
        S_z0_mach,             & ! major radius of machine center [m]
        S_r_min,               & ! major radius of inside of bounding box [m]
        S_r_max,               & ! major radius of outside of bounding box [m]
        S_z_min,               & ! Z of bottom of bounding box [m]
        S_z_max                  ! Z of top of bounding box [m]
            
    !-----------------------------------
    ! Particle Species
    !-----------------------------------
    
    INTEGER :: S_nspec       ! number of ion species = nspec_th + nspec_nonMax

    !-----------------------------------
    ! Main (thermal) Plasma Species
    !-----------------------------------m
    
    integer, parameter :: nrho_max = 220
!   integer, parameter :: n_spec_th_max = 5
    integer, parameter :: n_spec_th_max = 30
!   integer, parameter :: n_spec_max = 7
    integer, parameter :: n_spec_max = 32
    integer, parameter :: n_spec_nm_max = 2

    
    
    INTEGER :: S_nspec_th                      ! number of thermal ion species
    character(len = 32) ::  &
        S_s_name(0:n_spec_max)                 ! names of main species, (0:nspec_th)
    REAL (KIND = rspec) :: &
        S_q_s(0:n_spec_th_max),              & ! charge of species s [C], (0:nspec_th)
     &  S_m_s(0:n_spec_th_max)                 ! mass of species s [kg], (0:nspec_th)
    
    INTEGER :: S_nrho_n          ! number of rho values in thermal species density grid
    REAL (KIND = rspec) :: &
        S_rho_n_grid(nrho_max),            & ! rho values in density grid, (1:nrho_n)
     &  S_n_s(nrho_max, 0:n_spec_th_max),  & ! density profile of species s, (1:nrho_n, 0:nspec_th)
     &  S_q_impurity(nrho_max),            & ! effective impurity charge profile, (1:nrho_n)
     &  S_m_impurity(nrho_max)               ! effective impurity mass profile, (1:nrho_n)
 
    INTEGER :: S_nrho_T                  ! number of rho values in temperature grid
    REAL (KIND = rspec) :: &
        S_rho_T_grid(nrho_max),       &  ! rho values in temperature grid, (1:nrho_T)
      & S_T_s(nrho_max, 0:n_spec_th_max) ! Temperature profile of species s, (1:nrho_T, 0:nspec_th)
 
    INTEGER :: S_nrho_v_par      ! number of main rho values in parallel velocity grid
!    REAL (KIND = rspec), ALLOCATABLE :: &
!        PS_rho_v_par_grid(:),   & ! rho values in parallel velocity grid, (1:nrho_v_par)
!      & PS_v_par_s(:, :)        & ! v parallel profile of species s, 
                                   ! (1:nrho_v_par, 0:nspec_th)
 
    !-----------------------------------
    ! Non-Maxwellian Species
    !-----------------------------------
    
    INTEGER :: S_nspec_nonMax    ! number of non-Maxwellian species
    character(len=32), dimension(n_spec_nm_max ) :: &
        S_nonMax_name         ! names of non-Maxwellian species, (1:nspec_nonMax)
    
    REAL (KIND = rspec), dimension(n_spec_nm_max ) :: &
        S_q_nonMax_s,       & ! charge of species s [C], (1:nspec_nonMax)
        S_m_nonMaX_s          ! mass of species s [kg], (1:nspec_nonMax)
    
    INTEGER :: S_ntheta_n        ! number of theta values in 2D density grid

    REAL (KIND = rspec), ALLOCATABLE :: &
        S_n_nonMax2D_s(:, :,:)  ! 2D density profile of non-Maxwellian species s,
                                 ! (1:nrho_n, 1:ntheta_n, 1:nspec_nonMax)

    REAL (KIND = rspec), ALLOCATABLE :: &
        S_n_nonMax_s(:, :)   ! Flux surface average density profile of 
                              ! non-Maxwellian species s, (1:nrho_n, 1:nspec_nonMax)
 
    character(len = swim_string_length) :: &
        S_dist_fun_s         ! distribution function of non-Maxwellian  
                              ! species s, (1:nspec_nonMax) N.B. For now a distribution
                              ! function type is a file name


    !-----------------------------------
    ! Magnetics
    !
    ! magnetics: B(x), magnetic field.  Like distribution_fn there is a
    ! user-defined type that contains a file with the data.
    ! AORSA and TORIC get all their magntics data by reading an eqdisk file.
    ! Eventually all the magnetics data will appear separately in the Plasma
    ! state.
    !-----------------------------------

    character(len = swim_string_length) :: S_eqdsk_file   ! eqdisk file
        
    REAL (KIND = rspec) ::  &
        S_B_axis               ! Field at magnetic axis [T]

    !--------------------------------------------------------------------------
    !
    ! RF Data
    !
    ! Allow multiple RF sources, ICF, ICRH. So RF_frequency, power, etc may
    ! not be a scalar in that case.
    !   Assumption 1: the RF component invocation will
    !       involve a loop over each RF source, each of which can have its own
    !       RF_frequency, etc.  (vs. adding another dimension to the arrays).
    !   Assumption 2 : each source involves invoking another executable.
    !--------------------------------------------------------------------------
   
    
    INTEGER :: S_nrf_src         ! number of RF sources
       ! names of rf sources, (1:nrf_src)
        
 
    character(len = swim_string_length) :: S_ant_model_src !file name for antenna model 
    !---------------------------------------------------------------------------
    ! Note:
    ! Antenna model is currently defined in a file. The PREPARE_CODE_INPUT program
    ! should extract from it the data needed to define the geometry and operation
    ! i.e. phasing or mode number spectrum
    !
    ! For these, see the example namelist attached to the end
    ! of the aorsa.doc file sent by Fred to Swim list.  For now the data in this
    ! file includes:
    !   nphi =toroidal mode number (find another name not conflicting toroidal angle)
    !   antlen = vertical height of antenna [m]
    !   dpsiant0 = radial thickness of antenna in rho
    !   rant = radial location of antenna in major radius [m]
    !   yant = vertical location of antenna center [m]
    !
    !  N.B. In this scheme the toroidal mode number comes in through the antenna model
    !  The antenna geometry should eventually come from one of the standard machine
    !  definition files. 
    !
    !---------------------------------------------------------------------------
  
  
    ! RF Outputs that go back into the Plasma State.  Profiles are flux surface averages.
    
    ! N.B. We will want to put in 2D power deposition profiles, but I don't think they
    ! are needed for our initial coupling
        
    INTEGER ::  &
        S_nrho_prf,    &   ! number rho values for RF power deposition grid
        S_ntheta_prf       ! number of theta values in 2D RF power dep grid
        
    REAL (KIND = rspec), ALLOCATABLE :: &
        S_rho_prf_grid(:), &        ! rho values in RF power deposition grid, (1:nrho__prf)
        S_prf2D_src_s(:,:,:,:),   & ! 2D Power deposition from each source into each 
                                        ! species, (1:nrho__prf, 1:nrf_src, 0:nspec)
        S_prf_src_s(:,:,:)    ! Power deposition profile from each source into each 
                              ! species, (1:nrho__prf, 1:nrf_src, 0:nspec)

! Total rf power deposition profile into each species summed over sources, (1:nrho__prf, 0:nspec)
    real(kind = rspec) :: S_prf_total_s(nrho_max,0:n_spec_max)

! # of rho values for RF current drive grid for each species
    integer :: S_nrho_cdrf(n_spec_max)
        
    REAL (KIND = rspec), ALLOCATABLE :: &
        S_rho_cdrf_grid(:),    & ! rho values in RF current drive grid, (1:nrho__cdrf)
        S_cdrf_src_s(:,:,:),   & ! Driven current profile from each source, in each species
                                 ! (1:nrho__cdrf, 1:nrf_src, 1:nspec_nonMax)
        S_cdrf_total_s(:,:)      ! Total current driven by all sources in each species
    
 
    character(len = swim_string_length), dimension(n_spec_nm_max)  :: &
        S_ql_operator          ! quasilinear operator for each non Maxwellian--file name 
                               ! species, (1:nspec_nonMax)
    character(len = swim_string_length), dimension(n_spec_nm_max) :: &
        S_distribution_fun     ! distribution function for each non Maxwellian species


    namelist /state/ S_t0, S_t1, S_r_axis, S_z_axis, S_r0_mach, &
       S_z0_mach, S_r_min, S_r_max, S_z_min, S_z_max, &
       S_nspec,S_nspec_th,S_s_name, S_q_s, S_m_s,  &
       S_nrho_n, S_rho_n_grid, S_n_s, S_q_impurity, S_m_impurity, &
       S_nrho_T, S_rho_T_grid , S_T_S !, S_ant_model_src, S_ql_operator, &
       !S_distribution_fun     


! END  PLASMA STATE DATA that will be given to toric via a namelist 
! BEGIN executable section
!--------------------------------------------------------------------------------------

      write(*,*) 'Prepare toric input  mm'
      call get_arg_count(iarg)
      SELECT CASE (iarg)
      case(0)
         cur_state_file="cur_state.cdf"

      case(1)
         call get_arg(1,cur_state_file)

      case(2:)
         write(0,*) 'Error. Illegal number of arguments.'
         write(0,*) 'prepare toric usage: '
	 write(0,*) 'prepare_toric_input cur_state_file'
	 stop 'incorrect command line arguments'

      end select

      print*, 'using ', trim(cur_state_file), ' as default state.'

      call getlun(inp_unit,ierr)  ;  call getlun(out_unit,ierr)

	!------------------------------------------------------------------------------------
	!  Get current plasma state 
	!------------------------------------------------------------------------------------
			
    call ps_get_plasma_state(ierr, trim(cur_state_file))
    if (ierr .ne. 0) then
       print*, 'model_RF_IC:failed to get_plasma_state'
       call exit(1)
    end if

!frequency presently not set in plasma state
      freqcy = ps%freq_ic(isrc)*MHz !Hz to MHz
      prfin =  ps%power_ic(isrc) !watts
!     nspec = ps%nspec_th !+1 !toric includes electrons in species count

! PTB - begins
! Use the abridged species list index count (nspec_alla) that includes:
! Electrons,1 fully stripped light ion, 2 model impurity ions, and up to 3 non-thermal ion species 
! (energetic NBI + energetic minority ion + fast fusion alpha)
! TORIC expects the ion species list to include all thermal ions plus the non-thermal species
! PTB - ends
!

! Read in partial namelist with override values such as resolution, or mode      

! Flexibility for case specific settings, only include values which you want to 
! override the defaults. 

! Do not set NSPEC in the machine.inp file if you want to use nspec = ps%nspec_alla + 1 because
! this will override the setting of NSPEC above !!!
      
      open(unit=inp_unit, file='machine.inp', status='old', &
              form='formatted')
      INQUIRE(inp_unit, exist=lex)
      IF (lex) THEN
         read(inp_unit, nml = toricainp)
         read(inp_unit, nml = equidata)
         read(inp_unit, nml = nonthermals)
      ELSE
         write(*,*) &
            'machine.inp does not exist or there was a read error'
      END IF
      close(inp_unit)
      IF (nspec > nspmx ) THEN
         write(*,*) "Error, nspec > nspmx in toric, reducing to nspmx"
         nspec=nspmx
      END IF

!radial profiles generation, these are output to equilequ_file
      s_nrho_n = ps%nrho
      s_nrho_t = ps%nrho
      s_rho_n_grid(1:s_nrho_n) = ps%rho  !grids are the same for both n and T  FIX
      s_rho_t_grid(1:s_nrho_t) = ps%rho

      ! now the species
      ! note that in the plasma state, index 0 is electrons.
      ! Toric has integer normalized units and charges
! PTB - begins
! read ion species charge states and masses from the plasma State 
! PTB atm(1:nspec) = NINT(ps%m_s(1:nspec)/ps_mp)
! PTB azi(1:nspec) = NINT(ps%q_s(1:nspec)/ps_xe)
      zeff = ps%zeff(1)
      atm(1:ps%nspec_tha) = NINT(ps%m_ALLA(1:ps%nspec_tha)/ps_mp)
      azi(1:ps%nspec_tha) = NINT(ps%q_ALLA(1:ps%nspec_tha)/ps_xe)

      isp = ps%nspec_tha 
      if (allocated(ps%nbeami)) then
         isp = isp + 1
         atm(isp) = NINT(ps%m_ALLA(ps%snbi_to_alla(1))/ps_mp)
         azi(isp) = NINT(ps%q_ALLA(ps%snbi_to_alla(1))/ps_xe)
      endif

      if (allocated(ps%nmini)) then
         isp = isp + 1
         atm(isp) = NINT(ps%m_ALLA(ps%rfmin_to_alla(1))/ps_mp)
         azi(isp) = NINT(ps%q_ALLA(ps%rfmin_to_alla(1))/ps_xe)
      endif

      if (allocated(ps%nfusi)) then
         isp = isp + 1
         atm(isp) = NINT(ps%m_ALLA(ps%sfus_to_alla(1))/ps_mp)
         azi(isp) = NINT(ps%q_ALLA(ps%sfus_to_alla(1))/ps_xe)
      endif

      nspec = isp 
      mainsp = 1
! PTB - ends

!toric doesn't use these right now, but output them in state namelist for
!possible future use
      s_s_name(0:ps%nspec_th) = ps%s_name(0:ps%nspec_th)
      s_m_s(0:ps%nspec_th)    = ps%m_s(0:ps%nspec_th)
      s_q_s(0:ps%nspec_th)    = ps%q_s(0:ps%nspec_th)

      idprof = 1  !use numerical profiles
      gfile = trim(ps%eqdsk_file)

      toricmode='toric'
      open(unit=out_unit, file='torica.inp',                &
        status = 'unknown', form = 'formatted',delim='quote')

	!------------------------------------------------------------------------------------
	!  Loop over nicrf_src and num_nphi to write multiple torica.inp files
	!------------------------------------------------------------------------------------
      
      DO i_source = 1, ps%nicrf_src
		  Do i_nphi = 1, num_nphi(i_source)
		  
			  freqcy = ps%freq_icrf(i_source)
			  nphi = ps%nphi(i_nphi, i_source)
			  
			  path = trim(subdirectory_path)||'torica.inp'
			  open(unit=out_unit, file=path,                &
				status = 'unknown', form = 'formatted',delim='quote')
			  write(out_unit, nml = toric_mode)
			  write(out_unit, nml = toricainp)
			  write(out_unit, nml = equidata)
		      close(out_unit)

		  END DO      
	  END DO

      nprodt = size(ps%rho)
      nproeq = size(ps%rho_eq)

! PTB begins -
      allocate( tmp_prof(nprodt))
      allocate( vol_int(nprodt))
      allocate( ns_tha(nprodt-1,ps%nspec_tha))
      allocate( ts_tha(nprodt-1,ps%nspec_tha))
      allocate( v_pars(nprodt-1,ps%nspec_tha))
      allocate( aa_prof(nprodt))
      allocate( bb_prof(nprodt))
      allocate( x_orig (nprodt-1))
      allocate( x_toric (nprodt))
! PTB ends 

      open(unit=out_unit, file=profnt_file,              &
         status = 'unknown', form = 'formatted')
      write(out_unit,'(A10,5i4)') 'PS-state', nprodt, nspec, mainsp, &
                                  kdiff_idens, kdiff_itemp
      do  isp=1,nspec
         write(out_unit,'(2i4)')  atm(isp), azi(isp)
!         write(*,*)  atm(isp), azi(isp)
      end do
! PTB begins
!
! Define mid-cell values from the orginal Plasma State radial grid - "ps%rho" and use this array
! for all interpolations
!
     do irho = 1,nprodt-1
        x_orig(irho) = 0.5 * (ps%rho(irho) + ps%rho(irho+1))
      end do
! PTB end

! TORIC uses only one radial mesh for density and temperature profiles that is
! that is defined in terms of the sqrt (Psi_pol) - normalized
!
      x_toric = sqrt(ps%psipol / ps%psipol(nprodt))

!     write(out_unit,'(A10)')  'rho'
!     write(out_unit,'(5E16.9)')  ps%rho !check units

! PTB begins
! Need to actually write out the sqrt (Psi_pol) - normalized
      write(out_unit,'(A20)')  'Sqrt(Psi / Psilim)'
      write(out_unit,'(5E16.9)') x_toric
! PTB ends

      write(out_unit,'(A10)')  'n_e'
      write(*,*) 'ps%ns(1,0)', ps%ns(1,0)
!
! Interpolate the electron density profile from the Plasma State grid to the Toric grid
!
         call ps_user_1dintrp_vec(x_toric, x_orig, ps%ns(:,0), &
               tmp_prof(:),ierr )
         if(ierr .ne. 0) stop 'error interpolating PS electron density profile onto Toric grid'
!
! Interpolate the volume profile from the Plasma State grid to the Toric grid
!
         call ps_user_1dintrp_vec(x_toric, ps%rho_eq, ps%vol(:), &
               vol_int(:),ierr )
         if(ierr .ne. 0) stop 'error interpolating PS volume onto Toric grid'
!
! PTB Compare the volume average of the orginal density profile from the Plasma State 
! and the volume average of the interpolated profile
!
         Q_ps = 0.0_rspec
         Q_int = 0.0_rspec
         do irho = 1,nprodt-1
            dVol = ps%vol(irho+1) - ps%vol(irho)
            Q_ps = Q_ps + 0.5_rspec * (ps%ns(irho,0) + &
                   ps%ns(irho+1,0)) * dVol
         end do
         do irho = 1,nproeq-1
            dVol_int = vol_int(irho+1) - vol_int(irho)
            Q_int = Q_int +  0.5_rspec * (tmp_prof(irho) + &
                   tmp_prof(irho+1)) * dVol_int
         end do
         Q_ps = Q_ps / ps%vol(nprodt)
         Q_int = Q_int / ps%vol(nprodt)
      write(*,*) '<n_e(m-3)-PS> =',  Q_ps
      write(*,*) '<n_e(m-3)-Interpolated> =', Q_int
      write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3

      write(out_unit,'(A10)')  't_e'
         call ps_user_1dintrp_vec(x_toric, x_orig, ps%Ts(:,0), &
               tmp_prof(:),ierr )
         if(ierr .ne. 0) stop 'error interpolating PS electron temperature profile onto Toric grid'
      write(out_unit,'(5E16.9)')  tmp_prof   !keV

! PTB - begins
        call ps_tha_fetch (ierr, ps, tol_zero, &
             ns_tha, ts_tha, v_pars, iwarn)
         if(ierr .ne. 0) stop 'error getting the density and temperature data from the abridged species list'
      do isp=1,ps%nspec_tha
         write(out_unit,'(A4,I2.2)')  'n_i_',isp
         call ps_user_1dintrp_vec(x_toric, x_orig, ns_tha(:,isp), &
               tmp_prof(:),ierr )
         if(ierr .ne. 0) stop 'error interpolating PS ion density profile onto Toric grid'
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3

         write(out_unit,'(A4,I2.2)')  't_i_',isp
         write(*,*) "Thermal ion name, A, Z, dens, temp = "
         write(*,*) trim(ps%alla_name(isp)), atm(isp),azi(isp),ns_tha(1,isp),ts_tha(1,isp)
         call ps_user_1dintrp_vec(x_toric, x_orig, Ts_tha(:,isp), &
               tmp_prof(:),ierr )
         if(ierr .ne. 0) stop 'error interpolating PS ion temperature profile onto Toric grid'
         write(out_unit,'(5E16.9)')  tmp_prof !keV
      end do
! PTB - ends

! PTB - begins
!
! Move on to the non-thermal ion species - in order of NBI, RF minority, and fusion alphas
!
! (1) Write the neutral beam injection data to the Toric equidt.data file:
!
! First interpolate the NBI density [nbeami(:,1)] from the NuBeam grid onto the Toric input data 
! grid. Then write it to the equidt.data file.
!
      if(allocated(ps%nbeami)) then 
       isp = ps%snbi_to_alla(1)
         write(*,*) "Fast ion name, A, Z = ", trim(ps%alla_name(isp)), &
     &              NINT(ps%m_alla(isp)/ps_mp), NINT(ps%q_alla(isp)/ps_xe)
         call ps_user_1dintrp_vec(x_toric, ps%rho_nbi, ps%nbeami(:,1), &
               tmp_prof(:),ierr )
         if(ierr .ne. 0) stop 'error interpolating PS NBI density profile onto Toric grid'
         write(out_unit,'(A4,I2.2)')  'n_i_',isp
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3
!
! Next interpolate the NBI energies [eperp_beami(:,1) and epll_beami(:,1)] from the NuBeam grid 
! onto the Toric input data grid. Then compute the equivalent temperature profile for the NBI and
! write it to the equidt.data file.
!
          call ps_user_1dintrp_vec(x_toric, ps%rho_nbi, ps%eperp_beami(:,1), &
               aa_prof(:),ierr )
          if(ierr .ne. 0) stop 'error interpolating PS NBI perp. energy profile onto Toric grid'
          call ps_user_1dintrp_vec(x_toric, ps%rho_nbi, ps%epll_beami(:,1),  &
               bb_prof(:),ierr )
          if(ierr .ne. 0) stop 'error interpolating PS NBI parallel energy profile onto Toric grid'
          tmp_prof = 0.667 * (aa_prof + bb_prof)
          write(out_unit,'(A4,I2.2)')  't_i_',isp
          write(out_unit,'(5E16.9)')  tmp_prof !keV
      endif
!
! (2) Write the rf minority tail data to the Toric equidt.data file:
!
      if(allocated(ps%nmini)) then
       isp = ps%rfmin_to_alla(1)
         write(*,*) "Fast ion name, A, Z = ", trim(ps%alla_name(isp)), &
     &              NINT(ps%m_alla(isp)/ps_mp), NINT(ps%q_alla(isp)/ps_xe)
!
! First interpolate the electron density [ps%ns(:,0] from the PS grid onto the Toric input data 
! grid. Then write it to the equidt.data file.
!
! If (kdens_rfmin .EQ. 'fraction') then assume PS data is not available for nmini and instead
! compute nmini = fracmin * ne
!
        if (trim(kdens_rfmin) .EQ. 'fraction') then
         call ps_user_1dintrp_vec(x_toric, ps%rho, ps%ns(:,0), &
               tmp_prof(:),ierr )
         if(ierr .ne. 0) stop 'error interpolating PS electron density profile onto Toric grid'
         tmp_prof(:) = fracmin * tmp_prof(:)
         write(out_unit,'(A4,I2.2)')  'n_rfmin_',isp
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3
! 
! Next update ps%nmini in the Plasma State with th new minority ion density profile, by mapping
! fracmin * ps%ns(:,0) from the PS grid to the ICRF rho grid:
!
         call ps_user_1dintrp_vec(ps%rho_icrf,x_orig, fracmin*ps%ns(:,0), &
               ps%nmini(:,1),ierr )
         if(ierr .ne. 0) stop 'error interpolating new minority desnity profile onto PS grid'     
        endif
!
! If (kdens_rfmin .EQ. 'data') then assume nmini is available in the PS, read it, and interpolate
! it from the rho-icrf grid onto the Toric radial grid
!
        if (trim(kdens_rfmin) .EQ. 'data') then
         call ps_user_1dintrp_vec(x_toric, ps%rho_icrf, ps%nmini(:,1), &
               tmp_prof(:),ierr )
         if(ierr .ne. 0) stop 'error interpolating PS minority density profile onto Toric grid'
         write(out_unit,'(A4,I2.2)')  'n_rfmin_',isp
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3
        endif        
!
! If isThermal=1 then assume the RF minority tail temperature is equal to the bulk ion temperature.
! Set the minority tail temperature equal to the temperature profile of the first bulk ion species
! of the abridged species list, by interpolating that array from the PS grid onto the Toric input data
! grid, and then writing it to the equidt.data file.
!
         write(out_unit,'(A4,I2.2)')  't_rfmin_',isp
         if (ps%isThermal(1) .eq. 1) then
         call ps_user_1dintrp_vec(x_toric, ps%rho, Ts_tha(:,1), &
               tmp_prof(:),ierr )
         if(ierr .ne. 0) stop 'error interpolating PS RF minority temperature profile onto Toric grid'
         write(out_unit,'(5E16.9)')  tmp_prof !keV
         endif 
!
! If isThermal=2 then assume the RF minority tail temperature data is available in the PS. In this 
! case just interpolate the ICRF minority energies [eperp_mini(:,1) and epll_mini(:,1)] from the ICRF 
! grid (rho_icrf) onto the Toric input data grid. Then compute the equivalent temperature profile for 
! the RF minority profile. Finally write this array to the equidt.data file.
!
         if (ps%isThermal(1) .eq. 2) then
          call ps_user_1dintrp_vec(x_toric, ps%rho_icrf, ps%eperp_mini(:,1), &
               aa_prof(:),ierr )
          if(ierr .ne. 0) stop 'error interpolating PS RF minority perp. energy profile onto Toric grid'
          call ps_user_1dintrp_vec(x_toric, ps%rho_icrf, ps%epll_mini(:,1),  &
               bb_prof(:),ierr )
          if(ierr .ne. 0) stop 'error interpolating PS RF minority parallel energy profile onto Toric grid'
          tmp_prof = 0.667 * (aa_prof + bb_prof)
          write(out_unit,'(5E16.9)')  tmp_prof !keV
         endif         
      endif
!
! (3) Write the fast fusion alpha data to the Toric equidt.data file:
!
! First interpolate the alpha density [nfusi(:,1)] from the fusion grid onto the Toric input data 
! grid. Then write it to the equidt.data file.
!
      if(allocated(ps%nfusi)) then 
       isp = ps%sfus_to_alla(1)
         write(*,*) "Fast ion name, A, Z = ", trim(ps%alla_name(isp)), &
     &              NINT(ps%m_alla(isp)/ps_mp), NINT(ps%q_alla(isp)/ps_xe)
         call ps_user_1dintrp_vec(x_toric, ps%rho_fus, ps%nfusi(:,1), &
               tmp_prof(:),ierr )
         if(ierr .ne. 0) stop 'error interpolating PS fusion alpha density profile onto Toric grid'
         write(out_unit,'(A4,I2.2)')  'n_i_',isp
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3
!
! Next interpolate the alpha energies [eperp_fusi(:,1) and epll_fusi(:,1)] from the fusion grid 
! onto the Toric input data grid. Then compute the equivalent temperature profile for the NBI and
! write it to the equidt.data file.
!
          call ps_user_1dintrp_vec(x_toric, ps%rho_fus, ps%eperp_fusi(:,1), &
               aa_prof(:),ierr )
          if(ierr .ne. 0) stop 'error interpolating PS fusion alpha perp. energy profile onto Toric grid'
          call ps_user_1dintrp_vec(x_toric, ps%rho_fus, ps%epll_fusi(:,1),  &
               bb_prof(:),ierr )
          if(ierr .ne. 0) stop 'error interpolating PS fusion alpha parallel energy profile onto Toric grid'
          tmp_prof = 0.667 * (aa_prof + bb_prof)
          write(out_unit,'(A4,I2.2)')  't_i_',isp
          write(out_unit,'(5E16.9)')  tmp_prof !keV
      endif  
!
! Update Plasma State with new information - ps%nmini, ps%eperp_mini, ps%epll_mini
!
      CALL ps_store_plasma_state(ierr , trim(cur_state_file))
      if (ierr .ne. 0) stop 'cannot open state in perpare toric input'

!     end do
! PTB - ends
      deallocate(tmp_prof)
      deallocate(vol_int)
      deallocate(ns_tha)
      deallocate(ts_tha)
      deallocate(v_pars)
      deallocate(aa_prof)
      deallocate(bb_prof)
      deallocate(x_orig)
      deallocate(x_toric)
      close(out_unit)


   contains

      subroutine zone_check(rho, x_out, x_in)
    
      ! unravel zone points to boundary points

      REAL(KIND=rspec), intent(in) :: rho(:)  !nrho grid values
      REAL(KIND=rspec), intent(in) :: x_in(:)  !nrho-1 zone values
      REAL(KIND=rspec), intent(out) :: x_out(:)  !nrho boundary values
      
      integer :: irho,nrho
 
      nrho=size(x_out)
!assume there is a ghost point x_in(-1)=x_in(1), now we have nrho pts in x_in
!we still have no condition on the wall, so make f''=0 there
       x_out(nrho) = 0.5_rspec*(3.0_rspec * x_in(nrho -1) - x_in(nrho-2))
!midpoints values are average of mesh values:
      do irho = nrho - 1, 1, -1
         x_out(irho) = 2.0_rspec * x_in(irho) - x_out(irho + 1)
      end do

      end subroutine zone_check

      SUBROUTINE getlun (ilun,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Return an unused logical unit identifier.
!
!-----------------------------------------------------------------------
!
! ****** Upon successful completion (IERR=0), the first
! ****** unused logical unit number between MINLUN and
! ****** MAXLUN, inclusive, is returned in variable ILUN.
! ****** If all units between these limits are busy,
! ****** IERR=1 is returned.
!
!-----------------------------------------------------------------------
!
      INTEGER, INTENT(OUT) :: ilun, ierr
!
!-----------------------------------------------------------------------
!
! ****** Range of valid units.
!
      INTEGER, PARAMETER :: minlun=30, maxlun=99
      LOGICAL :: busy
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Find an unused unit number.
!
      DO i=minlun,maxlun
        INQUIRE (unit=i,opened=busy)
        IF (.NOT.busy) THEN
           ilun=1
           RETURN
        END IF
      END DO
!
! ****** Fall through here if all units are busy.
!
      ierr=1
      RETURN

      END subroutine getlun

      end program prepare_toric_input_mm
