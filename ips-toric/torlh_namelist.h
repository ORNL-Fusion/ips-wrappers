!------Namelist inputs-------------
! see the rf component write-up in the svn repository for more description of
! the inputs

!  Mode of TORIC utilization "equil" to preprocess equilibrium, "toric" to run solver
      character(10):: toricmode="toric"

! Dimensions of the problem
      integer :: nvrb=3       ! Generally three vector components
! Poloidal resolution
      integer :: ntt = 512
      integer :: nmod =255     ! nmod will be set as nmod=ntt/2-1,
                              ! unless user enters a non-zero positive value
                              ! less than ntt/2-1 to conserve memory at
                              ! the cost of efficiency.
! Radial resolution: Number of radial elements in the plasma
      integer :: nelm = 540
!N.B. max memory usage for double precision is approx 16*3*(6*Nmod)^2*Nelm/Nprocs
! on a machine with 2GB of ram, the practical in core limit in Nm=255, Nelm=240, but
! this run will take several hours, so most likely, memory will not be a problem since
! many processors will be used to make run times reasonable

!Don't change these three unless you really know what you're doing
! Number of points in vacuum
      integer :: nptvac = 5
! Number of poloidal modes in vacuum
      integer :: mxmvac = 63
! Number of gaussian quadrature points in radial finite elements
      integer :: ngaus = 3

! Namelist inputs for choice of the model, rarely need to be changed
      integer ::   isol=1     ! to solve or not to solve
      integer ::   mastch=2   ! solver choice
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
      logical ::   use_incore=.false.
                              ! default is to use out-of-core memory
      integer ::   pcblock=4  !  number of processors used in poloidal pc mesh when MASTCH=2
! Namelist inputs for control of the output, most are off to avoid too much data dumped
      integer ::   iout=0, ipltht=0, idlout=1
! NetCDF output for plots on by default.  torica.sol file has solutions but not necessarily
! many of the parameters
      integer ::   io_ncdf = 1

! Wave and antenna parameters (default values for now, later
!   use ps%ant_model file for machine state)
!      integer     :: nphi=10    ! anzedg is used for TORLH
      real     :: anzedg = -1.6 ! toroidal refractive index at the edge
      real(rspec) :: freqcy=4.6e9_rspec
      integer  :: ibcant =1     !boundary condition for antenna. When ibcant<0, the grill antenna is used

!units are in cm for lengths
      real(rspec) :: antlen=6.0_rspec, antlc=1.0_rspec, &
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
      character(24) :: equil_file='"equigs.data"', profnt_file='equidt.data', &
                       qlmin_file='', gfile
! Path to the file containing the namelists
      character(24) :: inputpath = '"./"'
! Path to the buffer file for the solver using virtual memory
      character(24) :: scratchpath = '"/tmp"'

! Namelist for numerical magnetic equilibrium
      integer :: nmhd, intchb=7

! Namelist inputs for positioning the plasma with respect to the wall
      real(rspec), dimension(:) :: dist_plawall(4)=        &
     &                         (/10._rspec,10._rspec,10._rspec,10._rspec/)
      real(rspec) ::  dist_plafars=0._rspec, dist_plaant=0.5_rspec
      real(rspec) ::  so_thickness=0.5_rspec

! Namelist inputs for profile specification
!  Density and temperature profiles - namelist and file inputs
!
      integer, parameter :: nspmx =30
! NOTE: A maximum of 15 ion species allowed
! The place,nspec+1, is reserved for the electrons in mod_direl.
! mainsp is used to impose charge neutrality
      integer ::   nspec,  mainsp, iudsym=0

      integer, dimension(:) ::                                       &
     &             atm(nspmx)=0._rspec,       azi(nspmx)=0._rspec

      real(rspec) :: q_rfmin = 1.0_rspec, qatom_rfmin = 1.0_rspec, &
     &               m_rfmin = 1.0_rspec, fracmin = 0.04_rspec
      integer :: isThermal = 1
      character(24) :: rfmin_name = 'H_min', kdens_rfmin = 'fraction'

! Note that for variables like rfmin_name and kdens_rfmin that will be put into the
! Plasma State that we do not want to use the construct of rfmin_name = '"H_min"' since
! this will try to place the string "H_min" inot the PS variable and double quotes
! apparently are not allowed !!!


      real(rspec), dimension(:) ::                                   &
     &             aconc(nspmx)=0._rspec,     tempic(nspmx)=0._rspec,&
     &             tisepr(nspmx)=0._rspec,    glti(nspmx)=0._rspec,  &
     &             pptii(nspmx)=0._rspec,     pptie(nspmx)=0._rspec

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
!     &   nvrb,   nmod,   ntt,    nelm,   nptvac, mxmvac, &
!     &   freqcy, anzedg, ibcant, antlen, antlc,  theant, &
!     &   iflr,   ibpol,  iqtor,  icoll,  enhcol, &
!     &   imdedg, iezvac, ibweld, icosig, iregax, &
!     &   isol,   mastch,         iout,   idlout, io_ncdf, &
!     &   iwdisk, ipltht, zeff,   iclres, dnures, tnures, &
!     &   timing_on, scratchpath, bscale, use_incore, pcblock
     &   nmod,   ntt,    nelm,   nptvac, mxmvac, &
     &   freqcy, anzedg, ibcant, antlen, antlc,  theant, &
     &   iflr,   ibpol,  icoll,  enhcol, &
     &   iregax, &
     &   isol,   mastch,         iout,   idlout, &
     &   iwdisk, zeff, &
     &   timing_on, scratchpath, use_incore, pcblock, inputpath

      namelist /nonthermals/ &
     &   fracmin, q_rfmin, qatom_rfmin, m_rfmin, rfmin_name, &
     &   kdens_rfmin, isThermal


!originally in t4_mod_toi2mex.F
!specifies numerical equilibrium (EFIT usually) settings
      namelist /gfreadinp/ &
     &   gfile, inputFormat, nmhd, ntheta, &
     &   i2mex_remap, i2mex_mx, i2mex_npsi, i2mex_kb, i2mex_direct, &
     &   i2mex_last_norm_surface_is, ic1, inc, imom1

!originally in t4_mod_equil.F
!specifies plasma equilibrium settings and profiles
!      namelist/equidata/ igsmhd, intchb, idprof, &
!     &   nspec,  mainsp, &
!     &   atm,    azi, so_thickness, iudsym,&
!     &   dist_plafars,   dist_plaant,    dist_plawall, &
!     &   equil_file,     profnt_file

      namelist/equidata/ igsmhd, intchb, idprof, &
     &   nspec,  mainsp, &
     &   atm,    azi, so_thickness, &
     &   dist_plafars,   dist_plaant,    dist_plawall, &
     &   equil_file,     profnt_file