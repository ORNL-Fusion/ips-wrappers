!SF cleaned up 2023. Removed comments see version control if you want
!comment history going back to '06

      program prepare_input
      USE plasma_state_mod
!--------------------------------------------------------------------------

    use swim_global_data_mod, only : &
            & rspec, ispec, &                ! int: kind specification for real and integer
            & swim_error                     ! error routine

      implicit none

      integer, parameter :: swim_string_length = 256  !for compatibility LAB
      real(rspec), parameter :: GHz=1.e9_rspec, cubic_cm=1.e-6_rspec
      integer, parameter :: isrc=1 
      
      integer :: ierr, i, isp, iwarn, irho, toric_to_alla(8) ! PTB
      logical :: lex
      
      character(len =swim_string_length) :: cur_state_file, cur_geq_file


      real(rspec), dimension(:),allocatable :: psi_poloidal_rho !DBB 6-27-2917
      real(rspec), dimension(:),allocatable :: psi_poloidal_eq !DBB 6-27-2917
      real(rspec), dimension(:),allocatable :: vol_rho !DBB 6-27-2917
      real(rspec), dimension(:),allocatable :: tmp_prof
      real(rspec), dimension(:),allocatable :: vol_int
      real(rspec), dimension(:,:),allocatable :: ns_tha
      real(rspec), dimension(:,:),allocatable :: ts_tha
      real(rspec), dimension(:,:),allocatable :: v_pars
      real(rspec), dimension(:),allocatable :: aa_prof
      real(rspec), dimension(:),allocatable :: bb_prof
      real(rspec), dimension(:),allocatable :: x_orig
      real(rspec), dimension(:),allocatable :: x_torlh
      real(rspec) :: tol_zero = 1.0E-12_rspec, dVol, dVol_int, Q_ps, Q_int

      !other variables for output of profiles
      integer :: nprodt, nproeq, kdiff_idens=1, kdiff_itemp=1, N_tmp
      !I/O units
      integer :: inp_unit, out_unit, iarg
      
      !program namelist variables
      !----------------------------------------------------------------
      character(10):: arg_toric_Mode = 'toric'
      character(10):: arg_inumin_Mode = 'Maxwell'
      integer:: arg_isol_Mode = 1
      real(rspec):: arg_enorm = 0.0_rspec
      character(16):: arg_specs = 'MIN'
      character(16):: arg_nphi = 'None' 

      !program namelist block
      !----------------------------------------------------------------
      namelist /toric_prepare_nml/ &
           arg_toric_mode, arg_inumin_mode, arg_isol_mode, &
           arg_enorm, arg_specs, arg_nphi

      
!------TORIC Namelist inputs-------------

      !/toricmode/ namelist variables
      !----------------------------------------------------------------
      character(16):: toricmode='toric'

      !/toricainp/ namelist variables
      !----------------------------------------------------------------
      !defaulted variables, these do not need to be changed
      !in most situations
      integer :: nvrb=3       ! Number of dimensions
      integer :: nptvac=-1    ! Number of points in vacuum
      integer :: mxmvac=15    ! Number of poloidal modes in vacuum
      integer :: iflr=1       ! flr effects
      integer :: ibpol=1      ! poloidal field on
      integer :: iqtor=1      ! toroidal magnetic broadening in Z functions
      integer :: icoll=0      ! No collisions
      integer :: imdedg=2     ! JPW pollution scheme for vacuum layer modelling
      integer :: iezvac=1     ! Ez suppressed in vacuum
      integer :: ibweld=1     ! eld damping of ibw's on
      integer :: icosig=0     ! extra damping for IBW when (kperp rho_i)^2>2
      integer :: iregax=1     ! turn on regularization at the magnetic axis
      integer :: isol=1       ! to solve or not to solve
      integer :: mastch=2     ! solver choice (=2 uses the new 3D parallel solver)
      integer :: lenwrd=8     !
      integer :: iout=0       ! 
      integer :: idlout=1     !
      integer :: io_ncdf=1    !
      integer :: iwdisk=1     ! output results for plotting in file
      integer :: IPLTHT=0     !
      integer :: iclres=0     ! adds collisions around isolated ion-ion resonances
      integer :: bscale=12    ! blocksize scaling factor for parallel runs only
      
      real(rspec) :: enhcol=1.0_rspec ! Collision enhancement factor
      real(rspec) :: dnures=1.0_rspec !width of layer for iclres (cm)
      real(rspec) :: tnures=1.0_rspec !strength of layer for iclres (A.U.)
      
      logical :: use_incore=.true. !solve matrix in RAM no saving to disk
      logical :: timing_on =.true. !internal time profiling
      
      !free variables, user must set in input deck or some set automagically by IPS
      INTEGER :: &
           nmod, ntt, nelm, nphi, pcblock

      REAL(rspec) :: &
           freqcy, antlen, antlc, theant, zeff

      integer, dimension(:) :: &
           inumin(0:nspmx) = 0
      
      !/qldciinp/ namelist variables
      !----------------------------------------------------------------
      integer, parameter :: max_runs=50 !max number of nphi used for dims.

      !defaulted variables
      integer :: subsmth=0 !Subgrid smoothing
      
      real(rspec) :: wdelta=-1.0_rspec ! Width of triangular delta
      real(rpsec) :: u_extr=10.0_rspec !

      logical :: ascii_out = .TRUE.
      logical :: ncdf_out  = .FALSE.

      character(80) :: path="./" ! Input data path

      !free variables
      INTEGER :: &
           num_runs, npsi_qld, ispec_Dql, iH_Dql,

      real(rspec) :: &
           d_u, rho_min, rho_max, enorm, deltapsi, pwtot,

      real(rspec), dimensions(:) :: &
           pw_nphi(max_runs)
      
      character(80) :: &
           ascii_name

      character(80), dimension(:) :: &
           files_toric(max_runs)
      
      !/equidata/ namelist variables
      !----------------------------------------------------------------
      integer, parameter :: nspmx = 8

      !defaulted variables
      integer :: igsmhd=1 !use equigs file for equil
      integer :: intchb=0 !use cubic spline interpolation
      integer :: idprof=1 !use equidt file for profs
      integer :: maxmod=0 !equilibrium simplification option off
      integer :: iudsym=0 !up down symmetry in equilibrium
      
      real(rspec) :: rtor=0.0_rspec   !geo variables that are generally unused
      real(rspec) :: rplasm=0.0_rspec !set to zero
      real(rspec) :: ashift=0.0_rspec
      real(rspec) :: aellip=0.0_rspec
      real(rspec) :: sellip=0.0_rspec
      real(rspec) :: strain=0.0_rspec
      real(rspec) :: azvert=0.0_rspec
      real(rspec) :: strigz=0.0_rspec
      real(rspec) :: sp_xs1=0.0_rspec
      real(rspec) :: bzero=0.0_rspec
      real(rspec) :: aicurr=0.0_rspec
      real(rspec) :: qaxis=0.0_rspec
      real(rspec) :: qedge=0.0_rspec
      real(rspec) :: ppjte=0.0_rspec
      real(rspec) :: ppjti=0.0_rspec
      real(rspec) :: denec=0.0_rspec
      real(rspec) :: tempec=0.0_rspec
      real(rspec) :: ppnei=0.0_rspec
      real(rspec) :: ppnee=0.0_rspec
      real(rspec) :: pptei=0.0_rspec
      real(rspec) :: pptee=0.0_rspec
      real(rspec) :: dnsepr=0.0_rspec
      real(rspec) :: tesepr=0.0_rspec
      real(rspec) :: so_thickness=0.0_rspec
      real(rspec) :: gldn=0.0_rspec
      real(rspec) :: glte=0.0_rspec
      
      real(rspec), dimension(:) :: aconc(nspmx)=0.0_rspec
      real(rspec), dimension(:) :: tempic(nspmx)=0.0_rspec
      real(rspec), dimension(:) :: pptii(nspmx)=0.0_rspec
      real(rspec), dimension(:) :: pptie(nspmx)=0.0_rspec
      real(rspec), dimension(:) :: tisepr(nspmx)=0.0_rspec
      real(rspec), dimension(:) :: glti(nspmx)=0.0_rspec

      character(80) :: inputpath='./'
      character(80) :: equil_file='equigs.data'
      character(80) :: profnt_file='equidt.data'

      !TORIC namelist blocks
      !----------------------------------------------------------------
      namelist /toric_mode/ &
           toricmode

      namelist /toricainp/ &
           nvrb,   nmod,   ntt,    nelm,   nptvac, mxmvac, &
           freqcy, nphi,   antlen, antlc,  theant, &
           iflr,   ibpol,  iqtor,  icoll,  enhcol, &
           imdedg, iezvac, ibweld, icosig, iregax, &
           isol,   mastch, lenwrd, iout,   idlout, io_ncdf, &
           iwdisk, ipltht, zeff,   iclres, dnures, tnures, &
           timing_on, scratchpath, bscale, use_incore, pcblock, &
           inumin

      namelist /qldciinp/ &
           num_runs, path, files_toric, npsi_qld, d_u, &
           rho_min, rho_max, enorm, wdelta, deltapsi, &
           ispec_Dql, iH_Dql, pwtot, pw_nphi, subsmth, &
           ascii_out, ncdf_out, ascii_name
      
      namelist/equidata/ &
           igsmhd, intchb, idprof, maxmod, &
           rtor,   rplasm, ashift, aellip, sellip, strian, &
           iudsym, azvert, strigz, sp_xs1, bzero,  aicurr, &
           qaxis,  qedge,  ppjte,  ppjti,  nspec,  mainsp, &
           atm,    azi,    aconc,  denec,  tempec, tempic, & 
           ppnei,  ppnee,  pptei,  pptee,  pptii,  pptie, &
           dnsepr, tesepr, tisepr, so_thickness, &   
           gldn,   glte,   glti, &  
           dist_plafars,   dist_plaant,    dist_plawall, &
           inputpath,      equil_file,     profnt_file &

      namelist /ips/ &
           toric_to_alla

      !
      ! begin program logic
      !

      ! Read in the namelist generated by the python driver
      ! -----------------------------------------------------------------------------
      call getlun(inp_unit,ierr)
      
      OPEN (unit=inp_unit, file = 'gacode_init.nml', status = 'old', &
           form = 'formatted', iostat = ierr)
  
      IF (ierr .ne. 0) THEN
         CALL SWIM_error ('open', 'generic_ps_init.f90',ps_init_nml_file)
         WRITE (*,*) 'generic_ps_init.f90: Cannot open ', TRIM(ps_init_nml_file)
         call exit(1)
      END IF

      READ(inp_unit, nml=ps_init_nml)
      CLOSE (inp_unit)
      WRITE (*, nml = ps_init_nml)

      if (arg_isol_Mode == 0) then
         isol = 0
      else if (arg_isol_Mode == 1) then
         isol = 1
      else
         write (*,*) 'prepare_toric_input_abr: unknown arg_isol_Mode = ', arg_isol_Mode
         stop
      end if

      ! Read the torica input file given by the ips
      ! -----------------------------------------------------------------------------
      call getlun(inp_unit,ierr)  ;  call getlun(out_unit,ierr)

      call ps_get_plasma_state(ierr,trim(cur_state_file))
      if(ierr .ne. 0) stop 'cannot get plasma state to get profiles '

!     nspec = ps%nspec_th !+1 !torlh includes electrons in species count
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

      write(*,*) 'Prepare toric input reading machine.inp'
      open(unit=inp_unit, file='machine.inp', status='old', &
              form='formatted')
      INQUIRE(inp_unit, exist=lex)
      IF (lex) THEN
         read(inp_unit, nml = qldciinp)
       	 read(inp_unit, nml = toricainp)
         read(inp_unit, nml = equidata)
      ELSE
         write(*,*) &
            'machine.inp does not exist or there was a read error'
      END IF
      close(inp_unit)
      IF (nspec > nspmx ) THEN
         write(*,*) "Error, nspec > nspmx in toric, reducing to nspmx"
         nspec=nspmx
      END IF

! Overwrite freqcy, pwtot, and enorm found in machine.inp

      freqcy = ps%freq_ic(isrc) !Hz
      pwtot =  ps%power_ic(isrc) !watts
      
      if (trim(arg_enorm) /= 'None') then
	 read(arg_enorm,*) enorm
      end if

      if (trim(arg_specs) /= 'MIN') then
         arg_specs = trim(arg_specs)
      endif
      
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
! DBB 9-28-2016 in plasma state Zeff is profile.  If ps%Zeff = 0 use default or value from
! toricainp namelist.  If ps%zeff(1) is not zero use that
      if (ps%zeff(1) .ne. 0.) then
          zeff = ps%zeff(1)
      endif
      atm(1:ps%nspec_tha) = NINT(ps%m_ALLA(1:ps%nspec_tha)/ps_mp) 
      azi(1:ps%nspec_tha) = NINT(ps%q_ALLA(1:ps%nspec_tha)/ps_xe)

      !SF Set first species to general always should be this in current
      !version
      if ((trim(arg_specs).eq.'MIN+').and. &
          (trim(arg_inumin_mode) .eq. 'nonMaxwell'))then
         inumin(0) = 1 !turn on abj mode for general
      endif
      toricmode = trim(arg_toric_Mode)

      !SF set up dql in cases that multiple non-Maxwellian species
      !are being evolved
      if (toricmode == 'qldci1') then
         toricmode = 'qldci_CQL3D'
         iH_dql = 1
         ispec_dql = 3
         ascii_name = 'du0u0_input_1'
      elseif (toricmode=='qldci2')then
         toricmode = 'qldci_CQL3D'
         iH_dql = 2
         ispec_dql = 1
         ascii_name = 'du0u0_input_2'
      elseif (toricmode=='qldci')then
         toricmode = 'qldci_CQL3D'
         ascii_name = 'du0u0_input'
      endif
      
      
! JCW create torlh to alla mapping in file:toric_alla_map.txt
!     First nspec_tha entries are just 1:nspec_tha
!     Then add index for non-therma species
      isp = ps%nspec_tha
      toric_to_alla(1:isp)=(/ (i, i=1,isp) /)
      if (allocated(ps%nbeami)) then
         isp = isp + 1
         toric_to_alla(isp)=ps%snbi_to_alla(1)
         atm(isp) = NINT(ps%m_ALLA(ps%snbi_to_alla(1))/ps_mp)
         azi(isp) = NINT(ps%q_ALLA(ps%snbi_to_alla(1))/ps_xe)
      endif

      if (allocated(ps%nmini)) then
         isp = isp + 1
         toric_to_alla(isp)=ps%rfmin_to_alla(1)
         atm(isp) = NINT(ps%m_ALLA(ps%rfmin_to_alla(1))/ps_mp)
         azi(isp) = NINT(ps%q_ALLA(ps%rfmin_to_alla(1))/ps_xe)
         if(trim(arg_inumin_mode) .eq. 'nonMaxwell')then
            inumin(isp-1) = 1 !SF turn on abj mode for minority species
         end if
      endif

      if (allocated(ps%nfusi)) then
         isp = isp + 1
         toric_to_alla(isp)=ps%sfus_to_alla(1)
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

      open(unit=out_unit, file='torica.inp',                &
        status = 'unknown', form = 'formatted',delim='quote')

      WRITE (*,*) 'toricmode = ', toricmode
      WRITE (*,*) 'inumin = ', inumin
      WRITE (*,*) 'isol = ', isol
      WRITE (*,*) 'pwtot = ', pwtot
      WRITE (*,*) 'enorm = ', enorm
      WRITE (*,*) 'nphi = ', nphi

      write(out_unit, nml = toric_mode)
      write(out_unit, nml = qldciinp) !SF switched to qldciinp
      write(out_unit, nml = toricainp)
      write(out_unit, nml = equidata)
      write(out_unit, nml = ips)

      close(out_unit)


      nprodt = size(ps%rho)
      nproeq = size(ps%rho_eq)

! PTB begins -
      allocate( psi_poloidal_rho(nprodt-1))  !DBB 6-27_2017
      allocate( psi_poloidal_eq(nproeq))  !DBB 6-27_2017
      allocate( vol_rho(nprodt-1))  !DBB 6-27_2017
      allocate( tmp_prof(nprodt))
      allocate( vol_int(nprodt))
      allocate( ns_tha(nprodt-1,ps%nspec_tha))
      allocate( ts_tha(nprodt-1,ps%nspec_tha))
      allocate( v_pars(nprodt-1,ps%nspec_tha))
      allocate( aa_prof(nprodt))
      allocate( bb_prof(nprodt))
      allocate( x_orig (nprodt-1))
      allocate( x_torlh (nprodt))
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
! DBB begins: First interpolate psi_poloidal from rho_eq grid to rho grid.  Note, if you
!             try to interpolate ps%psipol with the ps_user_1dintrp_vec routine it will
!              detect that ps%psipol is a zone centered variable and will interpolate it
!              to length nrho-1. So fool it by introducing dummy non-PS varible psi_poloidal_eq

	  psi_poloidal_eq = ps%psipol
      call ps_user_1dintrp_vec(ps%rho, ps%rho_eq, psi_poloidal_eq, x_torlh, ierr )
      x_torlh = sqrt(x_torlh/x_torlh(nprodt))
	  !write (*,*) " "
	  !write (*,*) "x_torlh = "
	  !write (*,*) x_torlh
! DBB ends

!     write(out_unit,'(A10)')  'rho'
!     write(out_unit,'(5E16.9)')  ps%rho !check units

! PTB begins
! Need to actually write out the sqrt (Psi_pol) - normalized
      write(out_unit,'(A20)')  'Sqrt(Psi / Psilim)'
      write(out_unit,'(5E16.9)') x_torlh
! PTB ends

      write(out_unit,'(A10)')  'n_e'
      !write(*,*) 'ps%ns(1,0)', ps%ns(1,0)
!
! Interpolate the electron density profile from the Plasma State grid to the Toric grid
! N.B. Toric grid maps directly to rho grid ps%rho, same radii.
         call ps_user_1dintrp_vec(ps%rho, x_orig, ps%ns(:,0), tmp_prof(:),ierr )  !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS electron density profile onto Toric grid'
!
!
         call ps_user_1dintrp_vec(ps%rho, ps%rho_eq, ps%vol(:), vol_int(:),ierr ) !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS volume onto Toric grid'
	  !write (*,*) " "
	  !write (*,*) "interpolated vol_int = "
	  !write (*,*) vol_int

!
! PTB Compare the volume average of the orginal density profile from the Plasma State
! and the volume average of the interpolated profile
!
         Q_ps = 0.0_rspec
         Q_int = 0.0_rspec
         do irho = 1,nprodt-1
            dVol = ps%vol(irho+1) - ps%vol(irho)
            Q_ps = Q_ps + 0.5_rspec * (ps%ns(irho,0) + &
                   ps%ns(irho,0)) * dVol
         end do
         do irho = 1,nproeq-1
            dVol_int = vol_int(irho+1) - vol_int(irho)
            Q_int = Q_int +  0.5_rspec * (tmp_prof(irho) + &
                   tmp_prof(irho)) * dVol_int
         end do
         Q_ps = Q_ps / ps%vol(nprodt)
         Q_int = Q_int / ps%vol(nprodt)
      write(*,*) '<n_e(m-3)-PS> =',  Q_ps
      write(*,*) '<n_e(m-3)-Interpolated> =', Q_int
      write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3

      write(out_unit,'(A10)')  't_e'
         call ps_user_1dintrp_vec(ps%rho, x_orig, ps%Ts(:,0), tmp_prof(:),ierr ) !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS electron temperature profile onto Torlh grid'
      write(out_unit,'(5E16.9)')  tmp_prof   !keV
	  !write (*,*) " "
	  !write (*,*) "tmp_prof = "
	  !write (*,*) tmp_prof
! PTB - begins
        call ps_tha_fetch (ierr, ps, tol_zero, &
             ns_tha, ts_tha, v_pars, iwarn)
         if(ierr .ne. 0) stop 'error getting the density and temperature data from the abridged species list'
      do isp=1,ps%nspec_tha
         write(out_unit,'(A4,I2.2)')  'n_i_',isp
         call ps_user_1dintrp_vec(ps%rho, x_orig, ns_tha(:,isp), tmp_prof(:),ierr ) !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS ion density profile onto Torlh grid'
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3
	   !  write (*,*) " "
	   !  write (*,*) "density_prof = "
	   !  write (*,*) tmp_prof

         write(out_unit,'(A4,I2.2)')  't_i_',isp
         !write(*,*) "Thermal ion name, A, Z, dens, temp = "
         !write(*,*) trim(ps%alla_name(isp)), atm(isp),azi(isp),ns_tha(1,isp),ts_tha(1,isp)
         call ps_user_1dintrp_vec(ps%rho, x_orig, Ts_tha(:,isp), tmp_prof(:),ierr )  !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS ion temperature profile onto Torlh grid'
         write(out_unit,'(5E16.9)')  tmp_prof !keV
	    ! write (*,*) " "
	    ! write (*,*) "temperature_prof = "
	    ! write (*,*) tmp_prof
      end do
! PTB - ends

! PTB - begins
!
! Move on to the non-thermal ion species - in order of NBI, RF minority, and fusion alphas
!
! (1) Write the neutral beam injection data to the Torlh equidt.data file:
!
! First interpolate the NBI density [nbeami(:,1)] from the NuBeam grid onto the Torlh input data
! grid. Then write it to the equidt.data file.
!
      if(allocated(ps%nbeami)) then
       isp = ps%snbi_to_alla(1)
         write(*,*) "Fast ion name, A, Z = ", trim(ps%alla_name(isp)), &
     &              NINT(ps%m_alla(isp)/ps_mp), NINT(ps%q_alla(isp)/ps_xe)
         call ps_user_1dintrp_vec(ps%rho, ps%rho_nbi, ps%nbeami(:,1), tmp_prof(:),ierr )  !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS NBI density profile onto Torlh grid'
         write(out_unit,'(A4,I2.2)')  'n_i_',isp
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3
	     write (*,*) " "
	     write (*,*) "beam density profile = "
	     write (*,*) tmp_prof

!
! Next interpolate the NBI energies [eperp_beami(:,1) and epll_beami(:,1)] from the NuBeam grid
! onto the Torlh input data grid. Then compute the equivalent temperature profile for the NBI and
! write it to the equidt.data file.
!
          call ps_user_1dintrp_vec(ps%rho, ps%rho_nbi, ps%eperp_beami(:,1), &
               aa_prof(:),ierr )  !DBB 6-27_2017
          if(ierr .ne. 0) stop 'error interpolating PS NBI perp. energy profile onto Torlh grid'
          call ps_user_1dintrp_vec(ps%rho, ps%rho_nbi, ps%epll_beami(:,1),  &
               bb_prof(:),ierr )  !DBB 6-27_2017
          if(ierr .ne. 0) stop 'error interpolating PS NBI parallel energy profile onto Torlh grid'
          tmp_prof = 0.667 * (aa_prof + bb_prof)
          write(out_unit,'(A4,I2.2)')  't_i_',isp
          write(out_unit,'(5E16.9)')  tmp_prof !keV
 	      write (*,*) " "
	      write (*,*) "beam temperature profile = "
	      write (*,*) tmp_prof

      endif
!
! (2) Write the rf minority tail data to the Torlh equidt.data file:
!

      
      if(allocated(ps%nmini)) then
       isp = ps%rfmin_to_alla(1)
         write(*,*) "Fast ion name, A, Z = ", trim(ps%alla_name(isp)), &
     &              NINT(ps%m_alla(isp)/ps_mp), NINT(ps%q_alla(isp)/ps_xe)

!
! If (kdens_rfmin .EQ. 'fraction') then assume PS data override for nmini and instead
! compute nmini = fracmin * ne
!
        if (trim(ps%kdens_rfmin) .EQ. 'fraction') then
           
         call ps_user_1dintrp_vec(x_orig, ps%rho, ps%ns(:,0), &
     &          tmp_prof(:),ierr ) !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS electron density profile onto Torlh grid'
         tmp_prof(:) = ps%fracmin(1)*tmp_prof(:)
         write(out_unit,'(A4,I2.2)')  'n_rfmin_',isp
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3
 	      !write (*,*) " "
	      !write (*,*) "rf minority density = "
	      !write (*,*) tmp_prof

        !write the minority fraction profile into the PS
        call ps_user_1dintrp_vec(ps%rho,x_orig,tmp_prof, &
               ps%nmini(:,1),ierr ) !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating new minority desnity profile onto PS grid'
        endif
 	      !write (*,*) " "
	      !write (*,*) "rf minority density put back into plasma state = "
	      !write (*,*) ps%nmini(:,1)

!
! If (kdens_rfmin .EQ. 'data') then assume nmini is available in the PS, read it, and interpolate
! it from the rho-icrf grid onto the Toric radial grid
!
        if (trim(ps%kdens_rfmin) .EQ. 'data') then
         call ps_user_1dintrp_vec(ps%rho, ps%rho_icrf, ps%nmini(:,1), &
               tmp_prof(:),ierr ) !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS minority density profile onto Toric grid'
         write(out_unit,'(A4,I2.2)')  'n_rfmin_',isp
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3
  	      write (*,*) " "
	      write (*,*) "rf minority density interpolated from plasma state = "
	      write (*,*) tmp_prof
       endif
!
! If isThermal=1 then assume the RF minority tail temperature is equal to the bulk ion temperature.
! Set the minority tail temperature equal to the temperature profile of the first bulk ion species
! of the abridged species list, by interpolating that array from the PS grid onto the Torlh input data
! grid, and then writing it to the equidt.data file.
!
         write(out_unit,'(A4,I2.2)')  't_rfmin_',isp
         if (ps%isThermal(1) .eq. 1) then
			 call ps_user_1dintrp_vec(x_orig, ps%rho, Ts_tha(:,1), &
				   tmp_prof(:),ierr ) !DBB 6-27_2017
			 if(ierr .ne. 0) stop 'error interpolating PS RF minority temperature profile onto Toric grid'
			 write(out_unit,'(5E16.9)')  tmp_prof !keV
			  !write (*,*) " "
			  !write (*,*) "rf minority tail temperature interpolated from Ts_tha(:,1) = "
			  !write (*,*) tmp_prof
         endif
!
! If isThermal=2 then assume the RF minority tail temperature data is available in the PS. In this
! case just interpolate the ICRF minority energies [eperp_mini(:,1) and epll_mini(:,1)] from the ICRF
! grid (rho_lhrf) onto the Torlh input data grid. Then compute the equivalent temperature profile for
! the RF minority profile. Finally write this array to the equidt.data file.
!
         if (ps%isThermal(1) .eq. 2) then
          call ps_user_1dintrp_vec(x_orig, ps%rho_icrf, ps%eperp_mini(:,1), &
               aa_prof(:),ierr )  !DBB 6-27_2017
          if(ierr .ne. 0) stop 'error interpolating PS RF minority perp. energy profile onto Torlh grid'
          call ps_user_1dintrp_vec(x_orig, ps%rho_icrf, ps%epll_mini(:,1),  &
               bb_prof(:),ierr )  !DBB 6-27_2017
          if(ierr .ne. 0) stop 'error interpolating PS RF minority parallel energy profile onto Torlh grid'
          tmp_prof = 0.667 * (aa_prof + bb_prof)
          write(out_unit,'(5E16.9)')  tmp_prof !keV
			  write (*,*) " "
			  write (*,*) "rf minority tail temperature interpolated from plasma state = "
			  write (*,*) tmp_prof
         endif
      endif
!
! (3) Write the fast fusion alpha data to the Torlh equidt.data file:
!
! First interpolate the alpha density [nfusi(:,1)] from the fusion grid onto the Torlh input data
! grid. Then write it to the equidt.data file.
!
      if(allocated(ps%nfusi)) then
       isp = ps%sfus_to_alla(1)
         write(*,*) "Fast ion name, A, Z = ", trim(ps%alla_name(isp)), &
     &              NINT(ps%m_alla(isp)/ps_mp), NINT(ps%q_alla(isp)/ps_xe)
         call ps_user_1dintrp_vec(x_orig, ps%rho_fus, ps%nfusi(:,1), &
               tmp_prof(:),ierr )  !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS fusion alpha density profile onto Torlh grid'
         write(out_unit,'(A4,I2.2)')  'n_i_',isp
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3
		  write (*,*) " "
		  write (*,*) "nfusi interpolated from plasma state = "
		  write (*,*) tmp_prof

!
! Next interpolate the alpha energies [eperp_fusi(:,1) and epll_fusi(:,1)] from the fusion grid
! onto the Torlh input data grid. Then compute the equivalent temperature profile for the NBI and
! write it to the equidt.data file.
!
          call ps_user_1dintrp_vec(x_orig, ps%rho_fus, ps%eperp_fusi(:,1), &
               aa_prof(:),ierr )  !DBB 6-27_2017
          if(ierr .ne. 0) stop 'error interpolating PS fusion alpha perp. energy profile onto Torlh grid'
          call ps_user_1dintrp_vec(x_orig, ps%rho_fus, ps%epll_fusi(:,1),  &
               bb_prof(:),ierr )
          if(ierr .ne. 0) stop 'error interpolating PS fusion alpha parallel energy profile onto Torlh grid'
          tmp_prof = 0.667 * (aa_prof + bb_prof)
          write(out_unit,'(A4,I2.2)')  't_i_',isp
          write(out_unit,'(5E16.9)')  tmp_prof !keV
		  write (*,*) " "
		  write (*,*) "alpha energy interpolated from plasma state = "
		  write (*,*) tmp_prof

      endif
!
! Update Plasma State with new information - ps%nmini, ps%eperp_mini, ps%epll_mini
!
      CALL ps_store_plasma_state(ierr , trim(cur_state_file))
      if (ierr .ne. 0) stop 'cannot open state in perpare torlh input'

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
      deallocate(x_torlh)
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


    SUBROUTINE ckerr(sbrtn)
      character*(*), intent(in) :: sbrtn

      IF(ierr.NE.0) then
         write(6,*) ' ?plasma_state_test: error in call: '//trim(sbrtn)
         stop
      ENDIF
    END SUBROUTINE ckerr


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

      end program prepare_input
