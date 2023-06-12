! SF cleaned up 2023. Removed comments see version control if you want
! comment history going back to '06. This is a much simpler version of
! the routine that is going to be a lot easier to run, and all the
! functionality in this one actually works.

      program prepare_toric_input
      USE plasma_state_mod
      use swim_global_data_mod, only : &
            & rspec, ispec, &                ! int: kind specification for real and integer
            & swim_error                     ! error routine

      implicit none

      !program params
      !---------------------------------------------------------------
      integer, parameter :: swim_string_length = 256  !for compatibility LAB
      real(rspec), parameter :: cubic_cm=1.e-6_rspec
      integer, parameter :: nspmx = 8
      integer, parameter :: max_runs=50 !max number of nphi used for dims.

      !program variables
      !---------------------------------------------------------------
      integer :: ierr, i, isp, iwarn, irho, minspec, pairspec, toric_to_alla(8) ! PTB
      logical :: lex

      !other variables for output of profiles
      integer :: nprodt, nproeq, kdiff_idens=1, kdiff_itemp=1
      !I/O units
      integer :: inp_unit, out_unit
      
      ! ptb modification arrays for equid_output
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
      real(rspec), dimension(:),allocatable :: x_toric
      real(rspec) :: tol_zero = 1.0E-12_rspec, dVol, dVol_int, Q_ps, Q_int

      character(16):: outfile='torica.inp'
      
      !program namelist variables
      !----------------------------------------------------------------
      character(len =swim_string_length) :: cur_state_file, cur_geq_file
      character(10):: arg_toric_Mode = 'toric'
      character(10):: arg_inumin_Mode = 'Maxwell'
      real(rspec):: arg_enorm = 0.0_rspec
      character(16):: arg_specs = 'MIN'
      integer:: arg_src_indx = 0
      integer:: arg_nphi_indx = 0 !if set to not 0 override nphi settings in machine.inp
      integer:: force_defaults = 0 !idiot proofing setting forces cleanup of input file
      integer:: arg_nfpsurf = 0
      real(rspec):: arg_fplo=0.1_rspec
      real(rspec):: arg_fphi=0.99_rspec
      
      !program namelist block
      !----------------------------------------------------------------
      namelist /toric_prepare_nml/ &
           cur_state_file, cur_geq_file, &
           arg_toric_mode, arg_inumin_mode, &
           arg_enorm, arg_specs, arg_src_indx, arg_nphi_indx, &
           arg_nfpsurf, arg_fplo, arg_fphi, &
           force_defaults
      
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
           inumin(nspmx) = 0

      character(80) :: scratchpath="~/pscratch/tmp/"
      
      
      !/qldciinp/ namelist variables
      !----------------------------------------------------------------
      !defaulted variables
      integer :: subsmth=0 !Subgrid smoothing
      
      real(rspec) :: wdelta=-1.0_rspec ! Width of triangular delta
      real(rspec) :: u_extr=10.0_rspec !

      logical :: ascii_out = .TRUE.
      logical :: ncdf_out  = .FALSE.

      character(80) :: path="./" ! Input data path
      
      !free variables
      INTEGER :: &
           num_runs, npsi_qld, ispec_Dql, iH_Dql

      real(rspec) :: &
           d_u, rho_min, rho_max, enorm, deltapsi, pwtot

      real(rspec), dimension(:) :: &
           pw_nphi(max_runs) = 0.0
      
      character(80) :: &
           ascii_name = ''

      character(80), dimension(:) :: &
           files_toric(max_runs) = ''
      
      !/equidata/ namelist variables
      !----------------------------------------------------------------
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
      real(rspec) :: strian=0.0_rspec
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

      !free variables
      integer :: &
           nspec, mainsp

      real(rspec) :: &
           dist_plafars, dist_plaant
      
      real(rspec), dimension(:) :: &
           atm(nspmx), azi(nspmx), dist_plawall(4)
      

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
           inputpath,      equil_file,     profnt_file 

      namelist /ips/ &
           toric_to_alla

      
      ! --- Begin Program Logic ---
      
      ! Read in the namelist generated by the python driver
      ! -----------------------------------------------------------------------------
      call getlun(inp_unit,ierr)
      
      OPEN (unit=inp_unit, file = 'toric_prepare.nml', status = 'old', &
           form = 'formatted', iostat = ierr)
  
      IF (ierr .ne. 0) THEN
         CALL SWIM_error ('open', 'prepare_toric_input.f90','toric_prepare.nml')
         WRITE (*,*) 'prepare_toric_input.f90: Cannot open ', 'toric_prepare.nml'
         call exit(1)
      END IF

      READ(inp_unit, nml=toric_prepare_nml)
      CLOSE (inp_unit)
      
      WRITE (*, nml = toric_prepare_nml)

      ! Get the plasma state
      ! -----------------------------------------------------------------------------
      
      call ps_get_plasma_state(ierr,trim(cur_state_file))
      if(ierr .ne. 0) stop 'cannot get plasma state to get profiles '

      
      ! Read the torica input file given by the ips
      ! -----------------------------------------------------------------------------
      call getlun(inp_unit,ierr)  ;  call getlun(out_unit,ierr)

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

      ! Modify namelist values based on the plasma state
      ! -----------------------------------------------------------------------------

      !defaults that should always be forced
      !toricainp
      nvrb=3; iflr=1; ibpol=1; isol=1; mastch=2;
      lenwrd=8; iout=0; idlout=1; io_ncdf=1
      iwdisk=1; use_incore=.true.
      
      !qldciinp
      ascii_out=.true.

      !equidinp
      igsmhd=1; idprof=1; maxmod=0; iudsym=0
      equil_file='equigs.data'
      profnt_file='equidt.data'

      ! if end user is an idiot and provides bad input file
      ! force_defaults = 1 should save the day
      if (force_defaults > 0)then
         !force toricainp values
         nvrb=3; nptvac=-1; mxmvac=15; iflr=1; ibpol=1; iqtor=1
         icoll=0; imdedg=2; iezvac=1; ibweld=1; icosig=0; iregax=1
         isol=1; mastch=2; lenwrd=8; iout=0; idlout=1; io_ncdf=1
         iwdisk=1; IPLTHT=0; iclres=0; bscale=12
         use_incore=.true.
         !force qldciinp values
         subsmth=0; wdelta=-1.0_rspec; u_extr=10.0_rspec
         ascii_out = .TRUE.
         ncdf_out  = .FALSE.
         path="./"
         !force equidt values
         igsmhd=1; intchb=0; idprof=1; maxmod=0; iudsym=0
         so_thickness=0.0_rspec
         inputpath='./'
         equil_file='equigs.data'
         profnt_file='equidt.data'
      endif

      ! if end user really screwed up nuke resolution and physics settings
      ! in the input file. these settings here generally work for minority
      ! heating when run on 1 node on nersc perlmutter
      ! force_defaults = 2
      if (force_defaults > 1)then
         !force toricainp values
         nmod = 127; ntt=256; nelm=500; pcblock=4
         !force qldciinp values
         d_u = 0.005
         !force equidt values
         dist_plafars=0.0; dist_plaant=1.0; dist_plawall(:)=2.5
      endif

      !set toric mode
      toricmode = trim(arg_toric_Mode)
      
      ! overwrite icrf src dependent values in machine.inp
      ! if using multiple sources makes things /a lot/ easier
      if(arg_nphi_indx.gt.0)then
         !toricainp
         nphi = ps%nphi(1,arg_nphi_indx) 
      endif

      if(arg_src_indx.gt.0)then
         !toricainp
         freqcy = ps%freq_ic(arg_src_indx)
         if (arg_nphi_indx.gt.0)then
            nphi = ps%nphi(arg_src_indx,arg_nphi_indx)
         endif
      endif
      
      !qldciinp
      num_runs = ps%num_nphi(1)
      pwtot =  ps%power_ic(1) !watts
      pw_nphi(1:num_runs) = ps%wt_nphi(1,1:num_runs)
      if (num_runs.gt.1)then
         do i=1,num_runs
            WRITE(files_toric(i),"(A,I2.2)") "fort.9_",i
         enddo
      else
         files_toric(1) = "fort.9"
         pw_nphi(1) = 1.0
      endif
      
      ! overwrite enorm in qldciinp
      if (arg_enorm /= 0.0_rspec) then
	 enorm = arg_enorm
      else
         enorm=1000.0_rspec
      end if

      ! overwrite fp surfaces in qldciinp
      if (arg_nfpsurf.ne.0)then
         npsi_qld = arg_nfpsurf
         rho_min  = arg_fplo
         rho_max  = arg_fphi
      endif

      ! overwrite deltapsi in qldciinp
      deltapsi = ps%psipol(ps%nrho_eq)
      
      ! overwrite Zeff (not particularily important in TORIC so just set naively)
      if (ps%zeff(1) .ne. 0.) then
          zeff = ps%zeff(1)
      endif

      ! set mass/charge based on plasma state values
      toric_to_alla = 0
      isp = ps%nspec_tha

      do i=1,isp
         toric_to_alla(i) = i
      enddo
      atm(1:isp) = 1.0_rspec*NINT(ps%m_ALLA(1:ps%nspec_tha)/ps_mp) 
      azi(1:isp) = 1.0_rspec*NINT(ps%q_ALLA(1:ps%nspec_tha)/ps_xe)
      
      if (ps%nspec_rfmin>0) then
         isp = isp + 1
         minspec = isp
         toric_to_alla(isp)=ps%rfmin_to_alla(1)
         atm(isp) = NINT(ps%m_ALLA(ps%rfmin_to_alla(1))/ps_mp)
         azi(isp) = NINT(ps%q_ALLA(ps%rfmin_to_alla(1))/ps_xe)
      endif

      if (ps%nspec_fusion>0) then
         isp = isp + 1
         toric_to_alla(isp)=ps%sfus_to_alla(1)
         atm(isp) = NINT(ps%m_ALLA(ps%sfus_to_alla(1))/ps_mp)
         azi(isp) = NINT(ps%q_ALLA(ps%sfus_to_alla(1))/ps_xe)
      endif

      nspec = isp
      mainsp = 1
      ! if necessary reduces species count
      IF (nspec > nspmx ) THEN
         write(*,*) "Error, nspec > nspmx in toric, reducing to nspmx"
         nspec=nspmx
      END IF

      ! find the minority's bulk pair (defined q_min/m_min = 2 q_pair/m_pair)
      ! prevent accidentally tripping on impurity species by adding axi<6 req.
      pairspec = 0
      do i = 1,ps%nspec_tha
         if ((2.0_rspec*azi(i)/atm(i)==azi(minspec)/atm(minspec)).and.(azi(i)<6)) then
            pairspec = i
         endif
      enddo
      IF (pairspec == 0) WRITE(*,*) 'No minority/bulk pair found'

      ! if nonmaxwellian turn on minority nonmax dist
      if ((trim(arg_specs).eq.'MIN').and. &
           (trim(arg_inumin_mode) .eq. 'nonMaxwell'))then
         inumin(minspec) = 1
      endif
      
      ! turn on nonmax dielectric for both bulk and 
      ! minority
      if ((trim(arg_specs).eq.'MIN+').and. &
           (trim(arg_inumin_mode) .eq. 'nonMaxwell'))then
         inumin(minspec) = 1
         inumin(pairspec) = 1
      endif
     
      if (toricmode == 'qldci1') then
         toricmode = 'qldci_CQL3D'
         iH_dql = 1
         ispec_dql = minspec
         ascii_name = 'du0u0_input_1'
      elseif (toricmode=='qldci2')then
         toricmode = 'qldci_CQL3D'
         iH_dql = 2
         ispec_dql = pairspec
         ascii_name = 'du0u0_input_2'
      elseif (toricmode=='qldci')then
         toricmode = 'qldci_CQL3D'
         iH_dql = 1
         ispec_dql = minspec
         ascii_name = 'du0u0_input'
      endif

      !dump modified namelists
      !logic for simultaneous execution
      if ((arg_nphi_indx.gt.0).and.(toricmode.eq.'toric'))then
         WRITE(outfile,"(A,I2.2)") "torica.inp_",i
      else
         outfile = "torica.inp"
      endif
      
      open(unit=out_unit, file='torica.inp', status='unknown', &
              form='formatted')
      write(out_unit, nml = toric_mode)
      write(out_unit, nml = qldciinp) 
      write(out_unit, nml = toricainp)
      write(out_unit, nml = equidata)
      write(out_unit, nml = ips)
      close(out_unit)
      
      ! Write an equidt.data/profnt file with plasma state profiles
      ! -----------------------------------------------------------------------------
      
      nprodt = size(ps%rho)
      nproeq = size(ps%rho_eq)

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
      allocate( x_toric (nprodt))

      open(unit=out_unit, file=profnt_file,              &
         status = 'unknown', form = 'formatted')
      write(out_unit,'(A10,5i4)') 'PS-state', nprodt, nspec, mainsp, &
                                  kdiff_idens, kdiff_itemp
      do  isp=1,nspec
         write(out_unit,'(2i4)')  int(atm(isp)), int(azi(isp))
      end do

      !set up zone centered grid
      do irho = 1,nprodt-1
        x_orig(irho) = 0.5 * (ps%rho(irho) + ps%rho(irho+1))
      end do

      !map zone centered tor flux to zone centered pol flux
      psi_poloidal_eq = ps%psipol
      call ps_user_1dintrp_vec(ps%rho, ps%rho_eq, psi_poloidal_eq, x_toric, ierr )
      x_toric = sqrt(x_toric/x_toric(nprodt))

      !output zone centered pol flux grid
      write(out_unit,'(A20)')  'Sqrt(Psi / Psilim)'
      write(out_unit,'(5E16.9)') x_toric

      ! Interpolate the electron density profile from the Plasma State grid to the Toric grid
      call ps_user_1dintrp_vec(ps%rho, x_orig, ps%ns(:,0), tmp_prof(:),ierr )  !DBB 6-27_2017
      if(ierr .ne. 0) stop 'error interpolating PS electron density profile onto Toric grid'

      call ps_user_1dintrp_vec(ps%rho, ps%rho_eq, ps%vol(:), vol_int(:),ierr ) !DBB 6-27_2017
      if(ierr .ne. 0) stop 'error interpolating PS volume onto Toric grid'

      write(out_unit,'(A10)')  'n_e'
      write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3
      
      call ps_user_1dintrp_vec(ps%rho, x_orig, ps%Ts(:,0), tmp_prof(:),ierr ) !DBB 6-27_2017
      if(ierr .ne. 0) stop 'error interpolating PS electron temperature profile onto Torlh grid'

      write(out_unit,'(A10)')  't_e'
      write(out_unit,'(5E16.9)')  tmp_prof   !keV

      call ps_tha_fetch (ierr, ps, tol_zero, ns_tha, ts_tha, v_pars, iwarn)
      if(ierr .ne. 0) stop 'error getting the density and temperature data from the abridged species list'

      do isp=1,ps%nspec_tha
         call ps_user_1dintrp_vec(ps%rho, x_orig, ns_tha(:,isp), tmp_prof(:),ierr ) !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS ion density profile onto Torlh grid'

         write(out_unit,'(A4,I2.2)')  'n_i_',isp
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3
	
         call ps_user_1dintrp_vec(ps%rho, x_orig, Ts_tha(:,isp), tmp_prof(:),ierr )  !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS ion temperature profile onto Torlh grid'

         write(out_unit,'(A4,I2.2)')  't_i_',isp
         write(out_unit,'(5E16.9)')  tmp_prof !keV
      end do

      
      ! Move on to the non-thermal ion species - in order of NBI, RF minority, and fusion alphas

      ! SF USER BEWARE: NBI NOT VALIDATED
      !
      ! First interpolate the NBI density [nbeami(:,1)] from the NuBeam grid onto the Torlh input data
      ! grid. Then write it to the equidt.data file.
      !
      ! Next interpolate the NBI energies [eperp_beami(:,1) and epll_beami(:,1)] from the NuBeam grid
      ! onto the Torlh input data grid. Then compute the equivalent temperature profile for the NBI and
      ! write it to the equidt.data file.
      if(allocated(ps%nbeami)) then
       isp = ps%snbi_to_alla(1)
         write(*,*) "Fast ion name, A, Z = ", trim(ps%alla_name(isp)), &
                            NINT(ps%m_alla(isp)/ps_mp), NINT(ps%q_alla(isp)/ps_xe)
         
         call ps_user_1dintrp_vec(ps%rho_nbi, x_orig, ps%nbeami(:,1), tmp_prof(:),ierr )  !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS NBI density profile onto Torlh grid'

         write(out_unit,'(A4,I2.2)')  'n_i_',isp
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3

      
          call ps_user_1dintrp_vec(ps%rho_nbi, x_orig, ps%eperp_beami(:,1), &
               aa_prof(:),ierr )  !DBB 6-27_2017
          if(ierr .ne. 0) stop 'error interpolating PS NBI perp. energy profile onto Torlh grid'
          call ps_user_1dintrp_vec(ps%rho_nbi, x_orig, ps%epll_beami(:,1),  &
               bb_prof(:),ierr )  !DBB 6-27_2017
          if(ierr .ne. 0) stop 'error interpolating PS NBI parallel energy profile onto Torlh grid'
          tmp_prof = 0.667 * (aa_prof + bb_prof)

          write(out_unit,'(A4,I2.2)')  't_i_',isp
          write(out_unit,'(5E16.9)')  tmp_prof !keV
 	  
      endif

      
      ! Write the rf minority tail data to the Torlh equidt.data file:
      if(allocated(ps%nmini)) then
         isp = ps%rfmin_to_alla(1)
         !if minority density set using a fraction of the electrons
         if (trim(ps%kdens_rfmin) .EQ. 'fraction') then        
            call ps_user_1dintrp_vec(x_orig, ps%rho, ps%ns(:,0), tmp_prof(:),ierr ) 
            if(ierr .ne. 0) stop 'error interpolating PS minority density profile onto Torlh grid'
            tmp_prof(:) = ps%fracmin(1)*tmp_prof(:)

            write(out_unit,'(A4,I2.2)')  'n_rfmin_',isp
            write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3

            !put modified minority density profile back in the plasma state
            call ps_user_1dintrp_vec(x_orig,ps%rho_icrf,tmp_prof, ps%nmini(:,1),ierr ) 
            if(ierr .ne. 0) stop 'error interpolating new minority desnity profile onto PS grid'
         endif

         !if minority density set using plasma state profile 
         if (trim(ps%kdens_rfmin) .EQ. 'data') then
            call ps_user_1dintrp_vec(ps%rho_icrf,x_orig, ps%nmini(:,1), tmp_prof(:),ierr ) 
            if(ierr .ne. 0) stop 'error interpolating PS minority density profile onto Toric grid'

            write(out_unit,'(A4,I2.2)')  'n_rfmin_',isp
            write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3  	     
         endif
    
         ! If isThermal=1 then assume the RF minority tail temperature is equal to the bulk ions
         if (ps%isThermal(1) .eq. 1) then
            call ps_user_1dintrp_vec(ps%rho, x_orig, Ts_tha(:,1), tmp_prof(:),ierr )  
            if(ierr .ne. 0) stop 'error interpolating PS ion temperature profile onto Torlh grid'
	
            write(out_unit,'(A4,I2.2)')  't_rfmin_',isp
            write(out_unit,'(5E16.9)')  tmp_prof !keV
         ! If isThermal=0 then the RF minority tail temperature data in the PS
         else
            call ps_user_1dintrp_vec(ps%rho_icrf, x_orig, ps%eperp_mini(:,1), aa_prof(:),ierr ) 
            if(ierr .ne. 0) stop 'error interpolating PS RF minority perp. energy profile onto Torlh grid'
            call ps_user_1dintrp_vec(ps%rho_icrf, x_orig, ps%epll_mini(:,1), bb_prof(:),ierr )  
            if(ierr .ne. 0) stop 'error interpolating PS RF minority parallel energy profile onto Torlh grid'
            tmp_prof = 0.667 * (aa_prof + bb_prof)

            write(out_unit,'(A4,I2.2)')  't_rfmin_',isp
            write(out_unit,'(5E16.9)')  tmp_prof !keV
	 endif
      endif

      ! Write the fast fusion alpha data to the Torlh equidt.data file:
      if(allocated(ps%nfusi)) then
         isp = ps%sfus_to_alla(1)

         !interpolate alpha density profile
         call ps_user_1dintrp_vec(ps%rho_fus,x_orig, ps%nfusi(:,1),tmp_prof(:),ierr )  
         if(ierr .ne. 0) stop 'error interpolating PS fusion alpha density profile onto Torlh grid'
         write(out_unit,'(A4,I2.2)')  'n_fus_',isp
         write(out_unit,'(5E16.9)')  tmp_prof*cubic_cm !M^-3 to cm^-3

         !interpolate alpha temp profile from energies
         call ps_user_1dintrp_vec(ps%rho_fus,x_orig, ps%eperp_fusi(:,1),aa_prof(:),ierr )  !DBB 6-27_2017
         if(ierr .ne. 0) stop 'error interpolating PS fusion alpha perp. energy profile onto Torlh grid'

         call ps_user_1dintrp_vec(ps%rho_fus,x_orig, ps%epll_fusi(:,1),bb_prof(:),ierr )
         if(ierr .ne. 0) stop 'error interpolating PS fusion alpha parallel energy profile onto Torlh grid'

         tmp_prof = 0.667 * (aa_prof + bb_prof)
         write(out_unit,'(A4,I2.2)')  't_i_',isp
         write(out_unit,'(5E16.9)')  tmp_prof !keV
      endif

      ! Update Plasma State with new information - ps%nmini, ps%eperp_mini, ps%epll_mini
      CALL ps_store_plasma_state(ierr , trim(cur_state_file))
      if (ierr .ne. 0) stop 'cannot open state in perpare torlh input'

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

      end program prepare_toric_input
