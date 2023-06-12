! SF cleanup 2023

PROGRAM model_EPA_mdescr

  USE plasma_state_mod
    
  USE swim_global_data_mod, only : &
       & rspec, ispec, &               ! int: kind specification for real and integer
       & SWIM_name, SWIM_filename, &   ! derived data types: containing one character
                                                                        ! string
       & SWIM_error                    ! subroutine: a simple error handling routine
    
  IMPLICIT none

  !Command line args
  CHARACTER (len=256) :: cur_state_file
  CHARACTER(len=32) :: mode
  CHARACTER(len=32) :: time_stamp

  !Internal data
  INTEGER, PARAMETER :: maxDim = 10 ! To avoid a lot of allocates
  INTEGER :: nzone
  REAL(KIND=rspec), ALLOCATABLE :: zone_center(:)

  !State data
  INTEGER :: nrho
  INTEGER :: isThermal
  REAL(KIND=rspec) :: fracmin, power_ic, power_lh, power_ec, power_nbi
  CHARACTER*32 kdens_rfmin

  !Evolving model data
  CHARACTER(len=32) :: Te_profile_model_name = 'Power_Parabolic' 
  CHARACTER(len=32) :: ne_profile_model_name = 'Power_Parabolic'
  CHARACTER(len=32) :: Zeff_profile_model_name = 'Power_Parabolic'
  CHARACTER(len=32) :: V_loop_profile_model_name = 'Power_Parabolic'
  CHARACTER(len=32) :: Ti_profile_model_name = 'fraction_of_electron'
  CHARACTER(len=32) :: ni_profile_model_name = 'fraction_of_electron'
  CHARACTER(len=32) :: T_min_profile_model_name = 'fraction_of_electron'
  CHARACTER(len=32) :: n_min_profile_model_name = 'fraction_of_electron'
    
  REAL(KIND=rspec) :: Te_0, Te_edge, alpha_Te_1, alpha_Te_2
  REAL(KIND=rspec) :: ne_0, ne_edge, alpha_ne_1, alpha_ne_2
  REAL(KIND=rspec) :: Ti_0, Ti_edge, alpha_Ti_1, alpha_Ti_2
  REAL(KIND=rspec) :: ni_0, ni_edge, alpha_ni_1, alpha_ni_2
  REAL(KIND=rspec) :: Zeff_0, Zeff_edge, alpha_Zeff_1, alpha_Zeff_2
  REAL(KIND=rspec) :: V_loop_0, V_loop_edge, alpha_V_loop_1, alpha_V_loop_2

  REAL(KIND=rspec) :: frac_ni(maxDim) = 0.0, frac_Ti(maxDim)
  REAL(KIND=rspec) :: fracmin_T, fracmin_n = 0.0

  REAL(KIND=rspec) :: Te_ratio, alpha_Te
  REAL(KIND=rspec) :: ne_ratio, alpha_ne
  REAL(KIND=rspec) :: T_min_0, T_min_ratio, alpha_Tmin
    
  REAL(KIND=rspec), ALLOCATABLE :: numerical_rho, numerical_ne, numerical_Te, numerical_Ti
  CHARACTER(len = 80) :: equidt_file_name

  !Input namelists
  namelist /state_data/ &
          nrho, kdens_rfmin, isThermal, fracmin, power_ic, power_lh, power_ec, power_nbi
                       
  namelist /evolving_model_data/ &
          Te_profile_model_name, ne_profile_model_name, &
          Te_0, Te_edge, alpha_Te_1, alpha_Te_2, &
          ne_0, ne_edge, alpha_ne_1, alpha_ne_2, &
          Ti_profile_model_name, ni_profile_model_name, &
          Ti_0, Ti_edge, alpha_Ti_1, alpha_Ti_2, &
          ni_0, ni_edge, alpha_ni_1, alpha_ni_2, &
          frac_ni, frac_Ti, &
          fracmin_T, fracmin_n, &
          Te_ratio, alpha_Te, ne_ratio, alpha_ne, &
          T_min_0, T_min_ratio, alpha_Tmin, &
          Zeff_profile_model_name, V_loop_profile_model_name, &
          Zeff_0, Zeff_edge, alpha_Zeff_1, alpha_Zeff_2, &
          V_loop_0, V_loop_edge, alpha_V_loop_1, alpha_V_loop_2, &
          equidt_file_name

  !Setup
  !--------------------------------------------------------------------
  WRITE (*,*) 'model_EPA_mdescr'

  !get command line args
  call get_arg_count(iarg)
  if(iarg .ne. 3) then
     print*, 'model_EPA usage: '
     print*, ' command line args = cur_state_file mode time_stamp'
     stop 'incorrect command line arguments'
   end if
      
   call getarg(1,cur_state_file)
   call getarg(2,mode)
   call getarg(3,time_stamp)
    
   WRITE (*,*)
   print*, 'mode = ', trim(mode)
   print*, 'cur_state_file = ', trim(cur_state_file)
   print*, 'time_stamp = ', trim(time_stamp)

   !read input namelists
   OPEN (unit=21, file = 'model_EPA_mdescr_input.nml', status = 'old',   &
        form = 'formatted', iostat = ierr)
   IF (ierr .ne. 0) THEN
       CALL SWIM_error ('open', 'model_EPA_mdescr.f90','model_EPA_mdescr_input.nml')
       WRITE (*,*) 'model_EPA_mdescr.f90: Cannot open ', 'model_EPA_mdescr_input.nml'
       call exit(1)
   END IF
   
   read(21, nml = state_data)
   IF (TRIM(mode) == 'INIT') THEN
      WRITE (*, nml = state_data)
      WRITE (*,*)
   END IF

   ead(21, nml = evolving_model_data)
   CLOSE (21)
   IF (TRIM(mode) == 'INIT') THEN
      WRITE (*, nml = evolving_model_data)
      WRITE (*,*)
   END IF

   If (TRIM(ne_profile_model_name) == 'read_equidt_file') THEN		
      CALL equidt_to_PS(TRIM(equidt_file_name), 'INIT')
   END IF

   !get the plasma state
   call ps_get_plasma_state(ierr, trim(cur_state_file))
   if(ierr .ne. 0) then
      print*, 'model_EPA_mdescr:failed to get_plasma_state'
      stop 1
   end if

   IF (trim(mode)=='INIT') THEN
      ps%nrho = nrho

      IF (ALLOCATED(ps%power_ec)) THEN   		       
         ps%power_ec = power_ec
         ps%nrho_ecrf = nrho
      END IF
	
      IF (ALLOCATED(ps%power_ic)) THEN   		       
         ps%power_ic = power_ic
         ps%nrho_icrf = nrho
      END IF

      IF (ALLOCATED(ps%power_lh)) THEN   		       
         ps%power_lh = power_lh
         ps%nrho_lhrf = nrho
      END IF

      IF (ALLOCATED(ps%m_RFMIN)) THEN   
         ps%kdens_rfmin = "fraction"
         ps%fracmin = fracmin
         ps%isThermal = 1
      END IF

      IF (ALLOCATE(ps%m_SFUS)) THEN
         ps%nrho_fus = ps%nrho
      ENDIF
      
      WRITE (*,*) 'model_EPA_mdescr: About to allocate thermal profile arrays'
      CALL    ps_alloc_plasma_state(ierr)
      WRITE (*,*) 'model_EPA_mdescr:  Thermal profile arrays allocated ierr= ',ierr

      !SF Switched this to init. Probably don't want to evolve radial grids
      CALL rho_grid(ps%nrho,ps%rho)
      IF (ALLOCATED(ps%power_ec)) THEN   		       
         ps%rho_ecrf = ps%rho
      END IF
      IF (ALLOCATED(ps%power_ic)) THEN   		       
         ps%rho_icrf = ps%rho
      END IF
      IF (ALLOCATED(ps%power_lh)) THEN   		       
         ps%rho_lhrf = ps%rho
      END IF
      IF (ALLOCATE(ps%m_SFUS)) THEN
         ps%nrho_fus = ps%nrho
      ENDIF
      
   ENDIF

   !Set profiles
   nzone = ps%nrho -1       
   ALLOCATE( zone_center(nzone), stat=istat )
   IF (istat /= 0 ) THEN
      CALL SWIM_error ('allocation', 'model_epa' , 'zone_center')
      ierr = istat
   END IF
   zone_center = ( ps%rho(1:nrho-1) + ps%rho(2:nrho) )/2.
   
   IF (TRIM(Te_profile_model_name) == 'Power_Parabolic') THEN
      CALL Power_Parabolic(Te_0, Te_edge, alpha_Te_1, alpha_Te_2, zone_center, ps%Ts(:, 0))
   END IF

   IF (TRIM(ne_profile_model_name) == 'Power_Parabolic') THEN
      CALL Power_Parabolic(ne_0, ne_edge, alpha_ne_1, alpha_ne_2, zone_center, ps%ns(:, 0))
   END IF

   IF (TRIM(Ti_profile_model_name) == 'Power_Parabolic') THEN
      DO i = 1, ps%nspec_th
         CALL Power_Parabolic(Ti_0, Ti_edge, alpha_Ti_1, alpha_Ti_2, zone_center, ps%Ts(:, i))
      END DO
   END IF

   IF (TRIM(Zeff_profile_model_name) == 'Power_Parabolic') THEN
      CALL Power_Parabolic(Zeff_0, Zeff_edge, alpha_Zeff_1, alpha_Zeff_2, zone_center, ps%Zeff(:))
   END IF

   IF (TRIM(V_loop_profile_model_name) == 'Power_Parabolic') THEN
      CALL Power_Parabolic(V_loop_0, V_loop_edge, alpha_V_loop_1, alpha_V_loop_2, ps%rho(:), ps%V_loop(:))
   END IF

   IF (TRIM(Ti_profile_model_name) == 'fraction_of_electron') THEN
      DO i = 1, ps%nspec_th
         ps%Ts(:,i) = frac_Ti(i)*ps%Ts(:, 0)
         !WRITE (*,*) 'model_EPA_mdescr:  initial Ti profile = ', ps%Ts(:, i)
         !WRITE (*,*)
      END DO
   END IF
   
   IF (TRIM(ni_profile_model_name) == 'fraction_of_electron') THEN
      ! SF 12/06/2022
      ! Thermal species AKA: not RF minority species, must be set respecting
      ! quasineutrality and accounting for minority in the namelist. This
      ! is probably not the best way to do things because it doesn't account
      ! for user error. However, other methods are also not robust to user
      ! error and confusion and are more complicated.
      DO i = 1, ps%nspec_th
         ps%ns(:,i) = frac_ni(i)*ps%ns(:, 0)
      END DO
   END IF  ! fraction_of_electron
   
   If (TRIM(ne_profile_model_name) == 'read_equidt_file') THEN
      CALL equidt_to_PS(TRIM(equidt_file_name), 'STEP')
   END IF

   ps%Ti = ps%Ts(:, 1)
   
   !update source powers if needed
   IF (ALLOCATED(ps%power_ec)) THEN   		       
      ps%power_ec = power_ec
   END IF
	
   IF (ALLOCATED(ps%power_ic)) THEN   		       
      ps%power_ic = power_ic
      ps%nrho_icrf = nrho
   END IF

   IF (ALLOCATED(ps%power_lh)) THEN   		       
      ps%power_lh = power_lh
      ps%nrho_lhrf = nrho
   END IF

   CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)
   IF (ierr .ne. 0) THEN
      WRITE (*,*) 'model_EPA_mdescr: PS_STORE_PLASMA_STATE failed'
      CALL EXIT(1)
   ELSE
      WRITE (*,*) "model_EPA_mdescr: Stored Plasma State"    
   END IF

   IF (TRIM(mode) .eq. 'FINALIZE') THEN                                           
      WRITE (*,*) "model_EPA_mdescr: Finalize called"
   ENDIF

CONTAINS

    !------------------------------------
    SUBROUTINE rho_grid(nrho,rho)

      ! generate evenly spaced rho grid

      integer, intent(in) :: nrho    ! # of pts covering [0:1]
      REAL(KIND=rspec) :: rho(nrho)  ! the rho grid being generated

      REAL(KIND=rspec), parameter :: ONE = 1.0_rspec

      integer :: irho

      do irho = 1,nrho
         rho(irho) = (irho-1)*ONE/(nrho-1)
      enddo

    END SUBROUTINE rho_grid

    SUBROUTINE  Power_Parabolic(f0, f_edge, exp1, exp2, rho, f)
    
      !  Generalized parabolic profile generator, dimensional (i.e. un-normalized)
    
      REAL(KIND=rspec), intent(in) :: f0, f_edge  ! Peak and edge values
      REAL(KIND=rspec), intent(in) :: exp1, exp2  ! rho exponent and parabolic exponent
      REAL(KIND=rspec), dimension(:), intent(in) :: rho  ! normalized sqrt tor. flux
      REAL(KIND=rspec), dimension(:), intent(out) :: f  ! output profile
    
      f = (f0 - f_edge)*(1. - rho**exp1)**exp2 + f_edge
    
    END SUBROUTINE Power_Parabolic

    
    SUBROUTINE  Power_Parabolic_Offset(alpha, h, rho, f)
    
      !  Generalized parabolic profile generator, normalized (i.e. integrates to 1.0)
    
      REAL(KIND=rspec), intent(in) :: alpha, h  ! exponent and edge to peak ratio
      REAL(KIND=rspec), dimension(:), intent(in) :: rho  ! normalized sqrt tor. flux
      REAL(KIND=rspec), dimension(:), intent(out) :: f  ! output profile
    
      f = (2*(1 + alpha)*(h + (1 - h)*(1 - rho**2)**alpha))/(1 + h*alpha)
    
    END SUBROUTINE Power_Parabolic_Offset

    
    SUBROUTINE Lorentz_Linear(rho_max, w, f0, f1, rho, f)
    
      !  Quick & dirty profile generator with off-axis peaking
    
      REAL(KIND=rspec), intent(in) :: rho_max, w  ! Peak location and width of Lorentzian
      REAL(KIND=rspec), intent(in) :: f0,f1  ! axis and edge values
      REAL(KIND=rspec), dimension(:), intent(in) :: rho  ! normalized sqrt tor. flux
      REAL(KIND=rspec), dimension(:), intent(out) :: f  ! output profile
      REAL(KIND=rspec), dimension(size(rho)) :: lorentz  ! Lorentzian part
    
      lorentz = w**2/(w**2 + (rho - rho_max)**2)
             
      f = lorentz + f0 + (f1 - f0)*rho
    
    END SUBROUTINE Lorentz_Linear

    
    SUBROUTINE Lorentz_Linear_norm(rho_max, w, f0, f1, rho, f)
    
      !  Quick & dirty profile generator with off-axis peaking, normalized (i.e. integrates to 1.0)
    
      REAL(KIND=rspec), intent(in) :: rho_max, w  ! Peak location and width of Lorentzian
      REAL(KIND=rspec), intent(in) :: f0,f1  ! axis and edge values
      REAL(KIND=rspec), dimension(:), intent(in) :: rho  ! normalized sqrt tor. flux
      REAL(KIND=rspec), dimension(:), intent(out) :: f  ! output profile
      REAL(KIND=rspec) :: L0, L1  ! axis and edge values of Lorentzian
      REAL(KIND=rspec), dimension(size(rho)) :: lorentz  ! Lorentzian part
      REAL(KIND=rspec) :: I0  ! integral of Lorentzian
      REAL(KIND=rspec) :: a, b, c  ! coefficients

      REAL(KIND=rspec), parameter :: ONE = 1.0_rspec
    
      lorentz = w**2/(w**2 + (rho - rho_max)**2)
      
      L0 = w ** 2/(w ** 2 + rho_max ** 2)
      L1 = w**2/(w**2 + (1 - rho_max)**2)
      I0 = (w*(2*rho_max*ATAN((1 - rho_max)/w) + 2*rho_max*ATAN(rho_max/w) +  &
           w*Log(1 + (1 - 2*rho_max)/(w**2 + rho_max**2))))/2.
      
      a = -((6 - f0 - 2*f1)/(-6*I0 + L0 + 2*L1))
      b = (2*(-3*L0 + 3*I0*f0 - L1*f0 + L0*f1))/(6*I0 - L0 - 2*L1)
      c = (-3*(-2*L0 + 2*L1 + 2*I0*f0 - L1*f0 - 2*I0*f1 + L0*f1))/ &
          (6*I0 - L0 - 2*L1)
         
      f = a*lorentz + b + c*rho
    
    END SUBROUTINE Lorentz_Linear_norm


      SUBROUTINE equidt_to_PS(equidt_file, mode)
!
! 11/20/2017 (DBB)
! Reads a reduced version of equidt.data file, four variables - psi(rho), ne(rho) = ni(rho),
! Te(rho), and Ti(rho).  Loads these into plasma state by host association.  Note, equidt
! data is node based, PS data is zone centered.
! 
! +--------------------------------------------------------------------+
!  Density and temperature profiles
! +--------------------------------------------------------------------+
!
      implicit none
      
      character(len=*) :: equidt_file, mode
      integer, parameter :: lun22 = 22
      character(80) :: var_name
      integer, parameter :: nspmx = 8
                             ! A maximum of 8 ion species allowed

      integer ::   idprof, nspec,  mainsp, nprodt,isp,nproeq,irho
      integer :: kdiff_idens, kdiff_itemp, nsptmp, iatm, iazi,ierr

      REAL(KIND=rspec), dimension(:), allocatable :: psipro
      REAL(KIND=rspec), dimension(:), allocatable :: tbne, tbte
      REAL(KIND=rspec), dimension(:), allocatable :: tbti
      REAL(KIND=rspec), dimension(:), allocatable :: tbni !JPL for ATOM project
      REAL(KIND=rspec), dimension(:), allocatable :: tbpsi
      REAL(KIND=rspec), dimension(:), allocatable :: tbphi
      REAL(KIND=rspec), dimension(:), allocatable :: tmpar
      real(rspec), dimension(:),allocatable :: psi_poloidal_eq,x_orig


! +--------------------------------------------------------------------+

         OPEN(lun22,file=TRIM(equidt_file),status='old',iostat=ierr)
         if (ierr/=0) then
            write(*,*) "Error opening file: ", equidt_file
            ierr=0
         endif
!
!  Reading the first variable name and the number of radial mesh points
!
         read(lun22,'(A10,5i4)')  var_name, nprodt, nspec, mainsp, &
                                  kdiff_idens, kdiff_itemp
         write (*,*) 'var_name = ', var_name, '  nprodt = ', nprodt
!
! On INIT just set nrho to nprodt and return (DBB)
!
		If (TRIM(mode) == 'INIT') THEN  
			nrho = nprodt
			close(lun22)
			return
		ENDIF
         do  isp=1,nspec
            read(lun22,'(2i4)')  iatm, iazi
         enddo

!
! Allocations
!
         nproeq = size(ps%rho_eq)
         write(*,*) 'before allocations nprodt nproeq = ',nprodt,nproeq
         allocate(tbpsi(nprodt))
         allocate(psipro(nprodt))
         allocate(tbphi(nprodt))
         allocate(tbne(nprodt),stat=ierr)
         allocate(tbte(nprodt),stat=ierr)
         allocate(tbti(nprodt),stat=ierr)
         allocate(tbni(nprodt),stat=ierr)
         allocate(tmpar(nprodt),stat=ierr)
         allocate( x_orig (nprodt-1))
         allocate( psi_poloidal_eq(nproeq))
         write(*,*) 'after allocations'

!  Reading the radial mesh
!  NOTE: psi is 0 at the magnetic axis and 1 at the plasma
!  edge, and is linear in SQRT(Psi_poloidal). An equidistant
!  mesh is required for the interpolation in toric.
!
         read(lun22,'(A10)')  var_name
         read(lun22, *)  tbpsi
         write (*,*) 'tbpsi = ', tbpsi

         if(tbpsi(nprodt) .ne. 1._rspec)  then
            write(*,*) "warning profnt radial mesh does not end at 1"
            tbpsi(1:nprodt) = tbpsi(1:nprodt)/tbpsi(nprodt)
         endif

         psi_poloidal_eq = sqrt(ps%psipol / ps%psipol(nproeq))
         do irho = 1,nprodt-1
            x_orig(irho) = 0.5 * (ps%rho(irho) + ps%rho(irho+1))
         end do
         WRITE(*,*)'psipoleq ',psi_poloidal_eq 
         WRITE(*,*)'rhoeq ',ps%rho_eq
         call ps_user_1dintrp_vec(tbpsi,psi_poloidal_eq,ps%rho_eq &
              ,tbphi,ierr)
         WRITE(*,*) 'tbphi ' ,tbphi
         
!
!  Reading the particle densities (hardwired ni = ne and one ion species)
!
         read(lun22,*)  var_name
         read(lun22, *)  tbne
         write (*,*) 'tbne = ', tbne
         call ps_user_1dintrp_vec(x_orig,tbphi,tbne,tmpar,ierr)
         !call ps_user_1dzone_centered_profile(nprodt, tmparr*1.0e6)
         ps%ns(:,0) = tmpar*1.0e6
!
!  Reading the electron temperature (units: keV)
!
         read(lun22,*)  var_name
         read(lun22,*)  tbte(1:nprodt)
         write (*,*) 'tbte = ', tbte
         call ps_user_1dintrp_vec(x_orig,tbphi,tbte,tmpar,ierr)
         ps%Ts(:,0) = tmpar
!
!  Reading the particle densities (hardwired ni = ne and one ion species)
!
         read(lun22,*)  var_name
         read(lun22, *)  tbni
         write (*,*) 'tbni = ', tbni
         call ps_user_1dintrp_vec(x_orig,tbphi,tbni,tmpar,ierr)
         ps%ns(:,1) = tmpar
!
!  Reading the ion temperature
!
         read(lun22,*)  var_name
	     read(lun22,*)  tbti(1:nprodt)
         write (*,*) 'tbti = ', tbti
         
         call ps_user_1dintrp_vec(x_orig,tbphi,tbti,tmpar,ierr)
         ps%Ts(:,1) = tmpar

         write(*,*) 'finished reading profiles'
         close(lun22)

      return
      end subroutine equidt_to_PS

! ----------------------------------------------------------------------------------------      
      function zone_centered_profile(n, profile)      
		  implicit none
	! Argument declarations
		  INTEGER, INTENT(in) :: n
		  REAL(KIND=rspec), INTENT(in), dimension(n) :: profile
		  REAL(KIND=rspec), dimension(n-1) :: zone_centered_profile
			
		  zone_centered_profile(:) = 0.5*(profile(1:n-1) + profile(2:n))
		  RETURN
	  END function zone_centered_profile
      

END PROGRAM model_EPA_mdescr

