PROGRAM model_EPA_mdescr

! Version 2.0 9/19/2018 (Batchelor)
! Added source_powers_only flag so that profiles are not initialized or set,
! only source powers are set. For now all sources are assumed to have the same power as
! set in the /state_data/ namelist.  With source_powers_only == true 
! nothing is initialized.  Use with fully allocated plasma state.
!
! version 1.1 11/20/2017 (Batchelor)
! Reads cut down version of TORIC equidt.data file for thermal profiles

!--------------------------------------------------------------------------
!
!   A code that writes equilibrium and thermal profile data into a plasma state
!   file based on simple analytic models.  This version requires that the necessary
!   machine description and shot configuration data has been written into the initial
!   plasma state file, as for example by the generic_ps_file_init component in the mdescr
!   mode.  This code allocates and initializes the plasma profiles from data in an
!   input namelist file, model_EPA_mdescr_input.nml.  This code is much less ambitious than its
!   predecessor, model_EPA_mdescr_2.f90, in that it only provides models for the thermal plasma
!   species density and temperature.  In particular it does nothing with the MHD equilibrium.
!
!   These models could become arbitrarily complex but for right now (4/2016) things are
!   very simple.  Ion density and temperature profiles are taken as fractions of the 
!   electron profile.  If you want more, call me (DBB).
!
!   NB: The only models implemented now (4/2016) for thermal ion and minority ion profiles
!       are: "fraction_of_electron" and "read_equidt_file"
!
!   NB: For now assume only one ICRF source and one minority species.  Easy to
!       generalize. Also for now the minority is assumed to be thermalized.
!
!   NB: For now assume only one LH source.  Easy to generalize.
!
!       The code requires 3 command-line arguments
!       1) path to the current plasma state file
!       2) action mode, i.e. one of: "INIT", "STEP", or "FINALIZE"
!       3) time stamp the time set by the driver component to which the simulation is 
!          supposed to advance.
!   
!
!       Don Batchelor
!       ORNL
!       Oak Ridge, TN 37831
!       batchelordb@ornl.gov
!
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! Details:
!   The data to be loaded is derived from input namelist file model_EPA_mdescr_input.nml.
!   This file must contain two namelists:
!
!   /static_state_data/ -> Contains plasma state data that goes directly into the state, or
!                          data used to initialize the EPA part of the plasma state (e.g. a 
!                          plasma state file to copy)
!
!   /evolving_model_data/ -> Contains data that defines the ad hoc models used to give 
!                            profiles and time variations for the EPA data (e.g names of
!                            different models and parameters to be used in the models, like profiles shapes or
!                            like profiles shapes or parameters to define time variation
!
!
! Data owned by EPA that is initialized here:
!   
!   MACHINE_DESCRIPTION data:  none
!   SHOT_CONFIGURATION DATA: none
!   SIMULATION_INIT: 
!       nrho            ! dimension of rho grid
!       rho             ! rho grid
!   STATE_DATA:
!       vsur            ! toroidal voltage at surface
!   STATE_PROFILES:
!       curr_bootstrap  ! Neoclassical bootstrap current
!       curr_ohmic      ! Ohmic current
!       pohme           ! Ohmic heating of thermal electrons
!       ns              ! thermal species density ns(nrho-1, 0:nspc_th)
!       Ts              ! thermal species temperature ns(nrho-1, 0:nspc_th)
!       Zeff            ! Total Z effective
!		V_loop			! Loop voltage : Volts
!       
!
! Data owned by IC that is initialized here:
!
!   SIMULATION_INIT: 
!       kdens_rfmin 
!   STATE_DATA:
!       power_ic        ! power on each icrf source
!   STATE_PROFILES:
!       nmini           ! minority density profiles (note: IC sets nrho_icrf and rho_icrf)
!
! Data owned by LH that is initialized here:
!
!   SIMULATION_INIT: None
!   STATE_DATA:
!       power_lh        ! power on each icrf source
!   STATE_PROFILES: None
!       
!--------------------------------------------------------------------------

! Working notes:

! 11/20/2017 DBB: Notes on "read_equidt_file" model for thermal profiles
! Added coding to specify profiles from numerical data.  This is tailored to get data from
! a cut down version of the TORIC/TORLH equidt.data file.
!
! Note assumptions for now (these are easily generalized): 
! 1) Number of data values = dimension of rho grid (i.e. grid based, not zone based as  
!     in PS parlance)
! 2) Data is on sqrt(toroidal flux) grid (i.e. no re-gridding here)
! 3) To get data in zone_center form for PS (ns, Ts) I take the average of the bounding 
!    grid points
! 4) There is only one thermal ion species
!
! To use this model must have ne_profile_model_name set to 'read_equidt_file' and the input
! equidt input file must be defined in model_EPA_mdescr_input.nml. Nothing need be 
! specified about the Te, ni or Ti models.  The input equidt file name is arbitrary. e.g.
!
! ne_profile_model_name = 'read_equidt_file'
! equidt_file_name = equidt_test.data
!
! Both of these files must be specified as input files in the simulation config file.

! 10-10-2016 DBB
! Added coding for Zeff and V_loop


    USE plasma_state_mod
    
    USE swim_global_data_mod, only : &
            & rspec, ispec, &               ! int: kind specification for real and integer
            & SWIM_name, SWIM_filename, &   ! derived data types: containing one character
                                                                        ! string
            & SWIM_error                    ! subroutine: a simple error handling routine
    
    IMPLICIT none
    
    INTEGER :: ierr = 0, iwarn = 0
    INTEGER :: istat, iarg
    INTEGER :: i

    !--------------------------------------------------------------------------
    !   Command line args
    !--------------------------------------------------------------------------

    CHARACTER (len=256) :: cur_state_file
    CHARACTER(len=32) :: mode
    CHARACTER(len=32) :: time_stamp

    !--------------------------------------------------------------------------
    !   Internal data
    !--------------------------------------------------------------------------
    INTEGER, PARAMETER :: maxDim = 10 ! To avoid a lot of allocates
    INTEGER :: nzone
    REAL(KIND=rspec), ALLOCATABLE :: zone_center(:)

    !--------------------------------------------------------------------------
    !   State data
    !--------------------------------------------------------------------------

    INTEGER :: nrho
    INTEGER :: isThermal
    REAL(KIND=rspec) :: fracmin, power_ic, power_lh, power_ec, power_nbi
    CHARACTER*32 kdens_rfmin

    !--------------------------------------------------------------------------
    !   Evolving model data
    !--------------------------------------------------------------------------

    CHARACTER(len=32) :: Te_profile_model_name = 'Power_Parabolic' 
    CHARACTER(len=32) :: ne_profile_model_name = 'Power_Parabolic'
    CHARACTER(len=32) :: Zeff_profile_model_name = 'Power_Parabolic'
    CHARACTER(len=32) :: V_loop_profile_model_name = 'Power_Parabolic'
    CHARACTER(len=32) :: Ti_profile_model_name = 'fraction_of_electron'
    CHARACTER(len=32) :: ni_profile_model_name = 'fraction_of_electron'
    CHARACTER(len=32) :: T_min_profile_model_name = 'fraction_of_electron'
    CHARACTER(len=32) :: n_min_profile_model_name = 'fraction_of_electron'
    CHARACTER(len=32) :: source_powers_only = 'false'
    
    ! namelist parameters for Power_Parabolic model:
    REAL(KIND=rspec) :: Te_0, Te_edge, alpha_Te_1, alpha_Te_2
    REAL(KIND=rspec) :: ne_0, ne_edge, alpha_ne_1, alpha_ne_2
    REAL(KIND=rspec) :: Ti_0, Ti_edge, alpha_Ti_1, alpha_Ti_2
    REAL(KIND=rspec) :: ni_0, ni_edge, alpha_ni_1, alpha_ni_2
    REAL(KIND=rspec) :: Zeff_0, Zeff_edge, alpha_Zeff_1, alpha_Zeff_2
    REAL(KIND=rspec) :: V_loop_0, V_loop_edge, alpha_V_loop_1, alpha_V_loop_2

    ! Fractional ion parameters relative to electron profiles
    ! Temperature fractions can be arbitrary but density fractions need to be consistent
    ! with charge neutrality
    REAL(KIND=rspec) :: frac_ni(maxDim) = 0.0, frac_Ti(maxDim)
    REAL(KIND=rspec) :: fracmin_T, fracmin_n = 0.0

    ! namelist parameters for Power_Parabolic_Offset model:
    REAL(KIND=rspec) :: Te_ratio, alpha_Te
    REAL(KIND=rspec) :: ne_ratio, alpha_ne
    REAL(KIND=rspec) :: T_min_0, T_min_ratio, alpha_Tmin
    
    ! namelist parameters for Numerical_Data model:
    REAL(KIND=rspec), ALLOCATABLE :: numerical_rho, numerical_ne, numerical_Te, numerical_Ti
    CHARACTER(len = 80) :: equidt_file_name

    !------------------------------------------------------------------------------------
    !   Model input namelists
    !------------------------------------------------------------------------------------

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
          equidt_file_name, source_powers_only
    
    
!------------------------------------------------------------------------------------
!     
!  Section 0: Setup
!
!------------------------------------------------------------------------------------
    WRITE (*,*)
    WRITE (*,*) 'model_EPA_mdescr'
         
    !------------------------------------------------------------------------------------
    !   Get command line arguments
    !------------------------------------------------------------------------------------

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
     
    !---------------------------------------------------------------------------------
    !  Get state data from model_EPA_mdescr_input.nml
    !---------------------------------------------------------------------------------

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

    !---------------------------------------------------------------------------------
    !  Get model data from model_EPA_mdescr_input.nml
    !---------------------------------------------------------------------------------

	read(21, nml = evolving_model_data)
	CLOSE (21)
	IF (TRIM(mode) == 'INIT') THEN
		WRITE (*, nml = evolving_model_data)
		WRITE (*,*)
	END IF
    IF (TRIM(source_powers_only)=='FALSE') THEN  ! lower case it
    	source_powers_only='false'
    END IF
    !---------------------------------------------------------------------------------
    !  If getting thermal profile data from equidta.dat file (signaled by 
    !  ne_profile_model_name = 'equidt_profile' in evolving_model_data) read equidt.data file
    !  N.B. nrho from equidt.data overrides nrho from state_data namelist
    !---------------------------------------------------------------------------------

        
	If (TRIM(ne_profile_model_name) == 'read_equidt_file') THEN		
		CALL equidt_to_PS(TRIM(equidt_file_name), 'INIT')
	END IF
      
    !------------------------------------------------------------------------------------
    !  Get current plasma state 
    !------------------------------------------------------------------------------------
            
    call ps_get_plasma_state(ierr, trim(cur_state_file))
    if(ierr .ne. 0) then
       print*, 'model_EPA_mdescr:failed to get_plasma_state'
       stop 1
    end if
    
    
!------------------------------------------------------------------------------------
!     
! INIT function
!
!------------------------------------------------------------------------------------

	IF ((TRIM(mode) == 'INIT').and.(TRIM(source_powers_only)=='false')) THEN
		!--------------------------------------------------------------------------
		!   Initialize and allocate species arrays and thermal profile grids
		!--------------------------------------------------------------------------

			ps%nrho = nrho
	   
			WRITE (*,*) 'model_EPA_mdescr: About to allocate thermal profile arrays'
			CALL    ps_alloc_plasma_state(ierr)
			WRITE (*,*) 'model_EPA_mdescr:  Thermal profile arrays allocated'
			WRITE (*,*)
	
		!---------------------------------------------------------------------------------
		! ICRF minority ion profiles
		!---------------------------------------------------------------------------------
			! If (kdens_rfmin .EQ. 'fraction') TORIC computes nmini = fracmin * ne
			! If kdens_rfmin .EQ. 'data' (not implemented here as of 4/2016) then nmini must 
			! be available in the PS, TORIC interpolates it from the rho-icrf grid onto the 
			! Toric radial grid.  But the rho_icrf grid must be allocated and initialized here.  
			IF (ALLOCATED(ps%m_RFMIN)) THEN   
				ps%kdens_rfmin = "fraction"
				ps%fracmin(:) = fracmin
				ps%isThermal(:) = 1
				WRITE (*,*) 'model_EPA_mdescr: minority specification loaded'
				WRITE (*,*)
			END IF                    
	END IF  ! End INIT function

!------------------------------------------------------------------------------------
!     
!  STEP function - Change state data and store plasma state.  Time dependence is now
!                   implemented in the python
!
! N.B. This section is also executed during INIT
!------------------------------------------------------------------------------------
WRITE(*,*)

!---------------------------------------------------------------------------------
! Thermal species profiles
!---------------------------------------------------------------------------------
	IF (TRIM(source_powers_only) == 'false') THEN         
				
		CALL rho_grid(ps%nrho,ps%rho)
		nzone = ps%nrho -1       
		ALLOCATE( zone_center(nzone), stat=istat )
		IF (istat /= 0 ) THEN
			CALL SWIM_error ('allocation', 'model_epa' , 'zone_center')
			ierr = istat
		END IF  
		zone_center = ( ps%rho(1:nrho-1) + ps%rho(2:nrho) )/2.

		! PowerParabolic models

		IF (TRIM(Te_profile_model_name) == 'Power_Parabolic') THEN
			CALL Power_Parabolic(Te_0, Te_edge, alpha_Te_1, alpha_Te_2, zone_center, ps%Ts(:, 0))
			WRITE (*,*) 'model_EPA_mdescr:  initial Te profile = ', ps%Ts(:, 0)
			WRITE (*,*)
		END IF

		IF (TRIM(ne_profile_model_name) == 'Power_Parabolic') THEN
			CALL Power_Parabolic(ne_0, ne_edge, alpha_ne_1, alpha_ne_2, zone_center, ps%ns(:, 0))
			WRITE (*,*) 'model_EPA_mdescr:  initial ne profile = ', ps%ns(:, 0)
			WRITE (*,*)
		END IF

		! Thermal ion profiles  N.B.  All thermal ion species at same temperature
		IF (TRIM(Ti_profile_model_name) == 'Power_Parabolic') THEN
			DO i = 1, ps%nspec_th
				CALL Power_Parabolic(Ti_0, Ti_edge, alpha_Ti_1, alpha_Ti_2, zone_center, ps%Ts(:, i))
				WRITE (*,*) 'model_EPA_mdescr:  initial Ti profile = ', ps%Ts(:, i)
				WRITE (*,*)
			END DO
		END IF

		! Zeff profile
		IF (TRIM(Zeff_profile_model_name) == 'Power_Parabolic') THEN
			CALL Power_Parabolic(Zeff_0, Zeff_edge, alpha_Zeff_1, alpha_Zeff_2, zone_center, ps%Zeff(:))
			WRITE (*,*) 'model_EPA_mdescr:  Zeff profile = ', ps%Zeff(:)
			WRITE (*,*)
		END IF

		! V_loop profile
		IF (TRIM(V_loop_profile_model_name) == 'Power_Parabolic') THEN
			CALL Power_Parabolic(V_loop_0, V_loop_edge, alpha_V_loop_1, alpha_V_loop_2, zone_center, ps%V_loop(:))
			WRITE (*,*) 'model_EPA_mdescr:  V_loop profile = ', ps%V_loop(:)
			WRITE (*,*)
		END IF

		! Fraction of electron models

		IF (TRIM(Ti_profile_model_name) == 'fraction_of_electron') THEN
			DO i = 1, ps%nspec_th
				ps%Ts(:,i) = frac_Ti(i)*ps%Ts(:, 0)
				WRITE (*,*) 'model_EPA_mdescr:  initial Ti profile = ', ps%Ts(:, i)
				WRITE (*,*)
			END DO
		END IF

		IF (TRIM(ni_profile_model_name) == 'fraction_of_electron') THEN

			! NB: if minority ion density model is fraction of electron density  and if the 
			! thermal ion model is also fraction of electrons then adjust fraction of
			! thermal species 1 to give quasineutrality
			! Define fraction of species #1 = 1.0 then subtract off minorities and other
			! thermals.  (Also must subtract beams and fusion when they get put in)
			IF (ps%kdens_rfmin == "fraction") THEN
				frac_ni(1) = 1.0 
				DO i = 1, ps%nspec_rfmin
					frac_ni(1) = frac_ni(1) - ps%q_RFMIN(i)*ps%fracmin(i)
				END DO
				IF (ps%nspec_th .GE. 2) THEN
					DO i = 2, ps%nspec_rfmin
						frac_ni(1) = frac_ni(1) - ps%q_S(i)*frac_ni(i)
					END DO
				ENDIF
				IF (frac_ni(1) < 0.0) THEN 
					WRITE (*,*) 'model_EPA_mdescr INIT: frac_ni(1) < 0'
					CALL EXIT(1)
				ENDIF
		
			END IF ! kdens_rfmin "fraction"
	
			! NB: If minorities are not fraction_of_electron, the ion fractions are just those
			! in the evolving_model_data namelist.  Some other mechanism must enforce charge
			! neutrality.
			DO i = 1, ps%nspec_th
				ps%ns(:,i) = frac_ni(i)*ps%ns(:, 0)
				WRITE (*,*) 'model_EPA_mdescr:  initial density profile for thermal ion #',i ,' = ', ps%ns(:, i)
				WRITE (*,*)
			END DO
		
		END IF  ! fraction_of_electron

		! Load single Ti profile from multi-species Ts()
		ps%Ti = ps%Ts(:, 1)

		! read_equidt_file model

		If (TRIM(ne_profile_model_name) == 'read_equidt_file') THEN
			CALL equidt_to_PS(TRIM(equidt_file_name), 'STEP')
		END IF

		! Load single Ti profile from multi-species Ts()
		ps%Ti = ps%Ts(:, 1)
	END IF ! End thermal profiles
	
    !--------------------------------------------------------------------------    !
    ! Source powers
    !--------------------------------------------------------------------------
		
	! NB: For now assume power is the same on all EC sources
	IF (ALLOCATED(ps%power_ec)) THEN   		       
		ps%power_ec = power_ec
	END IF
	
	! NB: For now assume power is the same on all ICRF sources
	IF (ALLOCATED(ps%power_ic)) THEN   		       
		ps%power_ic = power_ic
	END IF

	! NB: For now assume power is the same on all LH sources
	IF (ALLOCATED(ps%power_lh)) THEN   		       
		ps%power_lh = power_lh
	END IF

	! NB: For now assume power is the same on all NBI sources
	IF (ALLOCATED(ps%power_nbi)) THEN   		       
		ps%power_ic = power_nbi
	END IF
		

    !-------------------------------------------------------------------------- 
    ! Store plasma state
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)
        IF (ierr .ne. 0) THEN
            WRITE (*,*) 'model_EPA_mdescr: PS_STORE_PLASMA_STATE failed'
            CALL EXIT(1)
        ELSE
            WRITE (*,*) "model_EPA_mdescr: Stored Plasma State"    
        END IF

!--------------------------------------------------------------------------
!
!   Finalize function
!
!--------------------------------------------------------------------------

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

    !------------------------------------
    SUBROUTINE th_grid(nth,th)

      ! generate evenly spaced theta grid [-pi:pi]

      integer, intent(in) :: nth    ! # of pts covering the range
      REAL(KIND=rspec) :: th(nth)   ! the theta grid being generated

      REAL(KIND=rspec), parameter :: PI = 3.1415926535897931_rspec

      integer :: ith

      do ith = 1,nth
         th(ith) = -PI + (ith-1)*2*PI/(nth-1)
      enddo

    END SUBROUTINE th_grid


!--------------------------------------------------------------------------
!
!   A set of simple models to generate radial profiles based on adjustable input
!   parameters.  Input rho is a vector of values assumed to run from 0 to 1.
!
!   Model 1: Power_Parabolic_Offset(alpha, h, rho, f) -> Parabolic to a power
!   plus an offset. 
!       f = A (1 - rho**2)**alpha + B where A and B are chosen so that
!
!   Model 2: Power_Parabolic_Offset(alpha, h, rho, f) -> Parabolic to a power
!   plus an offset. 
!       f = A (1 - rho**2)**alpha + B where A and B are chosen so that
!           Integral rho*f(rho) from 0 to 1 is unity, and f(1)/f(0) = h.  This is
!           useful for a source profile where the total power would be known.
!
!   Model 3: Lorentz_Linear(rho_max, w, f0, f1, rho, f)
!       f = Lorentzian(centered at rho = rho_max, width = w) + f0 + (f1 - f0)*rho
!
!   Model 4: Lorentz_Linear_norm(rho_max, w, f0, f1, rho, f)
!       f = A*Lorentzian(centered at rho = rho_max, width = w) + B +C*rho
!           where A, B, C are chosen to give Integral(rho*f(rho) from 0 to 1) = 1
!           f(0) = f0, and f(1) = f1
!
!
!       Don Batchelor
!       ORNL
!       Oak Ridge, TN 37831
!       batchelordb@ornl.gov
!
!--------------------------------------------------------------------------
  


!--------------------------------------------------------------------------
!
!   Profile functions
!
!--------------------------------------------------------------------------
    
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

      integer ::   idprof, nspec,  mainsp, nprodt
      integer :: kdiff_idens, kdiff_itemp, nsptmp, iatm, iazi

      REAL(KIND=rspec), dimension(:), allocatable :: psipro
      REAL(KIND=rspec), dimension(:), allocatable :: tbne, tbte
      REAL(KIND=rspec), dimension(:), allocatable :: tbti
      REAL(KIND=rspec), dimension(:), allocatable :: tbni !JPL for ATOM project
      REAL(KIND=rspec), dimension(:), allocatable :: tbpsi


! +--------------------------------------------------------------------+

         OPEN(lun22,file=TRIM(equidt_file),status='old',iostat=ierr)
         if (ierr/=0) then
            write(*,*) "Error opening file: ", equidt_file
            ierr=0
         endif
!
!  Reading the first variable name and the number of radial mesh points
!
!         read(lun22,'(A10,i4)')  var_name, nprodt
         read(lun22, *)  var_name, nprodt
         write (*,*) 'var_name = ', var_name, '  nprodt = ', nprodt
!
! On INIT just set nrho to nprodt and return (DBB)
!
		If (TRIM(mode) == 'INIT') THEN  
			nrho = nprodt
			close(lun22)
			return
		ENDIF
!
! Allocations
!
         write(*,*) 'before allocations nprodt = ',nprodt
         allocate(tbpsi(nprodt))
         allocate(psipro(nprodt))
         allocate(tbne(nprodt),stat=ierr)
         allocate(tbte(nprodt),stat=ierr)
         allocate(tbti(nprodt),stat=ierr)
         allocate(tbni(nprodt),stat=ierr)
         write(*,*) 'after allocations'

!  Reading the radial mesh
!  NOTE: psi is 0 at the magnetic axis and 1 at the plasma
!  edge, and is linear in SQRT(Psi_poloidal). An equidistant
!  mesh is required for the interpolation in toric.
!
         read(lun22, *)  tbpsi
         write (*,*) 'tbpsi = ', tbpsi

         if(tbpsi(nprodt) .ne. 1._rspec)  then
            write(*,*) "warning profnt radial mesh does not end at 1"
            tbpsi(1:nprodt) = tbpsi(1:nprodt)/tbpsi(nprodt)
         endif
         ps%rho = tbpsi(1:nprodt)         
!
!  Reading the particle densities (hardwired ni = ne and one ion species)
!
         read(lun22,*)  var_name, nprodt
         read(lun22, *)  tbne
         write (*,*) 'tbne = ', tbne
         ps%ns(:,0) = zone_centered_profile(nprodt, tbne)
         ps%ns(:,1) = ps%ns(:,0)
!
!  Reading the electron temperature (units: keV)
!
         read(lun22,*)  var_name, nprodt
         read(lun22,*)  tbte(1:nprodt)
         write (*,*) 'tbte = ', tbte
         ps%Ts(:,0) = zone_centered_profile(nprodt, tbte)
!
!  Reading the ion temperature
!
         read(lun22,*)  var_name, nprodt
	     read(lun22,*)  tbti(1:nprodt)
         write (*,*) 'tbti = ', tbti

         ps%Ts(:,1) = zone_centered_profile(nprodt, tbti)

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

