PROGRAM model_EPA_mdescr

! version 1.0 4/10/2016 (Batchelor)  

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
!       The code requires 5 command-line arguments
!       1) path to the current plasma state file
!       2) path to the current plasma eqdsk file
!       3) path to current jsdsk file - cur_jsdsk_file
!       4) action mode, i.e. one of: "INIT", "STEP", or "FINALIZE"
!       5) time stamp the time set by the driver component to which the simulation is 
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
!						   data used to initialize the EPA part of the plasma state (e.g. a 
!						   plasma state file to copy)
!
!   /evolving_model_data/ -> Contains data that defines the ad hoc models used to give time 
!                            varying numbers for the EPA data (e.g names of different models
!							 and parameters to be used in the models, like profiles shapes or
!							 parameters to define time variation
!
!
!    Note that there is a huge amount of EPA data which the model_EPA_mdescr models do not use.  
!    Therefore only a subset of the EPA data is initialized and modeled. 
!    Typically this is only what is used by other components or that is watched by the 
!    monitor component.  More could be added at any time.

! Data owned by EPA that is initialized here:
!	
! 	MACHINE_DESCRIPTION data:  none
!	SHOT_CONFIGURATION DATA: none
!	SIMULATION_INIT: 
!		nrho			! dimension of rho grid
!		rho				! rho grid
!	STATE_DATA:
!		vsur			! toroidal voltage at surface
!	STATE_PROFILES:
!		curr_bootstrap	! Neoclassical bootstrap current
!		curr_ohmic		! Ohmic current
!		pohme			! Ohmic heating of thermal electrons
!		ns				! thermal species density ns(nrho-1, 0:nspc_th)
!		Ts				! thermal species temperature ns(nrho-1, 0:nspc_th)
!		Zeff			! Total Z effective
!		
!
! Data owned by IC that is initialized here:
!	
!	STATE_DATA:
!		power_ic		! power on each icrf source
!	STATE_PROFILES:
!		nmini			! minority density profiles (maybe - note: IC sets nrho_icrf and rho_icrf)
!		
! Data owned by NBI that is initialized here:
!	STATE_DATA:
!		power_nbi		! power on each beam
!		
!--------------------------------------------------------------------------


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

    CHARACTER (len=256) :: cur_state_file
    CHARACTER(len=32) :: mode
    CHARACTER(len=32) :: time_stamp
    INTEGER :: nrho

	!--------------------------------------------------------------------------
	!   Internal data
	!--------------------------------------------------------------------------
    INTEGER :: nzone
    REAL(KIND=rspec), ALLOCATABLE :: zone_center(:)


    CHARACTER(len=32) :: EPA_profile_model_name
    
    ! namelist parameters for Lorentz_Linear model:
    !   rho_max = rho of peak of the Lorentzian (not exactly the peak of the profile)
    !   w = width of Lorentzian
    !   f0 = value of normalized profile on axis, rho = 0
    !   f1 = value of normalized profile at rho = 1
    REAL(KIND=rspec) :: nbeam_peak, rho_max_nbeami, w_nbeami, f0_nbeami, f1_nbeami
    REAL(KIND=rspec) :: FP_th_e_beam,rho_max_P_th_e, w_P_th_e, f0_P_th_e, f1_P_th_e
    REAL(KIND=rspec) :: FP_th_i_beam, rho_max_P_th_i, w_P_th_i, f0_P_th_i, f1_P_th_i
    REAL(KIND=rspec) :: I_beam_MA, rho_max_I_beam, w_I_beam, f0_I_beam, f1_I_beam
    REAL(KIND=rspec) :: Epll_peak, rho_max_Epll, w_Epll, f0_Epll, f1_Epll
    REAL(KIND=rspec) :: Eperp_peak, rho_max_Eperp, w_Eperp, f0_Eperp, f1_Eperp
	
	!------------------------------------------------------------------------------------
	!   Model input namelists
	!------------------------------------------------------------------------------------

    namelist /state_data/ &
          mode, cur_state_file, cur_eqdsk_file, time_stamp, nrho
                       
    namelist /evolving_model_data/ &
          EPA_profile_model_name, &
          nbeam_peak, rho_max_nbeami, w_nbeami, f0_nbeami, f1_nbeami, &
          FP_th_e_beam, rho_max_P_th_e, w_P_th_e, f0_P_th_e, f1_P_th_e, &
          FP_th_i_beam, rho_max_P_th_i, w_P_th_i, f0_P_th_i, f1_P_th_i, &
          I_beam_MA, rho_max_I_beam, w_I_beam, f0_I_beam, f1_I_beam, &
          Epll_peak, rho_max_Epll, w_Epll, f0_Epll, f1_Epll, &
          Eperp_peak, rho_max_Eperp, w_Eperp, f0_Eperp, f1_Eperp
    
    
!------------------------------------------------------------------------------------
!     
!  Section 0: Setup
!
!------------------------------------------------------------------------------------
    WRITE (*,*)
    WRITE (*,*) 'generic_ps_init'      

!---------------------------------------------------------------------------------
!     
!  Get init configuration data from ps_init_nml_file
!
!---------------------------------------------------------------------------------

    OPEN (unit=21, file = 'model_EPA_mdescr_input.nml', status = 'old',   &
         form = 'formatted', iostat = ierr)
    IF (ierr .ne. 0) THEN
		CALL SWIM_error ('open', 'model_EPA_mdescr.f90','model_EPA_mdescr_input.nml')
		WRITE (*,*) 'model_EPA_mdescr.f90: Cannot open ', 'model_EPA_mdescr_input.nml'
		call exit(1)
	END IF

	read(21, nml = state_data)
	WRITE (*, nml = state_data)
      
	WRITE (*,*)
	WRITE (*,*) 'model_EPA_mdescr'      
	print*, 'cur_state_file = ', trim(cur_state_file)
	print*, 'mode = ', trim(mode)
	print*, 'time_stamp = ', trim(time_stamp)
      
	!------------------------------------------------------------------------------------
	!  Get current plasma state 
	!------------------------------------------------------------------------------------
			
    call ps_get_plasma_state(ierr, trim(cur_state_file))
    if(ierr .ne. 0) then
       print*, 'model_EPA_mdescr:failed to get_plasma_state'
       stop 1
    end if


	!--------------------------------------------------------------------------
	!    Open input namelist file
	!--------------------------------------------------------------------------

	read(21, nml = evolving_model_data)
	CLOSE (21)
	IF TRIM(mode) == 'INIT' then
		WRITE (*, nml = state_data)
    
!------------------------------------------------------------------------------------
!     
! INIT function
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'INIT') THEN
              
    WRITE(*,*) 'model_EPA_mdescr: INIT'          
    
Do the work

    !-------------------------------------------------------------------------- 
    ! Store initial plasma state
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)
        IF (ierr .ne. 0) THEN
			WRITE (*,*) 'model_EPA_mdescr INIT: PS_STORE_PLASMA_STATE failed'
			CALL EXIT(1)
		ELSE
        	WRITE (*,*) "model_EPA_mdescr: Stored initial Plasma State"    
        END IF
                    
END IF  ! End INIT function

!------------------------------------------------------------------------------------
!     
!  Section 2: STEP function - Change state data and store plasma state
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'STEP') THEN

    WRITE(*,*)
    WRITE(*,*) 'model_EPA_mdescr: STEP'


    !-------------------------------------------------------------------------- 
    !  2.2  model = const.  Don't touch the EPA data.  Pass plasma state through.
    !--------------------------------------------------------------------------

	IF (TRIM(EPA_profile_model_name) == 'const') THEN
	
		CONTINUE
		
    !-------------------------------------------------------------------------- 
    !  2.3  model = zeros.  Set all outputs to zero
    ! --------------------------------------------------------------------------

	ELSE IF (TRIM(EPA_profile_model_name) == 'zeros') THEN

		ps%nbeami = 0.
		ps%pbe = 0.        
		ps%pbi = 0.
		ps%pbth = 0.        
		ps%curbeam = 0.
		ps%epll_beami = 0.
		ps%eperp_beami = 0.
	
		
    !-------------------------------------------------------------------------- 
    !  2.4  model = profiles.  Actually generate the profiles from the models
    !--------------------------------------------------------------------------

	ELSE IF (TRIM(EPA_profile_model_name) == 'profiles') THEN
         
        nzone = ps%nrho_nbi -1
        ALLOCATE( zone_center(nzone), stat=istat )
        IF (istat /= 0 ) THEN
            CALL SWIM_error ('allocation', 'model_EPA_mdescr' , 'zone_nbi_center')
            ierr = istat
            STOP
        END IF  

        zone_center = ( ps%rho_nbi(1:nrho_nbi-1) + ps%rho_nbi(2:nrho_nbi) )/2.
           write(*,*) 'zone_center = ', zone_center
 
		! direct nbi power to thermal electrons
		CALL Lorentz_Linear(rho_max_P_th_e, w_P_th_e, f0_P_th_e, f1_P_th_e, &
							zone_center, ps%pbe)
		ps%pbe = SUM(ps%power_EPAI) * FP_th_e_beam * ps%pbe/SUM(ps%pbe)
		
		! direct nbi power to thermal ions
		CALL Lorentz_Linear(rho_max_P_th_i, w_P_th_i, f0_P_th_i, f1_P_th_i, &
							zone_center, ps%pbi)
		ps%pbi = SUM(ps%power_EPAI) * FP_th_i_beam * ps%pbi/SUM(ps%pbi)
		
		! beam current profile
		CALL Lorentz_Linear(rho_max_I_beam, w_I_beam, f0_I_beam, f1_I_beam, &
							zone_center, ps%curbeam)
		ps%curbeam = 1.0e6_rspec * I_beam_MA * ps%curbeam/SUM(ps%curbeam)
		
		DO i = 1, ps%nspec_beam

	 		! beam species densities
			CALL Lorentz_Linear(rho_max_nbeami, w_nbeami, f0_nbeami, f1_nbeami, &
							    zone_center, ps%nbeami(:,i))
			ps%nbeami(:,i) =  nbeam_peak * ps%nbeami(:,i)/MAXVAL(ps%nbeami(:,i))
		
            !  changed multipler to make EPA come out in keV	
	 		! beam Epll profile
			CALL Lorentz_Linear(rho_max_Epll, w_Epll, f0_Epll, f1_Epll, &
							    zone_center, ps%epll_beami(:,i))
			ps%epll_beami(:,i) = Epll_peak*ps%epll_beami(:,i)/MAXVAL(ps%epll_beami(:,i))

			! beam Eperp profile
			CALL Lorentz_Linear(rho_max_Eperp, w_Eperp, f0_Eperp, f1_Eperp, &
							    zone_center, ps%eperp_beami(:,i))
			ps%eperp_beami(:,i) = Eperp_peak*ps%eperp_beami(:,i)/ &
									MAXVAL(ps%eperp_beami(:,i))
				
		END DO

    END IF ! End of cases of different profile models
    
    !--------------------------------------------------------------------------    !
    ! Store the data in plasma_state file
    !--------------------------------------------------------------------------

	CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)

	WRITE (*,*) "model_EPA_mdescr: Stored Plasma State"    

	CALL sleep(3)   ! Wait for elvis to pick up
	
END IF ! End of STEP function

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
!	Model 1: Power_Parabolic_Offset(alpha, h, rho, f) -> Parabolic to a power
!	plus an offset. 
!		f = A (1 - rho**2)**alpha + B where A and B are chosen so that
!			Integral rho*f(rho) from 0 to 1 is unity, and f(1)/f(0) = h
!
!	Model 2: Lorentz_Linear(rho_max, w, f0, f1, rho, f)
!		f = Lorentzian(centered at rho = rho_max, width = w) + f0 + (f1 - f0)*rho
!
!	Model 3: Lorentz_Linear_norm(rho_max, w, f0, f1, rho, f)
!		f = A*Lorentzian(centered at rho = rho_max, width = w) + B +C*rho
!			where A, B, C are chosen to give Integral(rho*f(rho) from 0 to 1) = 1
!			f(0) = f0, and f(1) = f1
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
!   Power_Parabolic_Offset
!
!--------------------------------------------------------------------------
	
	SUBROUTINE  Power_Parabolic_Offset(alpha, h, rho, f)
	
	  !  quadratic profile generator
	
	  REAL(KIND=rspec), intent(in) :: alpha, h  ! exponent and edge to peak ratio
	  REAL(KIND=rspec), dimension(:), intent(in) :: rho  ! normalized sqrt tor. flux
	  REAL(KIND=rspec), dimension(:), intent(out) :: f  ! output profile
	
	  f = (2*(1 + alpha)*(h + (1 - h)*(1 - rho**2)**alpha))/(1 + h*alpha)
	
	END SUBROUTINE Power_Parabolic_Offset
	
	
	SUBROUTINE Lorentz_Linear(rho_max, w, f0, f1, rho, f)
	
	  !  quick & dirty profile generator
	
	  REAL(KIND=rspec), intent(in) :: rho_max, w  ! Peak location and width of Lorentzian
	  REAL(KIND=rspec), intent(in) :: f0,f1  ! axis and edge values
	  REAL(KIND=rspec), dimension(:), intent(in) :: rho  ! normalized sqrt tor. flux
	  REAL(KIND=rspec), dimension(:), intent(out) :: f  ! output profile
	  REAL(KIND=rspec), dimension(size(rho)) :: lorentz  ! Lorentzian part
	
	  lorentz = w**2/(w**2 + (rho - rho_max)**2)
	  		 
	  f = lorentz + f0 + (f1 - f0)*rho
	
	END SUBROUTINE Lorentz_Linear

	
	SUBROUTINE Lorentz_Linear_norm(rho_max, w, f0, f1, rho, f)
	
	  !  quick & dirty profile generator
	
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


END PROGRAM model_EPA_mdescr
