PROGRAM model_EPA

! version 2.0 8/26/2009 (Batchelor)  Added capability to initialize from an input plasma state
!             file. Added a model routine that gives off-axis peaked profiles: Lorentz_Linear.
!             Generally cleaned up.
! version 0.2 8/5/2009 (Batchelor)  Added beam species densities.
! version 0.0 3/13/2009 (Batchelor)

!--------------------------------------------------------------------------
!
!   Simple mock EPA code for testing.
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
!   The data to be loaded is derived from input namelist file model_EPA_input.nml.
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
! There are two ways to allocate the plasma state arrays for the EPA component:
!
! 1) If an initialization plasma state file is specified in the /static_state_data/ namelist,
!    that file is read as an aux state file and the EPA component data is copied to the current
!    plasma state.  
!
! Nota Bene:  The next item is under construction.  There is so much data that would have to
! be initialized and modeled that I will delay the implementation.  For example I don't
! want to try to generate psi(R,Z) or enter data for Rlim, Zlim by hand.
!
! 2) If no initialization plasma state file is specified, the EPA initialization
!    is done here from data in the /static_state_data/ namelist file.  All of the data in
!    the plasma state belongs to a specific component.  Most components only put data into
!    the plasma state the it owns.  However the EPA is distinguished in that it must set
!    some data owned by other components Ð for example the number, names, and power of ICRF 
!    sources; the number, types and even density profiles of rf minority species; the
!    number and names of neutral beam sources and their energies.
!
!    Note that there is a huge amount of EPA data which the model_EPA models do not use.  
!    Therefore only a subset of the EPA data is initialized and modeled. 
!    Typically this is only what is used by other components or that is watched by the 
!    monitor component.  More could be added at any time.

! Data owned by EPA that is initialized here:
!	
! 	MACHINE_DESCRIPTION data:  none
!	SHOT_CONFIGURATION DATA:
!		nspec_th		! number of thermal species
!       m_s, qatom_S, q_s, s_name ! specification of thermal species
!	SIMULATION_INIT: 
!		ALL and ALLA 	! composite species lists
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
! Data owned by EQ that is initialized here:
!	
! 	MACHINE_DESCRIPTION data:  
!		rlim, zlim		! R, Z points defining limiter
!
!	SHOT_CONFIGURATION DATA: none
!	SIMULATION_INIT: 
!		R_grid, Z_grid	! grid for psi(R,Z)
!		nrho_eq			! dimension of equilibrium rho grid
!		rho_eq			! equilibrium rho grid
!	STATE_DATA:
!		R_axis, Z_axis	! magnetic axis location
!	STATE_PROFILES:
!		curt			! integrated toroidal current inside rho
!		elong			! elongation (b/a)
!		triang			! triangularity (symmetrized)
!		grho2r2i		! metric coefficient needed to calculate li(3)
!		PsiRZ			! Poloidal flux versus R,Z
!		P_eq			! equilibrium scalar pressure
!		vol				! volume enclosed by rho
!
! Data owned by IC that is initialized here:
!	
! 	MACHINE_DESCRIPTION data: 
!		nicrf_src		! number of icrf sources
!		icrf_src_name	! names of icrf sources
!	SHOT_CONFIGURATION DATA:
!		nspec_rfmin		! number of icrf minority species
!       m_rfmin, qatom_rfmin, q_rfmin, rfmin_name ! specification of minority species
!	SIMULATION_INIT: 
!		fracmin			! fraction of electron densith in minorities (if kdens_rfmin = 'fraction')
!		kdens_rfmin		! specifies method to determine minority densities
!	STATE_DATA:
!		power_ic		! power on each icrf source
!	STATE_PROFILES:
!		nmini			! minority density profiles (maybe - note: IC sets nrho_icrf and rho_icrf)
!		
! Data owned by NBI that is initialized here:
!	
! 	MACHINE_DESCRIPTION data: 
!		nbeam		! number of beams
!		nbi_src_name	! names of beam sources
!	SHOT_CONFIGURATION DATA:
!		nbion		! names of ion species injected by each beam
!	SIMULATION_INIT: 
!		nspec_beam		! number of beam ion species
!       m_beam, qatom_beam, q_beam, rfmin_beam ! specification of beam species
!	STATE_DATA:
!		power_nbi		! power on each beam
!	STATE_PROFILES: none
!		
!
!
! One other point: EPA constructs all of the species lists and calls ps_merge_species(), which
! can only be done once.  So nspec_beam and the beam species must already be defined either
! by the input plasma state file as specified in the namelist or by the current state
! file previously initialized by EPA.  This affects this mdodel component through things
! like: nbeami(nrho_nbi-1, nspec_beam) and eperp_beami(nrho_nbi-1, nspec_beam) which are
! actually allocated here.
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
    
    INTEGER :: cclist(ps_ccount)   ! component activation list 

    CHARACTER (len=256) :: cur_state_file, cur_eqdsk_file, cur_jsdsk_file
    CHARACTER(len=32) :: mode
    CHARACTER(len=32) :: time_stamp
    
        !--------------------------------------------------------------------
        ! Some maximum array sizes to avoid having to do a lot of tedious allocates
        !--------------------------------------------------------------------
        
        integer, parameter :: nrho_max = 120
        integer, parameter :: nrho_nbi_max = 120
        integer, parameter :: n_spec_th_max = 5
        integer, parameter :: n_spec_max = 7
        integer, parameter :: n_spec_nonMax_max = 2
        integer, parameter :: n_nbi_src_max = 2
        integer, parameter :: n_spec_beam_max = 2


	!--------------------------------------------------------------------------
	!   Internal data
	!--------------------------------------------------------------------------
    CHARACTER (len=256) :: input_state_file
    INTEGER :: nzone
    REAL(KIND=rspec), ALLOCATABLE :: zone_center(:)

    INTEGER :: nrho_nbi

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

    namelist /static_state_data/ input_state_file, nrho_nbi
                      
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

	!------------------------------------------------------------------------------------
	!   Get command line arguments
	!------------------------------------------------------------------------------------

      call get_arg_count(iarg)
      if(iarg .ne. 5) then
        print*, 'model_EPA usage: '
        print*, ' command line args = cur_state_file cur_eqdsk_file  &
            &         cur_jsdsk_filemode mode time_stamp'
        stop 'incorrect command line arguments'
      end if
      
      call getarg(1,cur_state_file)
      call getarg(2,cur_eqdsk_file)
      call getarg(3, cur_jsdsk_file)
      call getarg(4,mode)
      call getarg(5,time_stamp)
      
     WRITE (*,*)
     WRITE (*,*) 'model_EPA'      
     print*, 'cur_state_file = ', trim(cur_state_file)
     print*, 'cur_eqdsk_file = ', trim(cur_eqdsk_file)
     print*, 'cur_jsdsk_file = ', trim(cur_jsdsk_file)
     print*, 'mode = ', trim(mode)
     print*, 'time_stamp = ', trim(time_stamp)
      
	!------------------------------------------------------------------------------------
	!  Get current plasma state 
	!------------------------------------------------------------------------------------
			
    call ps_get_plasma_state(ierr, trim(cur_state_file))
    if(ierr .ne. 0) then
       print*, 'model_EPA:failed to get_plasma_state'
       stop 1
    end if

	!--------------------------------------------------------------------------
	!    Open input namelist file
	!--------------------------------------------------------------------------

    OPEN (unit=21, file = 'model_EPA_input.nml', status = 'old',   &
         form = 'formatted', iostat = ierr)
    IF (ierr .ne. 0) STOP 'cannot open EPA_model.nml'
    
!------------------------------------------------------------------------------------
!     
!  Section 1: INIT function
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'INIT') THEN
              
    WRITE(*,*) 'model_EPA: INIT'

    READ (21, nml = static_state_data)
    CLOSE (21)
    WRITE (*, nml = static_state_data)

	!--------------------------------------------------------------------------
	!   1.1 Case: initialization data taken from input_state_file
	!--------------------------------------------------------------------------

    IF (TRIM(input_state_file) .ne. ' ') THEN
    
    	print*, 'Getting EPA source arrays from plasma state file ', TRIM(input_state_file)
          
		!---------------------------------------------------------------------------------
		!  Get aux plasma state to copy EPA data from.  It comes in from input_state_file
		!---------------------------------------------------------------------------------
			
		  CALL ps_get_plasma_state(ierr, TRIM(input_state_file), aux)
		  if(ierr .ne. 0) stop 'call failed to ps_get_plasma_state for aux state'
	
		!--------------------------------------------------------------------------
		!   Copy data in EPA sections of 'aux' state to current (i.e. 'ps') state. Note that
		!   this retains the rho_nbi grid and initial profile data from the input plasma 
		!   state.
		!--------------------------------------------------------------------------
		
                WRITE (*,*) 'aux%nrho_nbi = ', aux%nrho_nbi
                WRITE (*,*) 'aux%rho_nbi = ', aux%rho_nbi
 
		 CALL ps_cclist_remove("*", cclist, iwarn)
		 CALL ps_cclist_add("EPAI" ,cclist, iwarn)
		 CALL PS_COPY_DIMS(aux, ps, 0 ,ierr, cclist)
		 if(ierr .ne. 0) stop 'call failed to PS_COPY_DIMS for aux state to ps state'
		 CALL PS_COPY_PLASMA_STATE(aux, ps, ierr,  cclist)
		 if(ierr .ne. 0) stop 'call failed to PS_COPY_PLASMA_STATE for aux state to ps state'
			
		 print*,   '   number of nbi sources = ', ps%nbeam
		 print*,   '   names of sources = ', ps%nbi_src_name
		 print*,   '   nbi source powers = ', ps%power_nbi
		 print*,   '   number of beam species = ', ps%nspec_beam
		 print*,   '   names of beam species = ', ps%nbion
		 print*,   '   beam energy (keV) = ', ps%kvolt_nbi
		 print*,   '   eperp_beam allocated = ', ALLOCATED(ps%eperp_beami)
		 print*,   '   epll_beam allocated = ', ALLOCATED(ps%epll_beami)

                 WRITE (*,*) 'ps%nrho_nbi = ', ps%nrho_nbi
                 WRITE (*,*) 'ps%rho_nbi = ', ps%rho_nbi
                 WRITE (*,*) 'cclist = ', cclist

	ELSE
	
	!--------------------------------------------------------------------------
	!   1.2 Case: initialization data taken from namelist
	!--------------------------------------------------------------------------

		!---------------------------------------------------------------------------------
		!  Check if nbeam dimensioned arrays are allocated in current plasma state
		!---------------------------------------------------------------------------------
					 
		IF ( allocated(ps%nbi_src_name) ) THEN
			 print*, 'EPA source arrays allocated in initial Plasma State'             
			
			 print*,   '   number of nbi sources = ', ps%nbeam
			 print*,   '   names of sources = ', ps%nbi_src_name
			 print*,   '   nbi source powers = ', ps%power_nbi
			 print*,   '   number of beam species = ', ps%nspec_beam
			 print*,   '   names of beam species = ', ps%nbion
			 print*,   '   beam energy (keV) = ', ps%kvolt_nbi
			 print*,   '   eperp_beam allocated = ', ALLOCATED(ps%eperp_beami)
			 print*,   '   epll_beam allocated = ', ALLOCATED(ps%epll_beami)
	
		!-------------------------------------------------------------------------- 
		! Otherwise quit
		!--------------------------------------------------------------------------
			
		ELSE		  
			 print*, 'EPA sources NOT allocated in initial Plasma State'
			 print*, 'Can''t be done in this component without stepping on species lists'
			 print*, 'model_nb_2.f90: I give up'
			 STOP
			 			
		END IF

			
		!-------------------------------------------------------------------------- 
		! Allocate rho_nbi dimensioned profile arrays in plasma state
		!--------------------------------------------------------------------------
				 
		IF ( allocated(ps%rho_nbi) ) THEN
		
			 print*, "model_EPA: rho_nbi(:) already allocated in plasma state &
			 &        Error - I'm supposed to do that"
			 STOP
		ELSE
			ps%nrho_nbi = nrho_nbi  
			CALL ps_alloc_plasma_state(ierr)
			WRITE (*,*) 'model_EPA: Allocated nbi profiles in plasma state'      

			CALL rho_grid(ps%nrho_nbi,ps%rho_nbi)
			write (*,*) 'ps%rho_nbi = ', ps%rho_nbi

			!--------------------------------------------------------------------------
			! Initialize EPA output data to zero
			!--------------------------------------------------------------------------
			
			ps%nbeami = 0.
			ps%pbe = 0.        
			ps%pbi = 0.
			ps%pbth = 0.        
			ps%curbeam = 0.
			ps%epll_beami = 0.
			ps%eperp_beami = 0.
				  
		END IF

	END IF  ! End of initialization from namelist
            
    
    !-------------------------------------------------------------------------- 
    !1.4  Store initial plasma state for EPA
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)
        IF (ierr .ne. 0) THEN
		STOP 'model_EPA INIT: PS_STORE_PLASMA_STATE failed'
	ELSE
        	WRITE (*,*) "model_EPA: Stored initial Plasma State"    
        END IF
                    
END IF  ! End INIT function

!------------------------------------------------------------------------------------
!     
!  Section 2: STEP function - Change state data and store plasma state
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'STEP') THEN

    WRITE(*,*)
    WRITE(*,*) 'model_EPA: STEP'

    !---------------------------------------------------------------------------------
    !  2.1 Get data for profile and time evolution models from input_namelist_file.
    !      Note: At present there is no time evolution model for EPA data.  Everything 
    !      is constant, or zero, or is proportional to the beam power, i.e. responds 
    !      instantly.
    !---------------------------------------------------------------------------------
	
	READ (21, nml=evolving_model_data)
	CLOSE (21)

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
            CALL SWIM_error ('allocation', 'model_EPA' , 'zone_nbi_center')
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

	WRITE (*,*) "model_EPA: Stored Plasma State"    

	CALL sleep(3)   ! Wait for elvis to pick up
	
END IF ! End of STEP function

!--------------------------------------------------------------------------
!
!   Finalize function
!
!--------------------------------------------------------------------------

IF (TRIM(mode) .eq. 'FINALIZE') THEN
                                                                                                
     WRITE (*,*) "model_EPA: Finalize called"
     
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


END PROGRAM model_EPA
