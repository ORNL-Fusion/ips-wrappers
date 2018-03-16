PROGRAM model_RF_IC

! version 2.0 8/31/2009 (Batchelor)  Added capability to initialize from an input plasma state
!             file. Added a model routine that gives off-axis peaked profiles: Lorentz_Linear.
!             Generally cleaned up.
! version 0.0 9/11/2008 (Batchelor)

!--------------------------------------------------------------------------
!
!   Simple mock RF_IC code for testing.
!
!       The code requires 6 command-line arguments
!       1) path to the current plasma state file
!       2) path to the current plasma eqdsk file
!       3) path to current cql distribution function file - cur_cql_file
!       4) path to current quasilinear operator file - cur_dql_file
!       5) action mode, i.e. one of: "INIT", "STEP", or "FINALIZE"
!       6) time stamp the time set by the driver component to which the simulation is 
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
!   The data to be loaded is derived from input namelist file model_RF_IC_input.nml.
!   This file must contain two namelists:
!
!   /static_state_data/ -> Contains plasma state data that goes directly into the state, or
!						   data used to initialize the RF_IC part of the plasma state (e.g. a 
!						   plasma state file to copy)
!
!   /evolving_model_data/ -> Contains data that defines the ad hoc models used to give time 
!                            varying numbers for the RF_IC data (e.g names of different models
!							 and parameters to be used in the models, like profiles shapes or
!							 parameters to define time variation
!
! There are two ways to allocate the plasma state arrays for the RF_IC component:
!
! 1) If an initialization plasma state file is specified in the /static_state_data/ namelist,
!    that file is read as an aux state file and the RF_IC component data is merged to the 
!    current plasma state.  If that input file is properly initialized then the current state 
!    will inherit a complete set of RF_IC data and allocated arrays.  The models must then
!    work with the allocated array sizes.
!
! 2) If no initialization plasma state file is specified, it is assumed that this is part of
!    a normal initialization sequence in which the EPA has set the arrays assiciated with
!    the icrf sources: 

!       nicrf_src, icrf_src_name(nicrf_src), power_ic(nicrf_src), nspec_rfmin, 
!       and the RFMIN species arrays. 

!    This model component then specifies nrho_icrf and proceeds to allocate and initialize
!    all the variables dimensioned by nrho_icrf:
!
!      nrho_icrf, rho_icrf(nrho_icrf), picrf_srcs(nrho_icrf-1, nicrf_src,0:nspc_alla)
!      picrt_totals(nrho_icrf-1, 0:nspc_alla), picth(nrho_icrf-1), curich(nrho_icrf-1)
!      eperp_mini(nrho_icrf-1, nspec_rfmin)
!
!    Note that there is a huge amount of RF_IC data which the model_RF_IC models do not use.  
!    Therefore only the subset of the RF_IC data listed above is initialized and modeled. 
!    Typically this is only what is used by the EPA component or that is watched by the 
!    monitor component.
!
! One other point: EPA constructs all of the species lists and calls ps_merge_species(), which
! can only be done once.  So nspec_rfmin and the minority species must already be defined, 
! either by the input plasma state file as specified in the namelist or by the current state
! file previously initialized by EPA.  This affects this model component through things
! like: nmini(nrho_icrf-1, nspec_rfmin) and eperp_mini(nrho_icrf-1, nspec_rfmin) which are
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
    INTEGER :: i, j
    
    INTEGER :: cclist(ps_ccount)   ! component activation list 

    CHARACTER (len=256) :: cur_state_file, cur_eqdsk_file, cur_cql_file, cur_dql_file
    CHARACTER(len=32) :: mode
    CHARACTER(len=32) :: time_stamp
    
        !--------------------------------------------------------------------
        ! Some maximum array sizes to avoid having to do a lot of tedious allocates
        !--------------------------------------------------------------------
        
        integer, parameter :: nrho_max = 120
        integer, parameter :: nrho_icrf_max = 120
        integer, parameter :: n_spec_th_max = 5
        integer, parameter :: n_spec_max = 7
        integer, parameter :: n_spec_nonMax_max = 2
        integer, parameter :: n_icrf_src_max = 2
        integer, parameter :: n_spec_rfmin_max = 2


	!--------------------------------------------------------------------------
	!   Internal data
	!--------------------------------------------------------------------------
    CHARACTER (len=256) :: input_state_file
    INTEGER :: nzone
    REAL(KIND=rspec), ALLOCATABLE :: zone_center(:)

    INTEGER :: nrho_icrf
     REAL(KIND=rspec) :: freq_ic(n_icrf_src_max), power_ic(n_icrf_src_max)

    CHARACTER(len=32) :: RF_IC_profile_model_name
    
    ! namelist parameters for Lorentz_Linear model:
    !   rho_max = rho of peak of the Lorentzian (not exactly the peak of the profile)
    !   w = width of Lorentzian
    !   f0 = value of normalized profile on axis, rho = 0
    !   f1 = value of normalized profile at rho = 1
    REAL(KIND=rspec) :: nmini_peak, rho_max_nmini, w_nmini, f0_nmini, f1_nmini
    REAL(KIND=rspec) :: FP_th_e_icrf,rho_max_P_th_e, w_P_th_e, f0_P_th_e, f1_P_th_e
    REAL(KIND=rspec) :: FP_th_i_icrf, rho_max_P_th_i, w_P_th_i, f0_P_th_i, f1_P_th_i
    REAL(KIND=rspec) :: I_icrf_MA, rho_max_I_icrf, w_I_icrf, f0_I_icrf, f1_I_icrf

    REAL(KIND=rspec) :: FP_min
	
	!------------------------------------------------------------------------------------
	!   Model input namelists
	!------------------------------------------------------------------------------------

    namelist /static_state_data/ input_state_file, freq_ic, nrho_icrf
                      
    namelist /evolving_model_data/ &
          RF_IC_profile_model_name, &
          nmini_peak, rho_max_nmini, w_nmini, f0_nmini, f1_nmini, &
          FP_th_e_icrf, rho_max_P_th_e, w_P_th_e, f0_P_th_e, f1_P_th_e, &
          FP_th_i_icrf, rho_max_P_th_i, w_P_th_i, f0_P_th_i, f1_P_th_i, &
          I_icrf_MA, rho_max_I_icrf, w_I_icrf, f0_I_icrf, f1_I_icrf   
    
!------------------------------------------------------------------------------------
!     
!  Section 0: Setup
!
!------------------------------------------------------------------------------------

	!------------------------------------------------------------------------------------
	!   Get command line arguments
	!------------------------------------------------------------------------------------

      call get_arg_count(iarg)
      if(iarg .ne. 6) then
         print*, 'model_RF_IC usage: '
         print*, ' command line args = cur_state_file cur_eqdsk_file cur_cql_file, &
             &         cur_dql_filemode mode time_stamp'
         stop 'incorrect command line arguments'
      end if
      
      call getarg(1,cur_state_file)
      call getarg(2,cur_eqdsk_file)
      call getarg(3, cur_cql_file)
      call getarg(4, cur_dql_file)
      call getarg(5,mode)
      call getarg(6,time_stamp)
      
     WRITE (*,*)
     WRITE (*,*) 'model_RF_IC'      
     print*, 'cur_state_file = ', trim(cur_state_file)
     print*, 'cur_eqdsk_file = ', trim(cur_eqdsk_file)
     print*, 'cur_cql_file = ', trim(cur_cql_file)
     print*, 'cur_dql_file = ', trim(cur_dql_file)
     print*, 'mode = ', trim(mode)
     print*, 'time_stamp = ', trim(time_stamp)
      
	!------------------------------------------------------------------------------------
	!  Get current plasma state 
	!------------------------------------------------------------------------------------
			
    call ps_get_plasma_state(ierr, trim(cur_state_file))
    if(ierr .ne. 0) then
       print*, 'model_RF_IC:failed to get_plasma_state'
       stop 1
    end if

	!--------------------------------------------------------------------------
	!    Open input namelist file
	!--------------------------------------------------------------------------

    OPEN (unit=21, file = 'model_RF_IC_input.nml', status = 'old',   &
         form = 'formatted', iostat = ierr)
    IF (ierr .ne. 0) STOP 'cannot open RF_IC_model.nml'
    
!------------------------------------------------------------------------------------
!     
!  Section 1: INIT function
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'INIT') THEN
              
    WRITE(*,*) 'model_RF_IC: INIT'

    READ (21, nml = static_state_data)
    CLOSE (21)
    WRITE (*, nml = static_state_data)

	!--------------------------------------------------------------------------
	!   1.1 Case: initialization data taken from input_state_file
	!--------------------------------------------------------------------------

    IF (TRIM(input_state_file) .ne. ' ') THEN
    
    	print*, 'Getting RF_IC source arrays from plasma state file ', TRIM(input_state_file)
          
		!---------------------------------------------------------------------------------
		!  Get aux plasma state to copy RF_IC data from.  It comes in from input_state_file
		!---------------------------------------------------------------------------------
			
		  CALL ps_get_plasma_state(ierr, TRIM(input_state_file), aux)
		  if(ierr .ne. 0) stop 'call failed to ps_get_plasma_state for aux state'
	
		!--------------------------------------------------------------------------
		!   Copy data in RF_IC sections of 'aux' state to current (i.e. 'ps') state. Note that
		!   this retains the rho_icrf grid and initial profile data from the input plasma 
		!   state.
		!--------------------------------------------------------------------------
		
                WRITE (*,*) 'aux%nrho_icrf = ', aux%nrho_icrf
                WRITE (*,*) 'aux%rho_icrf = ', aux%rho_icrf
 
		 CALL ps_cclist_remove("*", cclist, iwarn)
		 CALL ps_cclist_add("IC" ,cclist, iwarn)
		 CALL PS_COPY_DIMS(aux, ps, 0 ,ierr, cclist)
		 if(ierr .ne. 0) stop 'call failed to PS_COPY_DIMS for aux state to ps state'
		 CALL PS_COPY_PLASMA_STATE(aux, ps, ierr,  cclist)
		 if(ierr .ne. 0) stop 'call failed to PS_COPY_PLASMA_STATE for aux state to ps state'
			
         print*,   '   number of icrf sources = ', ps%nicrf_src
         print*,   '   names of sources = ', ps%icrf_src_name
         print*,   '   icrf frequencies = ', ps%freq_ic
         print*,   '   icrf source powers = ', ps%power_ic
		 print*,   '   number of minority species = ', ps%nspec_rfmin
		 print*,   '   names of minority species = ', ps%nbion
		 print*,   '   eperp_icrf allocated = ', ALLOCATED(ps%eperp_mini)
		 print*,   '   epll_icrf allocated = ', ALLOCATED(ps%epll_mini)

                 WRITE (*,*) 'ps%nrho_icrf = ', ps%nrho_icrf
                 WRITE (*,*) 'ps%rho_icrf = ', ps%rho_icrf
                 WRITE (*,*) 'cclist = ', cclist

	ELSE
	
	!--------------------------------------------------------------------------
	!   1.2 Case: initialization data taken from namelist
	!--------------------------------------------------------------------------

		!---------------------------------------------------------------------------------
		!  Check if nicrf dimensioned arrays are allocated in current plasma state
		!---------------------------------------------------------------------------------
					 
		IF ( allocated(ps%icrf_src_name) ) THEN
			 print*, 'RF_IC source arrays allocated in initial Plasma State'  
			 
			 ps%freq_ic = freq_ic(1:ps%nicrf_src) ! EPA doesn't set source frequency
				
			 print*,   '   number of icrf sources = ', ps%nicrf_src
			 print*,   '   names of sources = ', ps%icrf_src_name
			 print*,   '   icrf frequencies = ', ps%freq_ic
			 print*,   '   icrf source powers = ', ps%power_ic
			 print*,   '   number of minority species = ', ps%nspec_rfmin
			 print*,   '   names of minority species = ', ps%nbion
			 print*,   '   eperp_icrf allocated = ', ALLOCATED(ps%eperp_mini)
			 print*,   '   epll_icrf allocated = ', ALLOCATED(ps%epll_mini)
	
		!-------------------------------------------------------------------------- 
		! Otherwise quit
		!--------------------------------------------------------------------------
			
		ELSE		  
			 print*, 'RF_IC sources NOT allocated in initial Plasma State'
			 print*, 'Can''t be done in this component without stepping on species lists'
			 print*, 'model_RF_IC_2.f90: I give up'
			 STOP
			 			
		END IF

			
		!-------------------------------------------------------------------------- 
		! Allocate rho_icrf dimensioned profile arrays in plasma state
		!--------------------------------------------------------------------------
				 
		IF ( allocated(ps%rho_icrf) ) THEN
		
			 print*, "model_RF_IC: rho_icrf(:) already allocated in plasma state &
			 &        Error - I'm supposed to do that"
			 STOP
		ELSE
			ps%nrho_icrf = nrho_icrf  
			CALL ps_alloc_plasma_state(ierr)
			WRITE (*,*) 'model_RF_IC: Allocated nbi profiles in plasma state'      

			CALL rho_grid(ps%nrho_icrf,ps%rho_icrf)
			write (*,*) 'ps%rho_icrf = ', ps%rho_icrf

			!--------------------------------------------------------------------------
			! Initialize RF_IC output data to zero
			!--------------------------------------------------------------------------
			
			ps%picrf_srcs = 0.        
			ps%picrf_totals = 0.        
			ps%picth = 0.
			ps%curich = 0.
				  
		END IF

	END IF  ! End of initialization from namelist
            
    
    !-------------------------------------------------------------------------- 
    !1.4  Store initial plasma state for RF_IC
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)
        IF (ierr .ne. 0) THEN
		STOP 'model_RF_IC INIT: PS_STORE_PLASMA_STATE failed'
	ELSE
        	WRITE (*,*) "model_RF_IC: Stored initial Plasma State"    
        END IF
                    
END IF  ! End INIT function

!------------------------------------------------------------------------------------
!     
!  Section 2: STEP function - Change state data and store plasma state
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'STEP') THEN

    WRITE(*,*)
    WRITE(*,*) 'model_RF_IC: STEP'

    !---------------------------------------------------------------------------------
    !  2.1 Get data for profile and time evolution models from input_namelist_file.
    !      Note: At present there is no time evolution model for RF_IC data.  Everything 
    !      is constant, or zero, or is proportional to the icrf power, i.e. responds 
    !      instantly.
    !---------------------------------------------------------------------------------
	
	READ (21, nml=evolving_model_data)
	CLOSE (21)

    !-------------------------------------------------------------------------- 
    !  2.2  model = const.  Don't touch the RF_IC data.  Pass plasma state through.
    !--------------------------------------------------------------------------

	IF (TRIM(RF_IC_profile_model_name) == 'const') THEN
	
		CONTINUE
		
    !-------------------------------------------------------------------------- 
    !  2.3  model = zeros.  Set all outputs to zero
    ! --------------------------------------------------------------------------

	ELSE IF (TRIM(RF_IC_profile_model_name) == 'zeros') THEN

		ps%picrf_srcs = 0.        
		ps%picrf_totals = 0.        
		ps%picth = 0.
		ps%curich = 0.
	
		
    !-------------------------------------------------------------------------- 
    !  2.4  model = profiles.  Actually generate the profiles from the models
    !--------------------------------------------------------------------------

	ELSE IF (TRIM(RF_IC_profile_model_name) == 'profiles') THEN
         
        nzone = ps%nrho_icrf -1
        ALLOCATE( zone_center(nzone), stat=istat )
        IF (istat /= 0 ) THEN
            CALL SWIM_error ('allocation', 'model_RF_IC' , 'zone_nbi_center')
            ierr = istat
            STOP
        END IF  

        zone_center = ( ps%rho_icrf(1:nrho_icrf-1) + ps%rho_icrf(2:nrho_icrf) )/2.
           write(*,*) 'zone_center = ', zone_center

        ps%picrf_totals = 0.

        FP_min = 1.0 - FP_th_e_icrf - FP_th_i_icrf   ! minority gets what's left
        
        DO j = 1, ps%nicrf_src  ! Loop over ICRF sources
            
            WRITE (*,*) 'power_ic(', j, ') = ',  ps%nicrf_src
            
            ! direct ICRF power to thermal electrons
			CALL Lorentz_Linear(rho_max_P_th_e, w_P_th_e, f0_P_th_e, f1_P_th_e, &
			zone_center, ps%picrf_srcs(:, j, 0))
			ps%picrf_srcs(:, j, 0) = FP_th_e_icrf*ps%power_ic(j)*ps%picrf_srcs(:, j, 0)/ &
									 SUM(ps%picrf_srcs(:, j, 0) )
            ps%picrf_totals(:, 0) = ps%picrf_totals(:,0) + ps%picrf_srcs(:, j, 0)

            WRITE (*,*) 'ps%picrf_srcs(:,', j, '0) = ', ps%picrf_srcs(:, j, 0)
            WRITE (*,*) 'ps%picrf_totals(:, 0) = ', ps%picrf_srcs(:, j, 0)
            WRITE (*,*) 'SUM(ps%picrf_totals(:, 0)) = ', SUM(ps%picrf_totals(:, 0) )

            ! direct ICRF power to thermal ions (take at most lowest two)
            IF (ps%nspec_th >= 1) THEN
                DO i=1, MIN(ps%nspec_th, 2)
                    ! generate profiles
                    CALL Lorentz_Linear(rho_max_P_th_i, w_P_th_i, f0_P_th_i, f1_P_th_i, &
                    	zone_center, ps%picrf_srcs(:, j, i) )
                    ! take power evenly distributed among thermal ion species
                    ps%picrf_srcs(:, j, i) =FP_th_i_icrf *ps%power_ic(j)/MIN(ps%nspec_th, 2) &
                                    *ps%picrf_srcs(:, j, i)/SUM(ps%picrf_srcs(:, j, i) )

                    !sum over sources
                    ps%picrf_totals(:,i) = ps%picrf_totals(:,i) + ps%picrf_srcs(:,j,i)
                 END DO
	
				WRITE (*,*) 'ps%picrf_srcs(:,', j, ',', i,') = ', ps%picrf_srcs(:, j, i)
				WRITE (*,*) 'ps%picrf_totals(:,' , i,') = ', ps%picrf_srcs(:, j, i)
				WRITE (*,*) 'SUM(ps%picrf_totals(:,' , i,')) = ', SUM(ps%picrf_totals(:, i) )

            END IF
              
        END DO

        ps%picth(:) = 0. 
        DO i=1, MIN(ps%nspec_th, 2)
        	ps%picth(:) = ps%picth(:) + ps%picrf_totals(:,i)

            WRITE (*,*) 'SUM(ps%picth(:)) = ', SUM(ps%picth(:) )

        END DO
 
		! ICRF driven current
		CALL Lorentz_Linear(rho_max_I_icrf, w_I_icrf, f0_I_icrf, f1_I_icrf, &
							zone_center, ps%curich)
		ps%curich = 1.0e6_rspec * I_icrf_MA * ps%curich/SUM(ps%curich)

		
		DO i = 1, ps%nspec_rfmin

	 		! icrf species densities
			CALL Lorentz_Linear(rho_max_nmini, w_nmini, f0_nmini, f1_nmini, &
							    zone_center, ps%nmini(:,i))
			ps%nmini(:,i) =  nmini_peak * ps%nmini(:,i)/MAXVAL(ps%nmini(:,i))
						
		END DO

    END IF ! End of cases of different profile models
    
    !--------------------------------------------------------------------------    !
    ! Store the data in plasma_state file
    !--------------------------------------------------------------------------

	CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)

	WRITE (*,*) "model_RF_IC_2: Stored Plasma State"    

	CALL sleep(1)   ! Wait for elvis to pick up
	
END IF ! End of STEP function

!--------------------------------------------------------------------------
!
!   Finalize function
!
!--------------------------------------------------------------------------

IF (TRIM(mode) .eq. 'FINALIZE') THEN
                                                                                                
     WRITE (*,*) "model_RF_IC: Finalize called"
     
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


END PROGRAM model_RF_IC
