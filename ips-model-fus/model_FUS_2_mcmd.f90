PROGRAM model_FUS

! version 2.0 12/7/2009 (Batchelor)  Added capability to initialize from an input plasma state
!             file. Changed to write partial plasma state update file for compatibility with MCMD.
!             Added internal state data file so component can have memory from one call to the
!             next.
! version 0.0 8/5/2009 (Batchelor)

!--------------------------------------------------------------------------
!
!   Simple mock FUS code for testing.
!
!       The code requires 4 command-line arguments
!       1) path to the current plasma state file
!       2) path to the current plasma eqdsk file
!       3) action mode, i.e. one of: "INIT", "STEP", or "FINALIZE"
!       4) time stamp the time set by the driver component to which the simulation is 
!          supposed to advance.
!   
!		Note: NUBEAM calculates both NBI and fusion fast ion physics.  However the
!		Plasma State lists a fusion fast ion component (FUS).  For this model it's 
!		simpler just to write a new FUS component than to add the fusion stuff to the
!		model neutral beam component.  This code is almost identical to model_NB.f90
!		with the substitution nb -> fus.
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
!   The data to be loaded is derived from input namelist file model_NB_input.nml.
!   This file must contain two namelists:
!
!   /static_state_data/ -> Contains plasma state data that goes directly into the state, or
!						   data used to initialize the FUS part of the plasma state (e.g. 
!						   an input plasma state file to copy)
!
!   /evolving_model_data/ -> Contains data that defines the ad hoc models used to give time 
!                            varying numbers for the FUS data (e.g names of different models
!							 and parameters to be used in the models, like profiles shapes or
!							 parameters to define time variation
!
! So that the models for time dependence can have some memory a second namelist: 
! /internal_state_data/ is written to file 'internal_state_data.nml'
!
! There are two ways to allocate the plasma state arrays for the FUS component:
!
! 1) If an initialization plasma state file is specified in the /static_state_data/ namelist,
!    that file is read as an aux state file and the FUS component data is merged to the current
!    plasma state.  If that input file is properly initialized then the current state will 
!    inherit a complete set of FUS data and allocated arrays.  The models must then work with
!    the allocated array sizes - including everthing dimensioned by nrho_fus
!
! 2) If no initialization plasma state file is specified, it is assumed that this is part of
!    a normal initialization sequence in which the EPA has set the arrays assiciated with
!    the NB species.  This includes nspec_fusion and the FUS species arrays dimensioned with it.
!
!    This model component then specifies nrho_fus and proceeds to allocate and initialize
!    all the variables dimensioned by nrho_fus:
!
!      nrho_fus, rho_fus(nrho_fus), nfusi(nrho_fus-1,nspec_fus), pfuse(nrho_fus-1),
!      pfusi(nrho_fus-1), epll_fusi(nrho_fus-1,nspec_fus), eperp_fusi(nrho_fus-1, nspec_fus),
!      curfusn(nrho_fus-1)
!
!
! One other point: EPA constructs all of the species lists and calls ps_merge_species(), which
! can only be done once.  So nspec_fus and the fusion species must already be defined either
! by the input plasma state file as specified in the namelist or by the current state
! file previously initialized by EPA.  These can't be specified here.
!
! This affects this model component through things like: nfusi(nrho_nbi-1, nspec_fus) and 
! eperp_fusi(nrho_nbi-1, nspec_fus) which are actually allocated here.
!
!--------------------------------------------------------------------------
 

    USE plasma_state_mod
    

    use swim_global_data_mod, only : &
            & rspec, ispec, &               ! int: kind specification for real and integer
            & SWIM_name, SWIM_filename, &   ! derived data types: containing one character
                                                                        ! string
            & SWIM_error, &                    ! subroutine: a simple error handling routine
            & SWIM_filename_length, &       ! standard length for long file names = 256 now
            & SWIM_name_length              ! standard length for item names = 31 now
    
    IMPLICIT none
    
    INTEGER :: ierr = 0, iwarn = 0
    INTEGER :: istat, iarg
    INTEGER :: ii, iout
    
    INTEGER :: cclist(ps_ccount)   ! component activation list   
 
    CHARACTER (len=256) :: cur_state_file, cur_eqdsk_file
    CHARACTER(len=32) :: mode
    CHARACTER(len=32) :: time_stamp
    
	!--------------------------------------------------------------------
	! Some maximum array sizes to avoid having to do a lot of tedious allocates
	!--------------------------------------------------------------------
	
	integer, parameter :: n_spec_fus_max = 2


	!--------------------------------------------------------------------------
	!   Internal data
	!--------------------------------------------------------------------------
    CHARACTER (len = SWIM_filename_length) :: input_state_file ! to initialize FUS arrays
    INTEGER :: nzone
    REAL(KIND=rspec), ALLOCATABLE :: zone_center(:)
    INTEGER :: i, j, index
    REAL(KIND=rspec) :: t
    INTEGER :: nspec_fusion, nrho_fus

	!--------------------------------------------------------------------------
	!
	!   Data for evolving profile models
	!
	!--------------------------------------------------------------------------

    CHARACTER(len = SWIM_filename_length) :: input_namelist_file = 'model_FUS_input.nml'
    CHARACTER(len = SWIM_filename_length) :: internal_state_file = 'internal_state_data.nml'
    CHARACTER (len = SWIM_name_length) :: FUS_profile_model_name

    INTEGER :: n_step_calls = 0
   
    REAL(KIND=rspec) :: n_fus_peak, alpha_n_fus
    REAL(KIND=rspec) :: Pth_e_fus_MW, P_th_i_fus_MW, alpha_P_th_e, alpha_P_th_i
    REAL(KIND=rspec) :: I_fus_MA, alpha_I_fus
    REAL(KIND=rspec) :: Epll_peak, Eperp_peak, alpha_Epll, alpha_Eperp

    CHARACTER*32 :: SFUS_name(n_spec_fus_max)
    
	!------------------------------------------------------------------------------------     
	!   Model input namelist
	!------------------------------------------------------------------------------------

    namelist /static_state_data/ input_state_file, nrho_fus
    
    namelist/evolving_model_data/ &
            FUS_profile_model_name, &
    		n_fus_peak, alpha_n_fus, &
            Pth_e_fus_MW, P_th_i_fus_MW, alpha_P_th_e, alpha_P_th_i, &
            I_fus_MA, alpha_I_fus, Epll_peak, Eperp_peak, alpha_Epll, alpha_Eperp

    namelist /internal_state_data/ &
          n_step_calls

!------------------------------------------------------------------------------------
!     
!   Get command line arguments
!
!------------------------------------------------------------------------------------

      call get_arg_count(iarg)
    if(iarg .ne. 4) then
        print*, 'model_FUS usage: '
        print*, ' command line args = cur_state_file cur_eqdsk_file  &
            &  mode time_stamp'
        call exit(1)
    end if
      
      call getarg(1,cur_state_file)
      call getarg(2,cur_eqdsk_file)
      call getarg(3,mode)
      call getarg(4,time_stamp)
      
     WRITE (*,*)
     WRITE (*,*) 'model_FUS'      
     print*, 'cur_state_file = ', trim(cur_state_file)
     print*, 'cur_eqdsk_file = ', trim(cur_eqdsk_file)
     print*, 'mode = ', trim(mode)
     print*, 'time_stamp = ', trim(time_stamp)
      
!------------------------------------------------------------------------------------
!     
!  Get current plasma state 
!
!------------------------------------------------------------------------------------
        
      call ps_get_plasma_state(ierr, trim(cur_state_file))
      if(ierr .ne. 0) then
        print*, 'model_FUS:failed to get_plasma_state'
        call exit(1)
      end if

	!--------------------------------------------------------------------------
	!    Open input namelist file
	!--------------------------------------------------------------------------

    OPEN (unit=21, file = input_namelist_file, status = 'old',   &
         action='read', form = 'formatted', iostat = ierr)
            IF (ierr /= 0 ) THEN
                CALL SWIM_error ('open', 'model_FUS_2_mcmd.f90',TRIM(input_namelist_file))
                ierr = istat
                print*, 'cannot open input_namelist_file'
                call exit(1)
            END IF
    
!------------------------------------------------------------------------------------
!     
!  Section 1: INIT function
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'INIT') THEN
              
    WRITE(*,*) 'model_FUS: INIT'

    READ (21, nml = static_state_data)
    CLOSE (21)
    WRITE (*, nml = static_state_data)

	!--------------------------------------------------------------------------
	!   1.1 Case: initialization data taken from input_state_file
	!--------------------------------------------------------------------------

    IF (TRIM(input_state_file) .ne. ' ') THEN
    
    	print*, 'Getting NB source arrays from plasma state file ', TRIM(input_state_file)
          
		!---------------------------------------------------------------------------------
		!  Get aux plasma state to copy FUS data from.  It comes in from input_state_file
		!---------------------------------------------------------------------------------
			
		  CALL ps_get_plasma_state(ierr, TRIM(input_state_file), aux)
		  if(ierr .ne. 0) then
		      print*, 'call failed to ps_get_plasma_state for aux state'
		      call exit(1)
		  end if
	
		!--------------------------------------------------------------------------
		!   Copy data in FUS sections of 'aux' state to current (i.e. 'ps') state. Note that
		!   this retains the rho_fus grid and initial profile data from the input plasma 
		!   state.
		!--------------------------------------------------------------------------
		 
			CALL ps_cclist_remove("*", cclist, iwarn)
			CALL ps_cclist_add("FUS" ,cclist, iwarn)
			CALL PS_COPY_DIMS(aux, ps, 0 ,ierr, cclist)
			if(ierr .ne. 0) then
			    print*, 'call failed to PS_COPY_DIMS for aux state to ps state'
			    call exit(1)
			end if
			CALL PS_COPY_PLASMA_STATE(aux, ps, ierr,  cclist)
			if(ierr .ne. 0) then
			    print*, 'call failed to PS_COPY_PLASMA_STATE for aux state to ps state'
			    call exit(1)
			end if
			
			 WRITE (*,*) 'ps%nrho_fus = ', ps%nrho_fus
			 WRITE (*,*) 'ps%rho_fus = ', ps%rho_fus
			 WRITE (*,*) 'ps%SFUS_name =', ps%SFUS_name
			 WRITE (*,*) 'eperp_fus allocated = ', ALLOCATED(ps%eperp_fusi)
			 WRITE (*,*) 'epll_fus allocated = ', ALLOCATED(ps%epll_fusi)
			 WRITE (*,*) 'cclist = ', cclist
	ELSE
	
	!--------------------------------------------------------------------------
	!   1.2 Case: initialization data taken from namelist
	!--------------------------------------------------------------------------

		!---------------------------------------------------------------------------------
		!  Check if nspce_fusions dimensioned arrays are allocated in current plasma state
		!---------------------------------------------------------------------------------
					 
		IF ( allocated(ps%SFUS_name) ) THEN
			 print*, 'FUS species arrays allocated in initial Plasma State'             
			
			 WRITE (*,*) 'ps%nrho_fus = ', ps%nrho_fus
			 WRITE (*,*) 'ps%rho_fus = ', ps%rho_fus
			 WRITE (*,*) 'ps%SFUS_name =', ps%SFUS_name
			 WRITE (*,*) 'eperp_fus allocated = ', ALLOCATED(ps%eperp_fusi)
			 WRITE (*,*) 'epll_fus allocated = ', ALLOCATED(ps%epll_fusi)
	
		!-------------------------------------------------------------------------- 
		! Otherwise quit
		!--------------------------------------------------------------------------
			
		ELSE		  
			 print*, 'FUS species NOT allocated in initial Plasma State'
			 print*, 'Can''t be done in this component without stepping on species lists'
			 print*, 'model_FUS_2.f90: I give up'
			 call exit(1)
			 			
		END IF

			
		!-------------------------------------------------------------------------- 
		! Allocate rho_fus dimensioned profile arrays in plasma state
		!--------------------------------------------------------------------------
				 
		IF ( allocated(ps%rho_fus) ) THEN
		
			 print*, "model_FUS: rho_fus(:) already allocated in plasma state &
			 &        Error - I'm supposed to do that"
			 call exit(1)
		ELSE
			ps%nrho_fus = nrho_fus  
			CALL ps_alloc_plasma_state(ierr)
			WRITE (*,*) 'model_FUS: Allocated fus profiles in plasma state'      

			CALL rho_grid(ps%nrho_fus,ps%rho_fus)
			write (*,*) 'ps%rho_fus = ', ps%rho_fus
		END IF

   	END IF ! End of selection of source for initialization data 
    !--------------------------------------------------------------------------    !
    ! 1.3 Initialize FUS output data
    !--------------------------------------------------------------------------
       
        ps%nfusi = 0.
        ps%pfuse = 0.
        ps%pfusi = 0.
        ps%curfusn = 0.
        ps%epll_fusi = 0.
        ps%eperp_fusi = 0.
        
    
    !-------------------------------------------------------------------------- 
    !1.4  Store initial plasma state for FUS
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)
			IF (ierr .ne. 0) THEN
			    print*, 'model_FUS INIT: PS_STORE_PLASMA_STATE failed'
			    call exit(1)
		ELSE
        	WRITE (*,*) "model_FUS: Stored initial Plasma State"    
        END IF

    !---------------------------------------------------------------------------------
    !  Create internal state file for time dependent model data
    !---------------------------------------------------------------------------------

       OPEN (unit=21, file=TRIM(internal_state_file), status='unknown', &
            action='write', iostat=istat, form='formatted')
            IF (istat /= 0 ) THEN
                CALL SWIM_error ('open', 'model_FUS_2_mcmd.f90',TRIM(internal_state_file))
                ierr = istat
                print*, 'cannot open internal state namelist file for output'
                call exit(1)
            END IF

        write(21, nml=internal_state_data)
        CLOSE (21)
               
                    
END IF  ! End INIT function

!------------------------------------------------------------------------------------
!     
!  Section 2: STEP function - Change state data and store plasma state
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'STEP') THEN

    WRITE(*,*)
    WRITE(*,*) 'model_FUS: STEP'

    !---------------------------------------------------------------------------------
    !  2.1 Get data for profile and time evolution models from input_namelist_file.
    !      Note: At present there is no time evolution model for FUS data.  Everything 
    !      is constant, or zero, or is proportional to the beam power, i.e. responds 
    !      instantly.
    !---------------------------------------------------------------------------------
	
	READ (21, nml=evolving_model_data)
	CLOSE (21)
    
    !---------------------------------------------------------------------------------
    !  2.2 Get internal state data for time dependent model
    !---------------------------------------------------------------------------------

       OPEN (unit=21, file=TRIM(internal_state_file), status='old', &
            action='read', iostat=istat, form='formatted')
            IF (istat /= 0 ) THEN
                CALL SWIM_error ('open', 'model_FUS_2_mcmd.f90',TRIM(internal_state_file))
                ierr = istat
                print*, 'cannot open input_namelist_file'
                call exit(1)
            END IF

        read(21, nml=internal_state_data)
        CLOSE (21)
        
        n_step_calls = n_step_calls + 1 ! increment counter for number of calls
        WRITE (*,*) 'Call ', n_step_calls, ' to model_FUS_2_mcmd step'

    !-------------------------------------------------------------------------- 
    !  2.3  model = const.  Don't touch the FUS data.  Pass plasma state through.
    !--------------------------------------------------------------------------

	IF (TRIM(FUS_profile_model_name) == 'const') THEN
	
		CONTINUE
		
    !-------------------------------------------------------------------------- 
    !  2.4  model = zeros.  Set all outputs to zero
    ! --------------------------------------------------------------------------

	ELSE IF (TRIM(FUS_profile_model_name) == 'zeros') THEN

        ps%nfusi = 0.
        ps%pfuse = 0.
        ps%pfusi = 0.
        ps%curfusn = 0.
        ps%epll_fusi = 0.
        ps%eperp_fusi = 0.	
		
    !-------------------------------------------------------------------------- 
    !  2.5  model = profiles.  Actually generate the profiles from the models
    !--------------------------------------------------------------------------

	ELSE IF (TRIM(FUS_profile_model_name) == 'profiles') THEN
         
        nzone = ps%nrho_fus - 1
        ALLOCATE( zone_center(nzone), stat=istat )
        IF (istat /= 0 ) THEN
            CALL SWIM_error ('allocation', 'model_FUS' , 'zone_fus_center')
            ierr = istat
            call exit(1)
        END IF  
        nrho_fus = ps%nrho_fus
        zone_center = ( ps%rho_fus(1:nrho_fus-1) + ps%rho_fus(2:nrho_fus) )/2.
            
		! direct fus power to thermal electrons
		CALL profgen_norm(alpha_P_th_e, 0.0_rspec, zone_center, ps%pfuse)
		ps%pfuse = 1.0e6_rspec * Pth_e_fus_MW * ps%pfuse/SUM(ps%pfuse)
	
		! direct fus power to thermal ions
		CALL profgen_norm(alpha_P_th_i, 0.0_rspec, zone_center, ps%pfusi)
		ps%pfusi = 1.0e6_rspec * P_th_i_fus_MW * ps%pfusi/SUM(ps%pfusi)
		
		! fus current profile
		CALL profgen_norm(alpha_I_fus, 0.0_rspec, zone_center, ps%curfusn)
		ps%curfusn = 1.0e6_rspec * I_fus_MA * ps%pfusi/SUM(ps%pfusi)
		
		Do i = 1, ps%nspec_fusion

	 		! fus species densities
			CALL profgen_norm(alpha_n_fus, 0.0_rspec, zone_center, ps%nfusi(:,i))
			ps%nfusi(:,i) =  n_fus_peak * ps%nfusi(:,i)/SUM(ps%nfusi(:,i))
	                
	 		! fus Epll profile
			CALL profgen_norm(alpha_Epll, 0.0_rspec, zone_center, ps%epll_fusi(:,i))
			ps%epll_fusi(:,i) =  Epll_peak * ps%epll_fusi(:,i)/SUM(ps%epll_fusi(:,i))

			! fus Eperp profile
			CALL profgen_norm(alpha_Eperp, 0.0_rspec, zone_center, ps%eperp_fusi(:,i))
			ps%eperp_fusi(:,i) =  Eperp_peak * ps%eperp_fusi(:,i)/SUM(ps%eperp_fusi(:,i))
				
		END DO
		
    END IF ! End of cases of different profile models
    
    !--------------------------------------------------------------------------    !
    ! 2.6 Store the data in plasma_state file
    !--------------------------------------------------------------------------

	CALL PS_WRITE_UPDATE_FILE('FUS_'// cur_state_file, ierr)

	WRITE (*,*) "model_FUS_2_mcmd: Stored Partial Plasma State"    
        
    !---------------------------------------------------------------------------------
    !  Write internal state data for time dependent model
    !---------------------------------------------------------------------------------

       OPEN (unit=21, file=TRIM(internal_state_file), status='old', &
            action='write', iostat=istat, form='formatted')
            IF (istat /= 0 ) THEN
                CALL SWIM_error ('open', 'model_FUS_2_mcmd.f90',TRIM(internal_state_file))
                ierr = istat
                print*, 'cannot open internal_state_file for output'
                call exit(1)
            END IF

        write(21, nml=internal_state_data)
        CLOSE (21)

        CALL sleep(3)   ! Wait for elvis to pick up
	
END IF ! End of STEP function

!--------------------------------------------------------------------------
!
!  Section 3: Finalize function
!
!--------------------------------------------------------------------------

IF (TRIM(mode) .eq. 'FINALIZE') THEN
                                                                                                
     WRITE (*,*) "model_FUS: Finalize called"
     
ENDIF

!--------------------------------------------------------------------------
!
!  Section 4: Internal routines
!
!--------------------------------------------------------------------------

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


    !------------------------------------
    SUBROUTINE profgen_norm(f0,fa,rho,fans)

      !  quick & dirty profile generator

      REAL(KIND=rspec), intent(in) :: f0,fa  ! core and edge values
      REAL(KIND=rspec), dimension(:), intent(in) :: rho  ! normalized sqrt tor. flux
      REAL(KIND=rspec), dimension(:), intent(out) :: fans  ! formula output
      REAL(KIND=rspec) :: norm  ! normalization to make integral rho*d(rho) = one

      norm = (2*fa + f0)/6
      fans = (f0-fa)*(1.0_rspec-rho**2)**2 + fa
      fans = fans/norm

    END SUBROUTINE profgen_norm



END PROGRAM model_FUS
