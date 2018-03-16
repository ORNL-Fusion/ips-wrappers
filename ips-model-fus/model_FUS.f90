PROGRAM model_FUS

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
  

    USE plasma_state_mod
    

    use swim_global_data_mod, only : &
            & rspec, ispec, &               ! int: kind specification for real and integer
            & SWIM_name, SWIM_filename, &   ! derived data types: containing one character
                                                                        ! string
            & SWIM_error                    ! subroutine: a simple error handling routine
    
    IMPLICIT none
    
    INTEGER :: ierr = 0
    INTEGER :: istat, iarg
    INTEGER :: ii, iout
    
    INTEGER :: system
    
    CHARACTER (len=256) :: cur_state_file, cur_eqdsk_file, cur_cql_file
    CHARACTER(len=32) :: mode
    CHARACTER(len=32) :: time_stamp
    
        !--------------------------------------------------------------------
        ! Some maximum array sizes to avoid having to do a lot of tedious allocates
        !--------------------------------------------------------------------
        
        integer, parameter :: n_spec_fus_max = 2


!--------------------------------------------------------------------------
!
!   Internal data
!
!--------------------------------------------------------------------------

    INTEGER :: nzone
    REAL(KIND=rspec), ALLOCATABLE :: zone_center(:)
    INTEGER :: i, j, index
    REAL(KIND=rspec) :: t
    REAL(KIND=rspec) :: n_fus_peak, alpha_n_fus
    REAL(KIND=rspec) :: Pth_e_fus_MW, P_th_i_fus_MW, alpha_P_th_e, alpha_P_th_i
    REAL(KIND=rspec) :: I_fus_MA, alpha_I_fus
    REAL(KIND=rspec) :: Epll_peak, Eperp_peak, alpha_Epll, alpha_Eperp

    INTEGER :: nspec_fusion, nrho_fus

    CHARACTER*32 :: SFUS_name(n_spec_fus_max)
    
!------------------------------------------------------------------------------------
!     
!   Model input namelist
!
!------------------------------------------------------------------------------------

    namelist/FUS_model/ nspec_fusion,SFUS_name, nrho_fus, &
    		n_fus_peak, alpha_n_fus, &
            Pth_e_fus_MW, P_th_i_fus_MW, alpha_P_th_e, alpha_P_th_i, &
            I_fus_MA, alpha_I_fus, Epll_peak, Eperp_peak, alpha_Epll, alpha_Eperp
    

!------------------------------------------------------------------------------------
!     
!   Get command line arguments
!
!------------------------------------------------------------------------------------

      call get_arg_count(iarg)
    if(iarg .ne. 4) then
        print*, 'model_FUS usage: '
        print*, ' command line args = cur_state_file cur_eqdsk_file  &
            &         cur_jsdsk_filemode mode time_stamp'
        stop 'incorrect command line arguments'
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
    stop 1
     end if


!--------------------------------------------------------------------------
!    Get FUS model definition from namelist
!--------------------------------------------------------------------------

      OPEN (unit=21, file = 'model_FUS_input.nml', status = 'old',   &
         form = 'formatted', iostat = ierr)
      IF (ierr .ne. 0) STOP 'cannot open FUS_model.nml'

      READ (21, nml = FUS_model)
      CLOSE (21)

    
!------------------------------------------------------------------------------------
!     
!  INIT function
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'INIT') THEN
              
    WRITE(*,*) 'model_FUS: INIT'

    !--------------------------------------------------------------------------
    !
    !   Initialize plasma state data
    !
    !--------------------------------------------------------------------------
             
    IF ( allocated(ps%SFUS_name) ) THEN
        print*, 'FUS species allocated in initial Plasma State'             
        
    ELSE
      
         print*, 'FUS species NOT allocated in initial Plasma State: doing it in model_FUS'

    !-------------------------------------------------------------------------- 
    ! Allocate fus source arrays in plasma state and load from namelist
    !--------------------------------------------------------------------------
         
        ps%nspec_fusion = nspec_fusion
   
        CALL ps_alloc_plasma_state(ierr)
        WRITE (*,*) 'model_FUS: Allocated fus source arrays'
        
        ps%SFUS_name = SFUS_name(1:nspec_fusion)
        
    END IF
            
         print*,   '   number of fusion species = ', ps%nspec_fusion
         print*,   '   names of fusion species = ', ps%SFUS_name
         print*,   '   eperp_fus allocated = ', ALLOCATED(ps%eperp_fusi)
         print*,   '   epll_fus allocated = ', ALLOCATED(ps%epll_fusi)
            
    !-------------------------------------------------------------------------- 
    ! Allocate fus profile arrays in plasma state
    !--------------------------------------------------------------------------
             
    IF ( allocated(ps%rho_fus) ) THEN
    
         print*, "model_FUS: rho_fus(:) already allocated in plasma state &
         &        Error - I'm supposed to do that"
         STOP
    ELSE
        ps%nrho_fus = nrho_fus  
        CALL ps_alloc_plasma_state(ierr)
        WRITE (*,*) 'model_FUS: Allocated fus profiles in plasma state'      
        print*,   'model_FUS: nrho_fus = ', ps%nrho_fus
              
        END IF
    
    !--------------------------------------------------------------------------    !
    ! Initialize FUS output data
    !--------------------------------------------------------------------------
        
        CALL rho_grid(ps%nrho_fus,ps%rho_fus)
       
        ps%nfusi = 0.
       
        ps%pfuse = 0.
        
        ps%pfusi = 0.
        
        ps%pbth = 0.
        
        ps%curfusn = 0.
        
        ps%epll_fusi = 0.
        
        ps%eperp_fusi = 0.
        
    
    !-------------------------------------------------------------------------- 
    ! Store initial plasma state for RF
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)

        WRITE (*,*) "model_FUS: Stored initial Plasma State"    
               
                    
END IF  ! End INIT function

!------------------------------------------------------------------------------------
!     
!  STEP function - Change state data and store plasma state
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'STEP') THEN

    WRITE(*,*) 'model_FUS: STEP'


    !--------------------------------------------------------------------------
    !
    !    fus power deposition profiles 
    !
    !--------------------------------------------------------------------------
         
        nzone = nrho_fus -1
        ALLOCATE( zone_center(nzone), stat=istat )
        IF (istat /= 0 ) THEN
            CALL SWIM_error ('allocation', 'model_FUS' , 'zone_fus_center')
            ierr = istat
            STOP
        END IF  

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
		
                !  changed multipler to make FUS come out in keV	
	 		! fus Epll profile
			CALL profgen_norm(alpha_Epll, 0.0_rspec, zone_center, ps%epll_fusi(:,i))
			ps%epll_fusi(:,i) =  Epll_peak * ps%epll_fusi(:,i)/SUM(ps%epll_fusi(:,i))

			! fus Eperp profile
			CALL profgen_norm(alpha_Eperp, 0.0_rspec, zone_center, ps%eperp_fusi(:,i))
			ps%eperp_fusi(:,i) =  Eperp_peak * ps%eperp_fusi(:,i)/SUM(ps%eperp_fusi(:,i))
				
		END DO

    !--------------------------------------------------------------------------    !
    ! Store the data in plasma_state file
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)

        WRITE (*,*) "model_FUS: Stored Plasma State"    

        CALL sleep(3)   ! Wait for elvis to pick up
END IF

!--------------------------------------------------------------------------
!
!   Finalize function
!
!--------------------------------------------------------------------------

IF (TRIM(mode) .eq. 'FINALIZE') THEN
                                                                                                
     WRITE (*,*) "model_FUS: Finalize called"
     
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
