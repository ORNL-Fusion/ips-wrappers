PROGRAM model_NB

! version 0.2 8/5/2009 (Batchelor)  Added beam species densities.
! version 0.0 3/13/2009 (Batchelor)

!--------------------------------------------------------------------------
!
!   Simple mock NB code for testing.
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
    
    CHARACTER (len=256) :: cur_state_file, cur_eqdsk_file, cur_cql_file, cur_jsdsk_file
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
!
!   Internal data
!
!--------------------------------------------------------------------------

    INTEGER :: nzone
    REAL(KIND=rspec), ALLOCATABLE :: zone_center(:)
    INTEGER :: i, j, index
    REAL(KIND=rspec) :: t
    REAL(KIND=rspec) :: n_beam_peak, alpha_n_beam
    REAL(KIND=rspec) :: Pth_e_beam_MW, P_th_i_beam_MW, alpha_P_th_e, alpha_P_th_i
    REAL(KIND=rspec) :: I_beam_MA, alpha_I_beam
    REAL(KIND=rspec) :: Epll_peak, Eperp_peak, alpha_Epll, alpha_Eperp

    INTEGER :: nbeam, nspec_beam, nrho_nbi
    
    CHARACTER*32 :: nbi_src_name(n_nbi_src_max)
    REAL(KIND=rspec) :: power_nbi(n_nbi_src_max)
    REAL(KIND=rspec) ::kvolt_nbi(n_nbi_src_max)
!    REAL(KIND=rspec) :: rho_nbi(nrho_nbi_max)

!------------------------------------------------------------------------------------
!     
!   Model input namelist
!
!------------------------------------------------------------------------------------

    namelist/NB_model/ nbeam, nbi_src_name, nrho_nbi, power_nbi, &
    		n_beam_peak, alpha_n_beam, &
            Pth_e_beam_MW, P_th_i_beam_MW, alpha_P_th_e, alpha_P_th_i, &
            I_beam_MA, alpha_I_beam, Epll_peak, Eperp_peak, alpha_Epll, alpha_Eperp
    

!------------------------------------------------------------------------------------
!     
!   Get command line arguments
!
!------------------------------------------------------------------------------------

      call get_arg_count(iarg)
    if(iarg .ne. 5) then
        print*, 'model_NB usage: '
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
     WRITE (*,*) 'model_NB'      
     print*, 'cur_state_file = ', trim(cur_state_file)
     print*, 'cur_eqdsk_file = ', trim(cur_eqdsk_file)
     print*, 'cur_jsdsk_file = ', trim(cur_jsdsk_file)
     print*, 'mode = ', trim(mode)
     print*, 'time_stamp = ', trim(time_stamp)
      
!------------------------------------------------------------------------------------
!     
!  Get current plasma state 
!
!------------------------------------------------------------------------------------
        
      call ps_get_plasma_state(ierr, trim(cur_state_file))
      if(ierr .ne. 0) then
        print*, 'model_NB:failed to get_plasma_state'
    stop 1
     end if


!--------------------------------------------------------------------------
!    Get NB model definition from namelist
!--------------------------------------------------------------------------

      OPEN (unit=21, file = 'model_NB_input.nml', status = 'old',   &
         form = 'formatted', iostat = ierr)
      IF (ierr .ne. 0) STOP 'cannot open NB_model.nml'

      READ (21, nml = NB_model)
      CLOSE (21)

    
!------------------------------------------------------------------------------------
!     
!  INIT function
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'INIT') THEN
              
    WRITE(*,*) 'model_NB: INIT'

    !--------------------------------------------------------------------------
    !
    !   Initialize plasma state data
    !
    !--------------------------------------------------------------------------
             
    IF ( allocated(ps%power_nbi) ) THEN
        print*, 'NB sources allocated in initial Plasma State'             
        
    ELSE
      
         print*, 'NB sources NOT allocated in initial Plasma State: doing it in model_NB'

    !-------------------------------------------------------------------------- 
    ! Allocate nbi source arrays in plasma state and load from namelist
    !--------------------------------------------------------------------------
         
        ps%nbeam = nbeam
   
        CALL ps_alloc_plasma_state(ierr)
        WRITE (*,*) 'model_NB: Allocated nbi source arrays'
        
        ps%nbi_src_name = nbi_src_name(1:nbeam)
        ps%power_nbi = power_nbi(1:nbeam)
        ps%kvolt_nbi = kvolt_nbi
        
    END IF
            
         print*,   '   number of nbi sources = ', ps%nbeam
         print*,   '   names of sources = ', ps%nbi_src_name
         print*,   '   nbi source powers = ', ps%power_nbi
         print*,   '   names of beam species = ', ps%nbion
         print*,   '   beam energy (keV) = ', ps%kvolt_nbi
         print*,   '   eperp_beam allocated = ', ALLOCATED(ps%eperp_beami)
         print*,   '   epll_beam allocated = ', ALLOCATED(ps%epll_beami)

         print*,   '   number of beam species = ', ps%nspec_beam
            
    !-------------------------------------------------------------------------- 
    ! Allocate nbi profile arrays in plasma state
    !--------------------------------------------------------------------------
             
    IF ( allocated(ps%rho_nbi) ) THEN
    
         print*, "model_NB: rho_nbi(:) already allocated in plasma state &
         &        Error - I'm supposed to do that"
         STOP
    ELSE
        ps%nrho_nbi = nrho_nbi  
        CALL ps_alloc_plasma_state(ierr)
        WRITE (*,*) 'model_NB: Allocated nbi profiles in plasma state'      
        print*,   'model_NB: nrho_nbi = ', ps%nrho_nbi
              
        END IF
    
    !--------------------------------------------------------------------------    !
    ! Initialize NB output data
    !--------------------------------------------------------------------------
        
        CALL rho_grid(ps%nrho_nbi,ps%rho_nbi)
       
        ps%nbeami = 0.
       
        ps%pbe = 0.
        
        ps%pbi = 0.
        
        ps%pbth = 0.
        
        ps%curbeam = 0.
        
        ps%epll_beami = 0.
        
        ps%eperp_beami = 0.
        
    
    !-------------------------------------------------------------------------- 
    ! Store initial plasma state for RF
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)

        WRITE (*,*) "model_NB: Stored initial Plasma State"    
               
                    
END IF  ! End INIT function

!------------------------------------------------------------------------------------
!     
!  STEP function - Change state data and store plasma state
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'STEP') THEN

    WRITE(*,*) 'model_NB: STEP'


    !--------------------------------------------------------------------------
    !
    !    nbi power deposition profiles 
    !
    !--------------------------------------------------------------------------
         
        nzone = nrho_nbi -1
        ALLOCATE( zone_center(nzone), stat=istat )
        IF (istat /= 0 ) THEN
            CALL SWIM_error ('allocation', 'model_NB' , 'zone_nbi_center')
            ierr = istat
            STOP
        END IF  

        zone_center = ( ps%rho_nbi(1:nrho_nbi-1) + ps%rho_nbi(2:nrho_nbi) )/2.
            
		! direct nbi power to thermal electrons
		CALL profgen_norm(alpha_P_th_e, 0.0_rspec, zone_center, ps%pbe)
		ps%pbe = 1.0e6_rspec * Pth_e_beam_MW * ps%pbe
		
		! direct nbi power to thermal ions
		CALL profgen_norm(alpha_P_th_i, 0.0_rspec, zone_center, ps%pbi)
		ps%pbi = 1.0e6_rspec * P_th_i_beam_MW * ps%pbi
		
		! beam current profile
		CALL profgen_norm(alpha_I_beam, 0.0_rspec, zone_center, ps%curbeam)
		ps%curbeam = 1.0e6_rspec * I_beam_MA * ps%pbi
		
		Do i = 1, ps%nspec_beam

	 		! beam species densities
			CALL profgen_norm(alpha_n_beam, 0.0_rspec, zone_center, ps%nbeami(:,i))
			ps%nbeami(:,i) =  n_beam_peak * ps%nbeami(:,i)
		
                !  changed multipler to make NB come out in keV	
	 		! beam Epll profile
			CALL profgen_norm(alpha_Epll, 0.0_rspec, zone_center, ps%epll_beami(:,i))
			ps%epll_beami(:,i) =  Epll_peak * ps%epll_beami(:,i)

			! beam Eperp profile
			CALL profgen_norm(alpha_Eperp, 0.0_rspec, zone_center, ps%eperp_beami(:,i))
			ps%eperp_beami(:,i) =  Eperp_peak * ps%eperp_beami(:,i)
				
		END DO

    !--------------------------------------------------------------------------    !
    ! Store the data in plasma_state file
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)

        WRITE (*,*) "model_NB: Stored Plasma State"    

        CALL sleep(3)   ! Wait for elvis to pick up
END IF

!--------------------------------------------------------------------------
!
!   Finalize function
!
!--------------------------------------------------------------------------

IF (TRIM(mode) .eq. 'FINALIZE') THEN
                                                                                                
     WRITE (*,*) "model_NB: Finalize called"
     
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



END PROGRAM model_NB
