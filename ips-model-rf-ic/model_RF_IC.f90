PROGRAM model_RF_IC

! version 0.0 9/11/2008 (Batchelor)

!--------------------------------------------------------------------------
!
!   Simple mock ICRF wave code code for testing. This is based on an older, simpler
!   model RF component called change_profiles_RF_IC.
!
!       At present (3/7/08) it only loads data into:
!       1) ps%picrf_srsc(:,1,0) = power from ICRF source #1 into thermal electrons
!       2) ps%picrf_srsc(:,1,1) = power from ICRF source #1 into thermal ion species #1 
!          (usually D)
!       3) ps%picrf_totals(:, 0) = power from all ICRF sources into thermal electrons
!       4) ps%picrf_totals(:, 1) = power from all ICRF sources into thermal ion sepcies #1
!       5) ps%icth(:) = power from all ICRF sources into all thermal ion species
!
!       This code is lauched by the model_RF_IC.py component script.  Which is 
!       of class RF_IC
!
!
!       The code requires 6 command-line arguments
!       1) path to the current plasma state file
!       2) path to the current plasma eqdsk file
!       3) path to current cql distribution function file - cur_cql_file
!       4) path to current quasilinear operator file - cur_dql_file
!       5) action mode, i.e. one of: "INIT", "STEP", or "FINALIZE"
!       6) time stamp the time set by the driver component to which the simulation is 
!          supposed to advance.
!      supposed to advance.
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
!
!   Internal data
!
!--------------------------------------------------------------------------

    INTEGER :: nzone
    REAL(KIND=rspec), ALLOCATABLE :: zone_center(:)
    INTEGER :: i, j, index
    REAL(KIND=rspec) :: t
    REAL(KIND=rspec) :: Pfrac_th_e, Pfrac_th_i, Pfrac_min
        REAL(KIND=rspec) :: alpha_P_th_e, alpha_P_th_i, alpha_P_min

    INTEGER :: nicrf_src, n_spec_rfmin, nrho_icrf
    
    CHARACTER*32 :: icrf_src_name(n_icrf_src_max)
    REAL(KIND=rspec) :: freq_ic(n_icrf_src_max), power_ic(n_icrf_src_max)
    REAL(KIND=rspec) :: rho_icrf(nrho_icrf_max)

!------------------------------------------------------------------------------------
!     
!   Model input namelist
!
!------------------------------------------------------------------------------------

    namelist/rf_ic_model/ nicrf_src,icrf_src_name, freq_ic, nrho_icrf, &
                         Pfrac_th_i, Pfrac_th_e, alpha_P_th_e,alpha_P_th_i,alpha_P_min
    

!------------------------------------------------------------------------------------
!     
!   Get command line arguments
!
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
!     
!  Get current plasma state 
!
!------------------------------------------------------------------------------------
        
      call ps_get_plasma_state(ierr, trim(cur_state_file))
      if(ierr .ne. 0) then
        print*, 'model_RF_IC:failed to get_plasma_state'
    stop 1
     end if


!--------------------------------------------------------------------------
!    Get RF model definition from namelist
!--------------------------------------------------------------------------

      OPEN (unit=21, file = 'model_RF_IC_input.nml', status = 'old',   &
         form = 'formatted', iostat = ierr)
      IF (ierr .ne. 0) STOP 'cannot open rf_ic_model.nml'

      READ (21, nml = rf_ic_model)
      CLOSE (21)

    
!------------------------------------------------------------------------------------
!     
!  INIT function
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'INIT') THEN

    !--------------------------------------------------------------------------
    !
    !   Initialize plasma state data
    !
    !--------------------------------------------------------------------------
             
    IF ( allocated(ps%power_ic) ) THEN
        print*, 'RF sources allocated in initial Plasma State'
        ps%freq_ic = freq_ic(1:ps%nicrf_src) ! EPA doesn't set source frequency
        
    ELSE
      
         print*, 'RF sources NOT located in initial Plasma State: doing it in model_RF_IC'

    !-------------------------------------------------------------------------- 
    ! Allocate ICRF source arrays in plasma state and load from namelist
    !--------------------------------------------------------------------------
         
        ps%nicrf_src = nicrf_src
   
        CALL ps_alloc_plasma_state(ierr)
        WRITE (*,*) 'model_RF_IC: Allocated icrf source arrays'
        
        ps%icrf_src_name = icrf_src_name(1:nicrf_src)
        ps%freq_ic = freq_ic(1:nicrf_src)
        
    END IF
            
         print*,   '   number of icrf sources = ', ps%nicrf_src
         print*,   '   names of sources = ', ps%icrf_src_name
         print*,   '   icrf frequencies = ', ps%freq_ic
         print*,   '   icrf source powers = ', ps%power_ic
            
    !-------------------------------------------------------------------------- 
    ! Allocate ICRF profile arrays in plasma state
    !--------------------------------------------------------------------------
             
    IF ( allocated(ps%rho_icrf) ) THEN
    
         print*, "model_RF_IC: rho_icrf(:) already allocated in plasma state &
         &        Error - I'm supposed to do that"
         STOP
    ELSE
        ps%nrho_icrf = nrho_icrf  
        CALL ps_alloc_plasma_state(ierr)
        WRITE (*,*) 'model_RF_IC: Allocated icrf profiles in plasma state'      
        print*,   'model_RF_IC: nrho_icrf = ', ps%nrho_icrf
              
        END IF
    
    !--------------------------------------------------------------------------    !
    ! Initialize RF output data
    !--------------------------------------------------------------------------
        
        CALL rho_grid(ps%nrho_icrf,ps%rho_icrf)
       
        ps%picrf_srcs = 0.
        
        ps%picrf_totals = 0.
        
        ps%picth = 0.
        
        
    
    !-------------------------------------------------------------------------- 
    ! Store initial plasma state for RF
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)

        WRITE (*,*) "model_RF_IC: Stored initial Plasma State"    
               
                    
END IF  ! End INIT function

!------------------------------------------------------------------------------------
!     
!  STEP function - Change state data and store plasma state
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'STEP') THEN

    WRITE(*,*) 'change_profiles_RF_IC: STEP'
    write(*,*) 'step  ps%power_ic = ', ps%power_ic, ' ps%power_nbi = ', ps%power_nbi


    !--------------------------------------------------------------------------
    !
    !    ICRF power deposition profiles 
    !
    !--------------------------------------------------------------------------
        ps%picrf_totals(:, 0) = 0.
        ps%picth(:) = 0.            

        Pfrac_min = 1.0 - Pfrac_th_e - Pfrac_th_i   ! minority gets what's left
         
        nzone = nrho_icrf -1
        ALLOCATE( zone_center(nzone), stat=istat )
        IF (istat /= 0 ) THEN
            CALL SWIM_error ('allocation', 'model_RF_IC' , 'zone_icrf_center')
            ierr = istat
            STOP
        END IF  
        zone_center = ( ps%rho_icrf(1:nrho_icrf-1) + ps%rho_icrf(2:nrho_icrf) )/2.

!    write(*,*) 'top of loop   ps%power_ic = ', ps%power_ic, ' ps%power_nbi = ', ps%power_nbi

        
        DO j = 1, ps%nicrf_src  ! Loop over ICRF sources
            
            ! direct ICRF power to thermal electrons
            CALL profgen_norm(alpha_P_th_e, 0.0_rspec, zone_center, ps%picrf_srcs(:, j, 0))
            ps%picrf_srcs(:, j, 0) = Pfrac_th_e*ps%power_ic(j)*ps%picrf_srcs(:, j, 0)
            ps%picrf_totals(:, 0) = ps%picrf_totals(:,0) + ps%picrf_srcs(:, j, 0)
            
!    write(*,*) 'middle of loop   ps%power_ic = ', ps%power_ic, ' ps%power_nbi = ', ps%power_nbi

            ! direct ICRF power to thermal ions (take at most lowest two)
            IF (ps%nspec_th >= 1) THEN
                DO i=1, MIN(ps%nspec_th, 2)
                    ! generate profiles
                    CALL profgen_norm(alpha_P_th_i, 0.0_rspec, zone_center, ps%picrf_srcs(:, j, i) )
                    ! take power evenly distributed among thermal ion species
                    ps%picrf_srcs(:, j, i) = Pfrac_th_i*ps%power_ic(j)/MIN(ps%nspec_th, 2) &
                                             *ps%picrf_srcs(:, j, i)

                    !sum over sources
                    ps%picrf_totals(:,index) = ps%picrf_totals(:,index) + ps%picrf_srcs(:,j,index)
 !   write(*,*) 'species name = ', ps%S_name(i), ' ps%power_ic = ', ps%power_ic, ' ps%power_nbi = ', ps%power_nbi
                END DO
              END IF
              
          end do

    !--------------------------------------------------------------------------    !
    ! Store the data in plasma_state file
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)

        WRITE (*,*) "model_RF_IC: Stored Plasma State"    

        CALL sleep(3)   ! Wait for elvis to pick up
END IF

!--------------------------------------------------------------------------
!
!   Finalize function
!
!--------------------------------------------------------------------------

IF (TRIM(mode) .eq. 'FINALIZE') THEN
                                                                                                
     WRITE (*,*) "change_n_T_EPA: Finalize called"
     
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



END PROGRAM model_RF_IC
