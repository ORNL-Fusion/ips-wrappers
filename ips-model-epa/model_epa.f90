PROGRAM model_epa

! version 0.0 7/25/2008 (Batchelor)

!   A EPA component that puts specified data into the plasma state for testing purposes and to
!   provide a simple driver for other components.  It is intended to be incremental in that
!   initially it only loads data needed for the RF_IC and Fokker Planck components.  Other 
!   data can be added later.
!
!   This code is driven by the model_epa.py component script
!
!   The data to be loaded is derived and input namelist file model_epa_input.nml.  This file
!   must contain two namelists:
!   /static_state_data/ -> contains data that goes directly into the state
!   /evolving_model_data/ -> contains data that tells this code how to generate profiles and
!      how they change in time
!
!    This version requires 5 commandline arguments: 
!    1) current state file
!    2) Next state file
!    3) current eqdsk file
!    4) mode = one of "INIT", "STEP", "FINALIZE"
!    5) timeStamp = initial time for "INIT", or = time at end of time stamp for "STEP"
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
            & SWIM_error, &                 ! subroutine: a simple error handling routine
            & SWIM_filename_length, &       ! standard length for long file names = 256 now
            & SWIM_name_length              ! standard length for item names = 31 now

    
    IMPLICIT none

    INTEGER :: ierr = 0
    INTEGER :: istat, iarg
    INTEGER :: i, ii, iout
    
    INTEGER :: system
    
    CHARACTER (len=256) :: cur_state_file, next_state_file, cur_eqdsk_file
    CHARACTER (len=32) :: mode
    CHARACTER (len=32) :: time_stamp
    
        !--------------------------------------------------------------------
        ! Some maximum array sizes to avoid having to do a lot of tedious allocates
        !--------------------------------------------------------------------
        
        integer, parameter :: nrho_max = 120
        integer, parameter :: n_spec_th_max = 5
        integer, parameter :: n_spec_max = 7
        integer, parameter :: n_spec_nonMax_max = 2
        integer, parameter :: n_icrf_src_max = 2
        integer, parameter :: n_spec_RFmin_max = 2
    

!--------------------------------------------------------------------------
!
!   Time invariant PLASMA STATE DATA (Set at INIT only)
!
!--------------------------------------------------------------------------

   
    !-----------------------------------
    ! Machine description
    !-----------------------------------
    
        INTEGER :: nicrf_src                    ! number of RF sources
        
        CHARACTER (len = SWIM_name_length) :: &
                icrf_src_name(n_icrf_src_max)   ! names of ICRF sources
       
        CHARACTER (len = SWIM_filename_length) :: &
                ant_model_src(n_icrf_src_max)   !file names for antenna models

        REAL (KIND = rspec) :: r_min_box, r_max_box, z_min_box, z_max_box  ! location 
                                                                          ! of bounding box
   
    !-----------------------------------
    ! Shot configuration data (see relevant sections of Plasma State spec.dat)
    !-----------------------------------

        !-----------------------------------
        ! Particle Species
        !-----------------------------------
        
        INTEGER :: nspec       ! number of ion species = nspec_th + nspec_beam, rfmin, fus
    
        !-----------------------------------
        ! Main (thermal) Plasma Species
        !-----------------------------------
        
        INTEGER :: nspec_th                   ! number of thermal ion species
        
        CHARACTER (len = SWIM_name_length) ::  &
            s_name(n_spec_max)              ! names of thermal ion species, (1:nspec_th)
        
        INTEGER :: &
               q_s(n_spec_th_max),        & ! Integer charge of ion species s
               m_s(n_spec_th_max)           ! Integer mass of ion species s  
        
        INTEGER ::   nrho     ! number of rho values in thermal species density grid
        
        REAL (KIND = rspec) ::  rho(nrho_max) ! rho values in density grid, (1:nrho_n)
    
        !-----------------------------------
        ! ICRF minority ions
        !-----------------------------------
        INTEGER ::    nspec_rfmin         ! number of rf minority species (on rf grid)
    
       
        CHARACTER (len = SWIM_name_length) ::  &
            rfmin_name(n_spec_RFmin_max)           ! names of rf min species
    
        INTEGER ::           &
            q_RFMIN(n_spec_RFmin_max),           & ! Integer charge of species s, (0:nspec_th)
            m_RFMIN(n_spec_RFmin_max)              ! Integer mass of species s, (0:nspec_th)

    
        REAL (KIND = rspec) ::  rho_icrf(nrho_max)  ! icrf rho grid
   
    !-----------------------------------
    ! Simulation init data (see relevant sections of Plasma State spec.dat)
    !-----------------------------------
        
        CHARACTER (len = 32) :: kdens_rfmin = 'fraction'
        REAL (KIND = rspec) ::  fracmin(n_spec_RFmin_max) ! fraction of electron density
                                                        ! for eacn minority species
        

!--------------------------------------------------------------------------
!
!   STATE DATA (Time varying PLASMA data)
!
!--------------------------------------------------------------------------

       
    !-----------------------------------
    ! Time at beginning and end of time step
    !-----------------------------------
    REAL (KIND = rspec) ::  &
           t0,                  &   ! time at beginning of step [msec]
           t1                       ! time at end of step [msec]

    !-----------------------------------
    ! Main (thermal) Plasma Species
    !-----------------------------------
 
        REAL (KIND = rspec) :: &
            n_s(nrho_max, 0:n_spec_th_max), & ! density profile of species s
            t_s(nrho_max, 0:n_spec_th_max), & ! temperature profiles
            zeff(nrho_max),           & ! effective charge state profile
            m_impurity(nrho_max)        ! effective impurity mass profile

    
    !-----------------------------------
    ! ICRF minority ions
    !-----------------------------------
 
        REAL (KIND = rspec) ::      &
            nmini(nrho_max, n_spec_RFmin_max),        & ! minority ion densities
            eperp_mini(nrho_max, n_spec_RFmin_max),   & ! minority ion perp energy moments
            epll_mini(nrho_max, n_spec_RFmin_max)       ! minority ion parallel energy moment
    
        CHARACTER (len = SWIM_filename_length) :: &
            dist_fun_RFmin_file, & ! distribution function of non-Maxwellian  
                                   ! species s, (1:nspec_nonMax) N.B. For now a distribution
                                   ! function type is a file name
            ql_operator_file       ! Quasilinear operator file  

    !--------------------------------------------------------------------------
    ! Equilibrium
    !--------------------------------------------------------------------------

        CHARACTER (len = SWIM_filename_length) :: model_eqdsk_file ! to be added to state

        REAL (kind = rspec) :: r_axis, z_axis
        
    !--------------------------------------------------------------------------
    ! RF Control Data
    !--------------------------------------------------------------------------
   
        REAL (kind = rspec) ::    power_ic(n_icrf_src_max)
    


!--------------------------------------------------------------------------
!
!   Data for evolving profile models
!
!--------------------------------------------------------------------------

    CHARACTER(len = SWIM_filename_length) :: input_namelist_file = 'model_epa_input.nml'

    REAL(KIND=rspec) :: ne0, ne_edge, alpha_ne, delta_ne 
    REAL(KIND=rspec) :: Te0, Te_edge, alpha_Te, delta_Te
    REAL(KIND=rspec) :: Ti0, Ti_edge, alpha_Ti, delta_Ti
    
    INTEGER :: nzone
    REAL(KIND=rspec), ALLOCATABLE :: zone_center(:)
    
    REAL(KIND=rspec) :: t
    INTEGER :: i_t
    
    REAL(KIND=rspec) ::  P_icrf_0(n_icrf_src_max), P_icrf_1(n_icrf_src_max), t_icrf_1
       
!------------------------------------------------------------------------------------
!     
!   Namelists
!
!------------------------------------------------------------------------------------

    namelist /static_state_data/ &
          r_min_box, r_max_box, z_min_box, z_max_box, &
          nspec_th, s_name, q_s, m_s, nrho, &
          zeff, &
          nspec_rfmin, q_RFMIN, m_RFMIN, &
          nicrf_src, icrf_src_name, &
          kdens_rfmin, fracmin, &
          model_eqdsk_file
    
    namelist /evolving_model_data/ &
          r_axis, z_axis,  &
          ne0, ne_edge, alpha_ne, delta_ne, &
          Te0, Te_edge, alpha_Te, delta_Te, &
          Ti0, Ti_edge, alpha_Ti, delta_Ti, &
          P_icrf_0, P_icrf_1, t_icrf_1
       
!------------------------------------------------------------------------------------
!     
!   Get command line arguments
!
!------------------------------------------------------------------------------------

    iarg = command_argument_count()
    if(iarg .ne. 5) then
        print*, 'model_epa: '
         print*, ' command line args = cur_state_file next_state_file cur_eqdsk_file  mode time_stamp'
         stop 'incorrect command line arguments'
    end if
      
      call get_command_argument(1,cur_state_file)
      call get_command_argument(2,next_state_file)
      call get_command_argument(3,cur_eqdsk_file)
      call get_command_argument(4,mode)
      call get_command_argument(5,time_stamp)
      
    WRITE (*,*)
    WRITE (*,*) 'model_epa'      
      print*, 'cur_state_file = ', trim(cur_state_file)
      print*, 'next_state_file = ', trim(next_state_file)
      print*, 'cur_eqdsk_file = ', trim(cur_eqdsk_file)
      print*, 'mode = ', trim(mode)
      print*, 'time_stamp = ', trim(time_stamp)
      


    
!------------------------------------------------------------------------------------
!     
!  INIT function
!
!------------------------------------------------------------------------------------

IF (TRIM(mode) == 'INIT') THEN
          
    WRITE (*,*)
    WRITE (*,*) 'model_epa INIT'
          
    !---------------------------------------------------------------------------------
    !     
    !  Get current plasma state 
    !
    !---------------------------------------------------------------------------------
        
      call   ps_get_plasma_state(ierr, trim(cur_state_file))
      if(ierr .ne. 0) stop 'call failed to ps_get_plasma_state to load profiles '
          
    !---------------------------------------------------------------------------------
    !     
    !  Get static EPA plasma state data from input_namelist_file
    !
    !---------------------------------------------------------------------------------

       OPEN (unit=21, file=TRIM(input_namelist_file), status='old', &
            action='read', iostat=istat, form='formatted')
            IF (istat /= 0 ) THEN
                CALL SWIM_error ('open', 'model_epa.f90',TRIM(input_namelist_file))
                ierr = istat
                stop 'cannot open input_namelist_file'

            END IF
        ierr = 0

        read(21, nml=static_state_data)
        CLOSE (21)
      
        
        WRITE (*, nml = static_state_data)

    !--------------------------------------------------------------------------
    !
    !   Initialize and allocate species arrays and thermal profile grids
    !
    !--------------------------------------------------------------------------
    
        ps%nspec_th = nspec_th
        ps%nrho = nrho
        ps%nspec_rfmin = nspec_rfmin
           
        WRITE (*,*) 'model_epa: About to allocate species and thermal profile arrays'
        CALL    ps_alloc_plasma_state(ierr)
        WRITE (*,*) 'model_epa: Allocated species and thermal profile arrays'
            
        !-----------------------------------
        ! Particle Species
        !-----------------------------------
        
            !-----------------------------------
            ! Main (thermal) Plasma Species
            !-----------------------------------

         WRITE (*,*) "Converting thermal electrons"
 
                CALL    ps_species_convert(-1, -1, 0, &
                      ps%qatom_s(0), ps%q_s(0), ps%m_s(0), ierr)  ! electron (species 0)
                CALL ckerr(' species_convert')

       WRITE (*,*) "Converted thermal electrons"
 
                DO i = 1, nspec_th ! Thermal ions
                    CALL    ps_species_convert(q_s(i), q_s(i), m_s(i), &
                                ps%qatom_s(i), ps%q_s(i), ps%m_s(i), ierr)  
                    CALL ckerr(' species_convert')
                END DO
 
       WRITE (*,*) "Converted thermal species"
           
            !-----------------------------------
            ! ICRF minority ion species
            !-----------------------------------
                 
                DO i = 1, ps%nspec_rfmin ! ICRF minority ions
                    CALL    ps_species_convert(q_rfmin(i), q_rfmin(i), m_rfmin(i), &
                                ps%qatom_rfmin(i), ps%q_rfmin(i), ps%m_rfmin(i), ierr)  
                    CALL ckerr(' species_convert')         
                END DO
                ps%kdens_rfmin = kdens_rfmin 
                ps%fracmin = fracmin

       WRITE (*,*) "Converted minority ion species"
            
            !-----------------------------------
            ! Label species and merge lists
            !-----------------------------------
        
                CALL ps_label_species(ierr)
                CALL ckerr('ps_label_species')
    
                CALL ps_merge_species_lists(ierr)
                CALL ckerr('ps_merge_species_lists')

          !  Print out combined species list...
    
      do ii=0, ps%nspec_all
         write(*, 1001) ii,ps%all_type(ii), &
              trim(ps%all_name(ii)),ps%q_all(ii),ps%m_all(ii)
      enddo
    1001 format(' Specie index & type: ',i2,1x,i2,1x, &
              '"',a,'" charge & mass: ',2(1pe12.5,1x))
              
    !--------------------------------------------------------------------------
    ! Initialize and allocate RF source arrayw
    !--------------------------------------------------------------------------
            
        ps%nicrf_src = nicrf_src               
        WRITE (*,*) 'model_epa: ps%nicrf_src     = ', ps%nicrf_src
       
        WRITE (*,*) 'model_epa: About to allocate ICRF source arrays' 
        CALL    ps_alloc_plasma_state(ierr)
        WRITE (*,*) 'model_epa: Allocated ICRF source arrays'
 
 
    !-----------------------------------
    !
    ! Insert initial plasma state data
    !
    !-----------------------------------
   
        ps%R_axis = r_axis
        ps%Z_axis = z_axis
        ps%R_min_box = r_min_box
        ps%R_max_box = r_max_box
        ps%Z_min_box = z_min_box
        ps%Z_max_box = z_max_box

          
    !---------------------------------------------------------------------------------
    !     
    !  Get data for profile and time evolution models from input_namelist_file
    !
    !---------------------------------------------------------------------------------

        OPEN (unit=21, file=TRIM(input_namelist_file), status='old', &
            action='read', iostat=istat, form='formatted')
            IF (istat /= 0 ) THEN
                CALL SWIM_error ('open', 'model_epa.f90',TRIM(input_namelist_file))
                ierr = istat
                stop 'cannot open input_namelist_file'

            END IF
        ierr = 0
        read(21, nml=evolving_model_data)
        CLOSE (21)
        WRITE (*, nml=evolving_model_data)
        !---------------------------------------------------------------------------------
        ! Put initial data in the thermal profile state arrays
        !---------------------------------------------------------------------------------
                    
            CALL rho_grid(ps%nrho,ps%rho)
            
            nzone = ps%nrho -1       
            ALLOCATE( zone_center(nzone), stat=istat )
            IF (istat /= 0 ) THEN
                CALL SWIM_error ('allocation', 'model_epa' , 'zone_center')
                ierr = istat
            END IF
        
            zone_center = ( ps%rho(1:nrho-1) + ps%rho(2:nrho) )/2.
                    
            ps%ns = 0.0
            CALL profgen( ne0, ne_edge, alpha_ne, zone_center, ps%ns(:,0) ) ! ne
            ps%ns(:,1) = ps%ns(:,0) ! thermal ion 1 density = ne, other ion densities are 0

            IF (nspec_th > 1) THEN
                WRITE (*,*) 'model_epa: INIT'
                WRITE (*,*) 'nspec_th = ', nspec_th, &
                & ' present model only has nonzero density in ion species 1'
            END IF
            
            CALL profgen( te0, te_edge, alpha_Te, zone_center, ps%Ts(:,0) )
            
            DO i = 1, nspec_th
                CALL profgen( ti0, ti_edge, alpha_Ti, zone_center, ps%Ts(:,i) )
            END DO
      
           
        !--------------------------------------------------------------------------
        ! RF input data
        !--------------------------------------------------------------------------

            IF (nicrf_src >= 0) THEN
                ps%icrf_src_name = icrf_src_name(1:nicrf_src)    
                ps%power_ic = P_icrf_0(1:nicrf_src)

                WRITE (*,*) 'model_epa: ps%power_ic  = ', ps%power_ic, ' src_name = ', &
                         trim(ps%icrf_src_name(1))
            END IF
          
    !-------------------------------------------------------------------------- 
    ! Load equilibrium data from model_eqdsk_file if one is provided
    !--------------------------------------------------------------------------

    IF (TRIM(model_eqdsk_file) .ne. ' ') THEN
        CALL ps_update_equilibrium(ierr, TRIM(model_eqdsk_file) )
        IF (ierr .ne. 0) STOP 'model_epa INIT: ps_update_equilibrium call failed'

        CALL ps_wr_geqdsk(ierr, cur_eqdsk_file)
        IF (ierr .ne. 0) STOP 'model_epa INIT: eqdsk write failed'
        
    END IF

    !-------------------------------------------------------------------------- 
    ! Store initial plasma state
    !--------------------------------------------------------------------------

        CALL    ps_STORE_PLASMA_STATE(ierr, cur_state_file)
        IF (ierr .ne. 0) THEN
            STOP 'model_epa INIT: STORE_PLASMA_STATE call failed'
        ELSE 
            WRITE (*,*) "model_epa: Stored initial Plasma State" 
        END IF
        
END IF  ! End of INIT function
        
!------------------------------------------------------------------------------------
!     
!  STEP function - Change state data and store plasma state
!
!------------------------------------------------------------------------------------
        
IF (TRIM(mode) == 'STEP') THEN    
        
    WRITE (*,*)
    WRITE (*,*) 'model_epa STEP'
    
    !------------------------------------------------------------------------------------
    !     
    !  Get current plasma state 
    !
    !------------------------------------------------------------------------------------
                    
        call    ps_get_plasma_state(ierr, trim(cur_state_file))
        if(ierr .ne. 0) stop 'call failed to ps_get_plasma_state to load profiles '
          
    !---------------------------------------------------------------------------------
    !     
    !  Get data for profile and time evolution models from input_namelist_file
    !
    !---------------------------------------------------------------------------------

        OPEN (unit=21, file=TRIM(input_namelist_file), status='old', &
            action='read', iostat=istat, form='formatted')
            IF (istat /= 0 ) THEN
                CALL SWIM_error ('open', 'model_epa.f90',TRIM(input_namelist_file))
                ierr = istat
                stop 'cannot open input_namelist_file'

            END IF
        ierr = 0
        read(21, nml=evolving_model_data)
        CLOSE (21)
        
        
    !--------------------------------------------------------------------------
    !
    !   Evolve state data
    !
    !--------------------------------------------------------------------------
        
        READ (time_stamp, '(f12.6)' ) t
        WRITE (*,*) 'model_epa STEP: time stamp = ', t, '   ps%t1 = ', ps%t1

        !---------------------------------------------------------------------------------
        ! Evolve data in the thermal profile state arrays
        !---------------------------------------------------------------------------------

            nzone = ps%nrho -1       
            ALLOCATE( zone_center(nzone), stat=istat )
            IF (istat /= 0 ) THEN
                CALL SWIM_error ('allocation', 'model_epa' , 'zone_center')
                ierr = istat
            END IF
        
            zone_center = ( ps%rho(1:ps%nrho-1) + ps%rho(2:ps%nrho) )/2.
       WRITE (*,*) 'ps%nrho = ', ps% nrho, '     zone_center = ', zone_center
         
            ne0 =ps%ns(1,0) *(1.0_rspec + delta_ne)
            ne_edge = ps%ns(nzone,0)
            te0 = ps%Ts(1,0)*(1.0_rspec + delta_Te)
            te_edge = ps%Ts(nzone,0)
            ti0 = ps%Ts(1,1)*(1.0_rspec + delta_Ti)
            ti_edge = ps%ns(nzone,1)
            
       WRITE(*,*) 'model_epa STEP: about to generate ne profile' 
        
            CALL profgen( ne0, ne_edge, alpha_ne, zone_center, ps%ns(:,0) ) ! ne
            ps%ns(:,1) = ps%ns(:,0) ! thermal ion 1 density = ne

            IF (ps%nspec_th > 1) THEN
                ps%ns(:, 2:ps%nspec_th) = 0.0_rspec ! other ion densities are 0
                WRITE (*,*) 'model_epa: STEP'
                WRITE (*,*) 'ps%nspec_th = ', ps%nspec_th, &
                & ' present model only has nonzero density in ion species 1'
            END IF
           
       WRITE(*,*) 'model_epa STEP: about to generate Te profile'

            CALL profgen( te0, te_edge, alpha_Te, zone_center, ps%Ts(:,0) )
            
            DO i = 1, ps%nspec_th
                CALL profgen( ti0, ti_edge, alpha_Ti, zone_center, ps%Ts(:,i) )
            END DO
      
      
    !--------------------------------------------------------------------------
    ! Store the data in plasma_state file
    !--------------------------------------------------------------------------
                
        CALL    ps_STORE_PLASMA_STATE(ierr, cur_state_file)
        WRITE (*,*) "model_epa STEP: Stored Plasma State"    

    !--------------------------------------------------------------------------
    !    Change total ICRF power for next time step
    !--------------------------------------------------------------------------
        
        IF (ps%t1 >= t_icrf_1) THEN
            ps%power_ic =  P_icrf_1(1:nicrf_src)
            WRITE(*,*) 'change_n_T_EPA (STEP): ps%power_ic = ', ps%power_ic
        END IF
    !--------------------------------------------------------------------------
    ! Copy current plasma state to next plasma state and store
    !--------------------------------------------------------------------------

      CALL    ps_STORE_PLASMA_STATE(ierr, trim(next_state_file) )
        if(ierr .ne. 0) stop 'model_epa STEP: failed to store_plasma_state next'

        WRITE (*,*) "model_epa STEP: Stored  next Plasma State"    
    
        
        CALL sleep(5)    ! Wait for elvis to pick up
END IF                
            
!--------------------------------------------------------------------------
!
!   Finalize function
!
!--------------------------------------------------------------------------

IF (TRIM(mode) .eq. 'FINALIZE') THEN    
                                                                                                
     WRITE (*,*) "model_epa: Finalize called"
     
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
    SUBROUTINE profgen(f0,f1, alpha, rho,f)

      !  quick & dirty parabolic to a power profile generator

      REAL(KIND=rspec), intent(in) :: f0,f1  ! core and edge values
      REAL(KIND=rspec), intent(in) :: alpha  ! power of parabolic
      REAL(KIND=rspec), dimension(:), intent(in) :: rho  ! normalized sqrt tor. flux
      REAL(KIND=rspec), dimension(:), intent(out) :: f  ! formula output
      
      f = 0.0_rspec
      WHERE (1.0_rspec-rho**2 > 0.0_rspec) f = (f0-f1) *(1.0_rspec-rho**2)**alpha + f1

    END SUBROUTINE profgen

    !------------------------------------
    SUBROUTINE ckerr(sbrtn)
      character*(*), intent(in) :: sbrtn

      IF(ierr.NE.0) then
         write(6,*) ' ?plasma_state_test: error in call: '//trim(sbrtn)
         stop
      ENDIF
    END SUBROUTINE ckerr


END PROGRAM model_epa
