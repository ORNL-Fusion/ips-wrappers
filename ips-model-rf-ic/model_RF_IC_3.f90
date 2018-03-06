PROGRAM model_RF_IC

!--------------------------------------------------------------------------
! Version 3.0 11/7/2012 (Batchelor) Added capability to initialize ICRF antenna parameters
!             from machine description and shot configuration namelists using the plasma
!             state functions ps_mdescr_read() and ps_sconfig_read().

! version 2.0 8/31/2009 (Batchelor)  Added capability to initialize from an input plasma state
!             file. Added a model routine that gives off-axis peaked profiles: Lorentz_Linear.
!             Generally cleaned up.
! version 0.0 9/11/2008 (Batchelor)

!--------------------------------------------------------------------------
!
!   Simple mock RF_IC code for testing.
!
!       The code requires 4 command-line arguments
!       1) path to the current plasma state file
!       2) path to the current plasma eqdsk file
!       3) action mode, i.e. one of: "INIT", "STEP", or "FINALIZE"
!       4) time stamp the time set by the driver component to which the simulation is 
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
!    a normal initialization sequence in which the EPA has allocated the arrays 
!    associated with the icrf sources: 
!
!       nicrf_src, icrf_src_name(nicrf_src), power_ic(nicrf_src), nspec_rfmin, 
!       and the RFMIN species arrays. 
!
!    This model component then specifies nrho_icrf and proceeds to allocate and initialize
!    all the variables dimensioned by nrho_icrf:
!
!      nrho_icrf, rho_icrf(nrho_icrf), picrf_srcs(nrho_icrf-1, nicrf_src,0:nspc_alla)
!      picrt_totals(nrho_icrf-1, 0:nspc_alla), picth(nrho_icrf-1), curich(nrho_icrf-1)
!      eperp_mini(nrho_icrf-1, nspec_rfmin)
!
!   The data to be loaded is derived from 3 namelist files, two of which correspond to
!   data sections of the plasma state: 1) Machine description, and 2) Shot configuration.
!   The third file is specific to this model component and provides initialization data
!   that in a real physics component would come in from the standard TORIC or AORSA
!   namelist file.  It also has data for the evolving ICRF models.
!
!   1) Machine description namelist file: mdescr_namelist.dat.  mdescr files can be
!   written by plasma state function ps_mdescr_write().  Sample mdescr files are produced
!   by the plasma_state_test routine.  The ICRF relevant items in namelist mdescr are:
!        icrf_src_name  = number & name of ICRF sources
!        ant_model = antenna model filenames (1 per antenna source)
!        nrz_antgeo = number of (R,Z) points, antenna geometries
!        max_nrz_antgeo = Maximum size of variable length enumeration: NRZ_ANTGEO
!        R_antgeo = antenna geo: R pts 
!        Z_antgeo = antenna geo: Z pts 
!        dx_fshield = distance, antenna to Faraday shield
!
!   2) Shot configuration namelist file: sconfig_namelist.dat.  sconfig files can be
!   written by plasma state function ps_sconfig_write(). Sample sconfig files are produced
!   by the plasma_state_test routine.  The ICRF relevant items in namelist sconfig are:
!        RFMIN = ICRF minority species
!        freq_ic = frequency on each ICRF source
!        n_straps = number of straps in the antenna
!        max_n_straps = Maximum size of variable length enumeration: N_STRAPS
!        R_ant = major radius of antenna
!        Z_ant = height of antenna
!        Z_mid_ant = center of antenna relative to TF coil midplane
!        num_nphi_vac = number of non-zero n_phi values
!        max_num_nphi_vac = Maximum size of variable length enumeration: NUM_NPHI_VAC
!        Real_ant_coef = real part of Fourier Coef
!        Imag_ant_coef = imag part of Fourier Coef
!        num_nphi = number of non-zero n_phi values
!        max_num_nphi = Maximum size of variable length enumeration: NUM_NPHI
!        nphi = n_phi wave spectrum from antenna
!        wt_nphi = vacuum spectrum n_phi weight
!
!   3) input namelist file model_RF_IC_input.nml.
!   This file must contain two namelists:
!
!   /static_state_data/ -> data used to initialize the RF_IC part of the plasma state 
!        input_state_file = ' ' or name of the initialization plasma state file if used
!        freq_ic = freqency of each ICRF source
!        nrho_icrf = grid dimension of rho_icrf
!
!   /evolving_model_data/ -> Contains data that defines the ad hoc models used to give time 
!                            varying numbers for the RF_IC data (e.g names of different models
!                            and parameters to be used in the models, like profiles shapes or
!                            parameters to define time variation
!
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

    CHARACTER (len=256) :: cur_state_file, cur_eqdsk_file
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

    !--------------------------------------------------------------------
    ! Machine description data (copied from plasma_state_kernel/ps_mdescr_namelist_read.f90
    !--------------------------------------------------------------------

          CHARACTER*32, dimension(ps_max_static_2d) :: icrf_src_name
          CHARACTER*256, dimension(ps_max_static_2d) :: ant_model
          INTEGER, dimension(ps_max_static_2d) :: nrz_antgeo
          INTEGER :: max_nrz_antgeo
          REAL(KIND=rspec), dimension(ps_max_static_2d,ps_max_static_2d) :: R_antgeo
          REAL(KIND=rspec), dimension(ps_max_static_2d,ps_max_static_2d) :: Z_antgeo
          REAL(KIND=rspec), dimension(ps_max_static_2d) :: dx_fshield

    !--------------------------------------------------------------------
    ! Shot configuration data (copied from plasma_state/ps_sconfig_namelist_read.f90)
    !--------------------------------------------------------------------

          REAL(KIND=rspec), dimension(ps_max_static_2d) :: freq_ic
          INTEGER, dimension(ps_max_static_2d) :: n_straps
          INTEGER :: max_n_straps
          REAL(KIND=rspec), dimension(ps_max_static_2d) :: R_ant
          REAL(KIND=rspec), dimension(ps_max_static_2d) :: Z_ant
          REAL(KIND=rspec), dimension(ps_max_static_2d) :: Z_mid_ant
          INTEGER, dimension(ps_max_static_2d) :: num_nphi_vac
          INTEGER :: max_num_nphi_vac
          REAL(KIND=rspec), dimension(ps_max_static_2d,ps_max_static_2d) :: Real_ant_coef
          REAL(KIND=rspec), dimension(ps_max_static_2d,ps_max_static_2d) :: Imag_ant_coef
          INTEGER, dimension(ps_max_static_2d) :: num_nphi
          INTEGER :: max_num_nphi
          INTEGER, dimension(ps_max_static_2d,ps_max_static_2d) :: nphi
          REAL(KIND=rspec), dimension(ps_max_static_2d,ps_max_static_2d) :: wt_nphi
        
    !--------------------------------------------------------------------------
    !   Internal data
    !--------------------------------------------------------------------------
    CHARACTER (len=256) :: input_state_file
    INTEGER :: nzone
    REAL(KIND=rspec), ALLOCATABLE :: zone_center(:)

    INTEGER :: nicrf_src, nrho_icrf
    REAL(KIND=rspec) :: power_ic(ps_max_static_2d)

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
      if (iarg .ne. 4) then
         print*, 'model_RF_IC usage: '
         print*, ' command line args = cur_state_file cur_eqdsk_file mode time_stamp'
         call exit(1)
      end if
      
      call getarg(1,cur_state_file)
      call getarg(2,cur_eqdsk_file)
      call getarg(3,mode)
      call getarg(4,time_stamp)
      
     WRITE (*,*)
     WRITE (*,*) 'model_RF_IC'      
     print*, 'cur_state_file = ', trim(cur_state_file)
     print*, 'cur_eqdsk_file = ', trim(cur_eqdsk_file)
     print*, 'mode = ', trim(mode)
     print*, 'time_stamp = ', trim(time_stamp)
      
    !------------------------------------------------------------------------------------
    !  Get current plasma state 
    !------------------------------------------------------------------------------------
            
    call ps_get_plasma_state(ierr, trim(cur_state_file))
    if (ierr .ne. 0) then
       print*, 'model_RF_IC:failed to get_plasma_state'
       call exit(1)
    end if

    !--------------------------------------------------------------------------
    !    Open input namelist file
    !--------------------------------------------------------------------------

    OPEN (unit=21, file = 'model_RF_IC_input.nml', status = 'old',   &
         form = 'formatted', iostat = ierr)
    IF (ierr .ne. 0) then
        print*, 'cannot open RF_IC_model.nml'
        call exit(1)
    END IF
    
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

    IF (len(TRIM(input_state_file)) .ne. 0) THEN
    
        print*, 'Getting RF_IC source arrays from plasma state file ', TRIM(input_state_file)
          
        !---------------------------------------------------------------------------------
        !  Get aux plasma state to copy RF_IC data from. It comes in from input_state_file
        !---------------------------------------------------------------------------------
            
          CALL ps_get_plasma_state(ierr, TRIM(input_state_file), aux)
          if (ierr .ne. 0) then
              print*, 'call failed to ps_get_plasma_state for aux state'
              call exit(1)
          end if
    
        !--------------------------------------------------------------------------
        !   Copy data in RF_IC sections of 'aux' state to current (i.e. 'ps') state. Note
        !   that this retains the rho_icrf grid and initial profile data from the input 
        !   plasma state.
        !--------------------------------------------------------------------------
        

        CALL ps_store_plasma_state(ierr, 'initial_ps.cdf', ps)
        CALL ps_store_plasma_state(ierr, 'aux_ps.cdf', aux)
        
        !--------------------------------------------------------------------------
        ! Copy Zimp1 and Zimp2 from ps to aux to avoid these being clobbered in the
        ! plasma state copy from aux to ps
        !--------------------------------------------------------------------------
        
        aux%Zimp1 = ps%Zimp1
        aux%Zimp2 = ps%Zimp2
 
         CALL ps_cclist_remove("*", cclist, iwarn)
         CALL ps_cclist_add("IC" ,cclist, iwarn)
         CALL PS_COPY_DIMS(aux, ps, 0 ,ierr, cclist, 0)
         if (ierr .ne. 0) then
             print*, 'call failed to PS_COPY_DIMS for aux state to ps state'
             call exit(1)
         end if
         CALL PS_COPY_PLASMA_STATE(aux, ps, ierr, cclist, 1, 0)
         if (ierr .ne. 0) then
             print*, 'call failed to PS_COPY_PLASMA_STATE for aux state to ps state'
             call exit(1)
         end if
            
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

    ELSE
    
    !--------------------------------------------------------------------------
    !   1.2 Case: initialization data taken from namelists
    !--------------------------------------------------------------------------

        !---------------------------------------------------------------------------------
        !  Check if nicrf_src dimensioned arrays are already allocated in current plasma
        !  state.  EPA component should have done this in its init.
        !---------------------------------------------------------------------------------
                     
        IF ( allocated(ps%icrf_src_name) ) THEN
             print*, 'Checking ICRF data in current Plasma State'  
             
             !ps%freq_ic = freq_ic(1:ps%nicrf_src) ! EPA doesn't set source frequency
                
             print*,   '   number of icrf sources = ', ps%nicrf_src
             print*,   '   names of sources = ', ps%icrf_src_name
             print*,   '   icrf frequencies = ', ps%freq_ic
             print*,   '   icrf source powers = ', ps%power_ic
             print*,   '   number of minority species = ', ps%nspec_rfmin
             print*,   '   names of minority species = ', ps%RFMIN_name
             print*,   '   eperp_icrf allocated = ', ALLOCATED(ps%eperp_mini)
             print*,   '   epll_icrf allocated = ', ALLOCATED(ps%epll_mini)
    
        !-------------------------------------------------------------------------- 
        ! Otherwise quit
        !--------------------------------------------------------------------------
            
        ELSE          
             print*, 'RF_IC sources NOT allocated in initial Plasma State'
             print*, 'Can''t be done in this component without stepping on species lists'
             print*, 'model_RF_IC_3.f90: I give up'
             call exit(1)
                        
        END IF


        !-------------------------------------------------------------------------- 
        ! Get ICRF machine configuration data from machine description
        ! namelist file.  Put into aux plasma state
        !--------------------------------------------------------------------------

        CALL ps_store_plasma_state(ierr, 'initial_aux_ps.cdf', aux)

        call ps_mdescr_namelist_read(.False., 'mdescr_namelist.dat', ' ',  ' ', aux, ierr)
        IF (ierr .ne. 0) THEN
            print*, 'Could not get namelist mdescr'
            call exit(1)
        END IF
        
        ! Check that nicrf_src read from namelist is consistent with that from current
        ! plasma state.
        
        IF (ps%nicrf_src .ne. aux%nicrf_src) THEN
            print*, 'nicrf_src read from namelist is inconsistent with that from current&
                & plasma state'
            call exit(1)
        END IF      
             print*,   '   number of icrf sources = ', ps%nicrf_src

        !-------------------------------------------------------------------------- 
        ! Get ICRF shot configuration data from shot configuration description
        ! namelist file.  Put into aux plasma state
        !--------------------------------------------------------------------------

        call ps_sconfig_namelist_read(.False., 'sconfig_namelist.dat', ' ', aux, ierr)
        IF (ierr .ne. 0) THEN
            print*, 'Could not get namelist sconfig'
            call exit(1)
        END IF

             print*,   'aux%icrf_src_name = ', aux%icrf_src_name      
             print*,   'aux%ant_model = ', aux%ant_model          
             print*,   'aux%nrz_antgeo = ', aux%nrz_antgeo        
             print*,   'aux%freq_ic = ', aux%freq_ic          
             print*,   'aux%max_num_nphi = ', aux%max_num_nphi        
             print*,   'aux%nphi = ', aux%nphi        
             print*,   'aux%Real_ant_coef = ', aux%Real_ant_coef          

        CALL ps_store_plasma_state(ierr, 'initial_ps.cdf', ps)
        CALL ps_store_plasma_state(ierr, 'aux_ps.cdf', aux)
        
        !--------------------------------------------------------------------------
        ! Copy Zimp1 and Zimp2 from ps to aux to avoid these being clobbered in the
        ! plasma state copy from aux to ps
        !--------------------------------------------------------------------------
        
        aux%Zimp1 = ps%Zimp1
        aux%Zimp2 = ps%Zimp2
        
        !--------------------------------------------------------------------------
        ! Copy data in RF_IC sections of 'aux' state to current (i.e. 'ps') state. 
        !--------------------------------------------------------------------------
         
         CALL ps_cclist_remove("*", cclist, iwarn)
         CALL ps_cclist_add("IC" ,cclist, iwarn)
         
         CALL PS_COPY_DIMS(aux, ps, 0 ,ierr, cclist, 0)
         if (ierr .ne. 0) then
             print*, 'call failed to PS_COPY_DIMS for aux state to ps state'
             call exit(1)
         end if
         CALL PS_COPY_PLASMA_STATE(aux, ps, ierr, cclist, 1, 0)
         if (ierr .ne. 0) then
             print*, 'call failed to PS_COPY_PLASMA_STATE for aux state to ps state'
             call exit(1)
         end if

             print*,   'ps%icrf_src_name = ', ps%icrf_src_name        
             print*,   'ps%ant_model = ', ps%ant_model        
             print*,   'ps%nrz_antgeo = ', ps%nrz_antgeo          
             print*,   'ps%freq_ic = ', ps%freq_ic        
             print*,   'ps%max_num_nphi = ', ps%max_num_nphi          
             print*,   'ps%nphi = ', ps%nphi          
             print*,   'ps%Real_ant_coef = ', ps%Real_ant_coef        
            
        !-------------------------------------------------------------------------- 
        ! Allocate rho_icrf dimensioned profile arrays in plasma state.
        ! nrho_icrf is in /static_state_data/ namelist.  For a real RF component
        ! it would be in standard TORIC or AORSA namelist.
        !--------------------------------------------------------------------------
                 
        IF ( allocated(ps%rho_icrf) ) THEN
        
             print*, "model_RF_IC: rho_icrf(:) already allocated in plasma state &
             &        Error - I'm supposed to do that"
             call exit(1)
        ELSE
            ps%nrho_icrf = nrho_icrf  
            CALL ps_alloc_plasma_state(ierr)
            WRITE (*,*) 'model_RF_IC: Allocated icrf profiles in plasma state'      

            CALL rho_grid(ps%nrho_icrf,ps%rho_icrf)
            write (*,*) 'ps%rho_icrf = ', ps%rho_icrf

            !--------------------------------------------------------------------------
            ! Initialize RF_IC output data to zero if a time dependent model is provided.
            ! If the model is 'const' just keep the icrf data in the aux state file
            !--------------------------------------------------------------------------

            IF (TRIM(RF_IC_profile_model_name) /= 'const') THEN         
                ps%picrf_srcs = 0.        
                ps%picrf_totals = 0.        
                ps%picth = 0.
                ps%curich = 0.
            END IF
            
        END IF

    END IF  ! End of initialization from namelists
            
    
    !-------------------------------------------------------------------------- 
    !1.4  Store initial plasma state for RF_IC
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)
        IF (ierr .ne. 0) THEN
        print*,  'model_RF_IC INIT: PS_STORE_PLASMA_STATE failed'
        call exit(1)
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
            call exit(1)
        END IF  
        nrho_icrf = ps%nrho_icrf
        zone_center = ( ps%rho_icrf(1:nrho_icrf-1) + ps%rho_icrf(2:nrho_icrf) )/2.

        ps%picrf_totals = 0.

        FP_min = 1.0 - FP_th_e_icrf - FP_th_i_icrf   ! minority gets what's left
        
        DO j = 1, ps%nicrf_src  ! Loop over ICRF sources
            
            WRITE (*,*) 'power_ic(', j, ') = ',  ps%power_ic(j)
            
            ! direct ICRF power to thermal electrons
            CALL Lorentz_Linear(rho_max_P_th_e, w_P_th_e, f0_P_th_e, f1_P_th_e, &
            zone_center, ps%picrf_srcs(:, j, 0))
            ps%picrf_srcs(:, j, 0) = FP_th_e_icrf*ps%power_ic(j)*ps%picrf_srcs(:, j, 0)/ &
                                     SUM(ps%picrf_srcs(:, j, 0) )
            ps%picrf_totals(:, 0) = ps%picrf_totals(:,0) + ps%picrf_srcs(:, j, 0)

            ! WRITE (*,*) 'ps%picrf_srcs(:,', j, '0) = ', ps%picrf_srcs(:, j, 0)
            ! WRITE (*,*) 'ps%picrf_totals(:, 0) = ', ps%picrf_srcs(:, j, 0)
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
                    WRITE (*,*) 'SUM(ps%picrf_totals(:,' , i,')) = ', SUM(ps%picrf_totals(:, i) )
                 END DO

            END IF
              
        END DO

        ps%picth(:) = 0. 
        DO i=1, MIN(ps%nspec_th, 2)
            ps%picth(:) = ps%picth(:) + ps%picrf_totals(:,i)
        END DO
        WRITE (*,*) 'SUM(ps%picth(:)) = ', SUM(ps%picth(:) )

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

    CALL PS_WRITE_UPDATE_FILE('RF_IC_'//cur_state_file, ierr)

    WRITE (*,*) "model_RF_IC_3: Stored Partial RF Plasma State"    

    CALL sleep(3)   ! Wait for elvis to pick up
    
END IF ! End of STEP function

!--------------------------------------------------------------------------
!
!  Section 3: Finalize function
!
!--------------------------------------------------------------------------

IF (TRIM(mode) .eq. 'FINALIZE') THEN
                                                                                                
     WRITE (*,*) "model_RF_IC: Finalize called"
     
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


!--------------------------------------------------------------------------
!
!   A set of simple models to generate radial profiles based on adjustable input
!   parameters.  Input rho is a vector of values assumed to run from 0 to 1.
!
!   Model 1: Power_Parabolic_Offset(alpha, h, rho, f) -> Parabolic to a power
!   plus an offset. 
!       f = A (1 - rho**2)**alpha + B where A and B are chosen so that
!           Integral rho*f(rho) from 0 to 1 is unity, and f(1)/f(0) = h
!
!   Model 2: Lorentz_Linear(rho_max, w, f0, f1, rho, f)
!       f = Lorentzian(centered at rho = rho_max, width = w) + f0 + (f1 - f0)*rho
!
!   Model 3: Lorentz_Linear_norm(rho_max, w, f0, f1, rho, f)
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

    subroutine ps_namrd_ilist_chk(listname,nmax,list,nfound,iout,ierr)
    
      ! determine length of list by number of non-blank elements in the list
      ! once a blank element is detected, all subsequent elements must also
      ! be blank, or, an error flag is set.
    
      implicit NONE
    
      !-----------------
      ! arguments:
    
      character*(*), intent(in) :: listname   ! name of list (for error message)
      integer, intent(in) :: nmax             ! size of list array
      character*(*), intent(in) :: list(nmax) ! the list itself
    
      integer, intent(out) :: nfound          ! address of last non-blank element
    
      integer, intent(in) :: iout             ! I/O unit for error messages
      integer, intent(out) :: ierr            ! exit status (0=OK)
    
      !------------------
      ! local:
    
      integer :: ii,iblank
      !------------------
    
      nfound=0
      ierr=0
    
      iblank=0
      do ii=1,nmax
         if(list(ii).eq.' ') then
            if(iblank.eq.0) iblank=ii
         else
            if(iblank.gt.0) then
               write(iout,*) ' ?ps_namrd_ilist_chk: in list "'//trim(listname)//'":'
               write(iout,*) '  non-blank element #',ii, &
                    ' follows blank element #',iblank
               ierr=1
               exit
            endif
         endif
      enddo
    
      if(ierr.eq.0) then
         nfound = iblank - 1
      endif
    
    end subroutine ps_namrd_ilist_chk


END PROGRAM model_RF_IC
