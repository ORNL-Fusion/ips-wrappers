PROGRAM model_epa_ps_file_init_no_RF

! no_RF version 1.0 11/29/2012 (Batchelor)
! This is adaped from model_epa_ps_file_init.f90 version 3.3.  As before it generates the initial
! plasma state file from an input_state_file, but now it does not copy the machine description
! data for the IC component from the input state file.  Instead it reads a mdescr namelist file
! for that data.  This enables one to change the number of RF spources from what is in the original
! input_state_file.  The only application I have in mind is for testing the multi-toroidal-mode 
! ICRF component, but maybe others will come up.

!--------------------------------------------------------------------------
!
! Details:
!
!	1) Machine description namelist file: mdescr_namelist.dat.  mdescr files can be
!	written by plasma state function ps_mdescr_write().  Sample mdescr files are produced
!	by the plasma_state_test routine.  The ICRF relevant items in namelist mdescr are:
!        icrf_src_name  = number & name of ICRF sources
!        ant_model = antenna model filenames (1 per antenna source)
!        nrz_antgeo = number of (R,Z) points, antenna geometries
!        max_nrz_antgeo = Maximum size of variable length enumeration: NRZ_ANTGEO
!        R_antgeo = antenna geo: R pts 
!        Z_antgeo = antenna geo: Z pts 
!        dx_fshield = distance, antenna to Faraday shield
!
!	2)	The code takes one command-line argument, the name of the current plasma state
!		to be initialized.  The names of the namelist files are assumed to be
!		mdescr_namelist.dat and sconfig_namelist.dat.
!
! One other point: The EPA component is responsible for setting the number of ICRF
! sources, nicrf_src and allocating the arrays dimensioned by it.  This is because EPA
! controls the power on the different sources.  Note that nicrf_src is determined in
! the mdescr namelist.  Although it appears only implicitly as the number
! of elements in the namelist that are dimensioned nicrf_src.  For example icrf_src_name,
! ant_model and dx_fshield.
!
! For that matter many of the arrays in the sconfig namelist are also dimensioned
! nicrf_src, for example freq_ic.  So the dimensons in the two namelist files must be
! consistent.
!
!--------------------------------------------------------------------------

! Version 3.3 3/9/2011 (Batchelor)
! Eliminated the next_state_file. Now the source power levels as of the end of the time
! step are put into the full and partial plasmas states that get written.  This is to be
! consistent with the new generic_driver.py that no longer copies next_state to current_state
! at the beginning of each time step.

! version 3.2 1/11/2011 (Batchelor)

! Re-added the behavior that "input_eqdsk_file" overwrites the equilibrium
! data from the "input_state_file".  This works by attempting to do 
! ps_update_equilibrium(ierr, TRIM(input_eqdsk_file) ).  If there is no "input_eqdsk_file"
! present the update will fail (not fatally I hope) so the equilibrium data from
! "input_state_file" stays.  It could be that we will have input plasma state files that have 
! the equilibrium data we want to use and have no eqdsk files to go with them.  In that case we
! will want to go back to writing the eqdsk file from the plasma state.  I worry about that when! it
! happens.

! version 3.1 10/18/2010 (Batchelor)

! This version writes both a full plasma state and a partial plasma state.  This way the
! Python component in the step function can use either services.update_plasma_state() to 
! write a full plasma state or services.merge_current_plasma_state() to write a partial
! plasma state depending on how the component is being used.

! version 3.0 10/13/2010 (Batchelor)

! This version eliminates manipulation of the equilibrium data from an input eqdsk file.  
! It also eliminates input_eqdsk_file and input_state_file from the /static_state_data/
! namelist.  The input_eqdsk_file is not used and input_state_file is now just assumed to
! be generically named "input_state_file".  Now the python component copies the input eqdsk
! file to current_eqdsk_file and copies the input state file to the generic name.  
! In the [model_EPA] section of the simulation config file the paths
! to the input plasma state file and input eqdsk files should be specified as input files
! in the simulation config file.  Also model_EPA config variables "INPUT_PLASMA_STATE_FILE"
! and "INPUT_EQDSK_FILE"  should be defined  with the names of these input files.
!
! Note: The plasma state function pr_wr_geqdsk function writes slightly
! different eqdsk files from those written by TSC for what is supposed to be the same
! plasma state data.  The consensus is that the TSC versions are more trustworthy.
! Hence the decision, in the 'init', to eliminate generating the current eqdsk file from the plasma 
! state.  The Python component script just copies INPUT_PLASMA_STATE_FILE to current plasma state.

 


! version 2.0 9/15/2009 (Batchelor)

!   A EPA component that puts specified data into the plasma state for testing purposes and to
!   provide a simple driver for other components.  It is intended to be incremental in that
!   initially it only loads data needed for the RF_IC and neutral beam components.  Other 
!   data can be added later.
!
!   This code is driven by the model_epa_ps_init.py component script
!
!   This version gets the EPA state data out of existing plasma state file 'input_state_file'.
!   The file path for this input state is specified in the /static_state_data/ namelist of the
!   'model_epa_ps_init.nml' namlist file (See below)
!
! Details:
!   The data to be loaded is derived from input namelist file model_EPA_input.nml.
!   This file must contain two namelists:
!
! Note (10/14/10): /static_state_data/ is now empty since the input state file name is now generic
!   /static_state_data/ -> Contains plasma state data that goes directly into the state, or
!						   data used to initialize the EPA part of the plasma state 
!
!   /evolving_model_data/ -> Contains data that defines the ad hoc models used to give time 
!                            varying numbers for the EPA data (e.g names of different models
!							 and parameters to be used in the models, like profiles shapes or
!							 parameters to define time variation
!
! So that the models for time dependence can have some memory a third namelist: 
! /internal_state_data/ is written to file 'internal_state_data.nml'
!
! Data owned by EPA for which time dependent models are provided:
!	
! 	MACHINE_DESCRIPTION data:  none
!	SHOT_CONFIGURATION DATA:  none
!	SIMULATION_INIT: none
!	STATE_DATA:
!		vsur			! toroidal voltage at surface
!	STATE_PROFILES:
!		curr_bootstrap	! Neoclassical bootstrap current
!		curr_ohmic		! Ohmic current
!		pohme			! Ohmic heating of thermal electrons
!		ns				! thermal species density ns(nrho-1, 0:nspc_th)
!		Ts				! thermal species temperature ns(nrho-1, 0:nspc_th)
!		vsur			! toroidal voltage at surface
!		Zeff			! Total Z effective
!		
! Data owned by EQ for which time dependent models are provided:
!	
! 	MACHINE_DESCRIPTION data:  none
!	SHOT_CONFIGURATION DATA:  none
!	SIMULATION_INIT: none
!	STATE_DATA: none
!	STATE_PROFILES:
!		curt			! integrated toroidal current inside rho
!		elong			! elongation (b/a)
!		triang			! triangularity (symmetrized)
!		P_eq			! equilibrium scalar pressure
!		q_eq			! equilibrium q profile
!
! Data owned by IC that is initialized by EPA:
!	
! 	MACHINE_DESCRIPTION data: 
!		nicrf_src		! number of icrf sources
!		icrf_src_name	! names of icrf sources
!	SHOT_CONFIGURATION DATA:
!		nspec_rfmin		! number of icrf minority species
!       m_rfmin, qatom_rfmin, q_rfmin, rfmin_name ! specification of minority species
!	SIMULATION_INIT: 
!		fracmin			! fraction of electron density in minorities (if kdens_rfmin = 'fraction')
!		kdens_rfmin		! specifies method to determine minority densities
!
! Data owned by IC for which time dependent models are provided:
!	STATE_DATA:
!		power_ic		! power on each icrf source
!	STATE_PROFILES:
!		nmini			! minority density profiles (maybe: IC sets nrho_icrf and rho_icrf)
!		
! Data owned by NBI that is initialized by EPA:
!	
! 	MACHINE_DESCRIPTION data: 
!		nbeam		! number of beams
!		nbi_src_name	! names of beam sources
!	SHOT_CONFIGURATION DATA:
!		nbion		! names of ion species injected by each beam
!	SIMULATION_INIT: 
!		nspec_beam		! number of beam ion species
!       m_beam, qatom_beam, q_beam, rfmin_beam ! specification of beam species
!
! Data owned by NBI for which time dependent models are provided:
!	STATE_DATA:
!		power_nbi		! power on each beam
!	STATE_PROFILES: none

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

    INTEGER :: ierr = 0, iwarn = 0
    INTEGER :: istat, iarg
    INTEGER :: i, ii, iout
    
    INTEGER :: system
    INTEGER :: cclist(ps_ccount)   ! component activation list 

    CHARACTER (len=256) :: cur_state_file, cur_eqdsk_file
    CHARACTER (len = SWIM_filename_length) :: input_state_file = "input_state_file"
    CHARACTER (len = SWIM_filename_length) :: input_eqdsk_file = "input_eqdsk_file"
    CHARACTER (len=32) :: mode
    CHARACTER (len=32) :: time_stamp
        

!--------------------------------------------------------------------------
!
!   Data for evolving profile models
!
!--------------------------------------------------------------------------

    CHARACTER(len = SWIM_filename_length) :: input_namelist_file = 'model_epa_input.nml'
    CHARACTER(len = SWIM_filename_length) :: internal_state_file = 'internal_state_data.nml'
    INTEGER :: n_step_calls = 0
    
    CHARACTER (len = SWIM_name_length) :: Te_evolution_model 
    CHARACTER (len = SWIM_name_length) :: Ti_evolution_model 
    CHARACTER (len = SWIM_name_length) :: n_evolution_model 
    CHARACTER (len = SWIM_name_length) :: curr_bootstrap_evolution_model 
    CHARACTER (len = SWIM_name_length) :: curr_ohmic_evolution_model 
    CHARACTER (len = SWIM_name_length) :: pohme_evolution_model 
    CHARACTER (len = SWIM_name_length) :: vsur_evolution_model 
    CHARACTER (len = SWIM_name_length) :: Zeff_evolution_model 
    CHARACTER (len = SWIM_name_length) :: curt_evolution_model 
    CHARACTER (len = SWIM_name_length) :: P_eq_evolution_model 
    CHARACTER (len = SWIM_name_length) :: q_eq_evolution_model 
    CHARACTER (len = SWIM_name_length) :: elong_evolution_model 
    CHARACTER (len = SWIM_name_length) :: triang_evolution_model 
    CHARACTER (len = SWIM_name_length) :: power_ic_evolution_model 
    CHARACTER (len = SWIM_name_length) :: power_nbi_evolution_model 

    REAL(KIND=rspec) :: delta_Te, f_Te, x_Te = 1
    REAL(KIND=rspec) :: delta_Ti, f_Ti, x_Ti = 1
    REAL(KIND=rspec) :: delta_n, f_n, x_n = 1
    REAL(KIND=rspec) :: delta_curr_bootstrap, f_curr_bootstrap, x_curr_bootstrap  = 1
    REAL(KIND=rspec) :: delta_curr_ohmic, f_curr_ohmic, x_curr_ohmic = 1
    REAL(KIND=rspec) :: delta_pohme, f_pohme, x_pohme = 1 
    REAL(KIND=rspec) :: delta_vsur, f_vsur, x_vsur = 1
    REAL(KIND=rspec) :: delta_Zeff, f_Zeff, x_Zeff = 1 
    REAL(KIND=rspec) :: delta_curt, f_curt, x_curt = 1 
    REAL(KIND=rspec) :: delta_P_eq, f_P_eq, x_P_eq = 1 
    REAL(KIND=rspec) :: delta_q_eq, f_q_eq, x_q_eq = 1 
    REAL(KIND=rspec) :: delta_elong, f_elong, x_elong = 1 
    REAL(KIND=rspec) :: delta_triang, f_triang, x_triang = 1 
    REAL(KIND=rspec) :: P_icrf_0, P_icrf_1, t_icrf_1
    REAL(KIND=rspec) :: P_nbi_0, P_nbi_1, t_nbi_1
        
    REAL(KIND=rspec) :: t, delta
    INTEGER :: i_t
       
!------------------------------------------------------------------------------------
!     
!   Namelists
!
!------------------------------------------------------------------------------------

!    namelist /static_state_data/
    
    namelist /evolving_model_data/ &
          Te_evolution_model, delta_Te, f_Te, &
          Ti_evolution_model, delta_Ti, f_Ti, &
          n_evolution_model, delta_n, f_n, &
          curr_bootstrap_evolution_model, delta_curr_bootstrap, f_curr_bootstrap, &
          curr_ohmic_evolution_model, delta_curr_ohmic, f_curr_ohmic, &
          pohme_evolution_model, delta_pohme, f_pohme, &
          vsur_evolution_model, delta_vsur, f_vsur, &
          Zeff_evolution_model, delta_Zeff, f_Zeff, &
          curt_evolution_model, delta_curt, f_curt, &
          P_eq_evolution_model, delta_P_eq, f_P_eq, &
          q_eq_evolution_model, delta_q_eq, f_q_eq, &
          elong_evolution_model, delta_elong, f_elong, &
          triang_evolution_model, delta_triang, f_triang, &
          power_ic_evolution_model, P_icrf_0, P_icrf_1, t_icrf_1, &
          power_nbi_evolution_model, P_nbi_0, P_nbi_1, t_nbi_1

    namelist /internal_state_data/ &
          n_step_calls, &
    	  x_Te,  x_Ti,  x_n, x_curr_bootstrap, x_curr_ohmic, x_pohme, x_vsur, x_Zeff, &
    	  x_curt,  x_P_eq, x_q_eq, x_elong, x_triang
!------------------------------------------------------------------------------------
!     
!   Get command line arguments
!
!------------------------------------------------------------------------------------

    iarg = command_argument_count()
    if(iarg .ne. 4) then
        print*, 'model_epa: '
        print*, ' command line args = cur_state_file next_state_file cur_eqdsk_file  mode time_stamp'
        print*, 'incorrect command line arguments'
        call exit(1)
    end if
      
      call get_command_argument(1,cur_state_file)
      call get_command_argument(2,cur_eqdsk_file)
      call get_command_argument(3,mode)
      call get_command_argument(4,time_stamp)
      
    WRITE (*,*)
    WRITE (*,*) 'model_epa_ps_file_init_no_RF'      
      print*, 'cur_state_file = ', trim(cur_state_file)
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
    WRITE (*,*) 'model_epa_ps_file_init_no_RF INIT'
          
    !---------------------------------------------------------------------------------
    !     
    !  Get current plasma state 
    !
    !---------------------------------------------------------------------------------
        
      call   ps_get_plasma_state(ierr, trim(cur_state_file))
      if(ierr .ne. 0) then 
          WRITE (*,*) 'model_epa_ps_file_init_no_RF.f90: call ps_get_plasma_state failed'
          call exit(1)
      end if
          
          
    !---------------------------------------------------------------------------------
    !     
    !  Get static EPA plasma state data and model data from input_namelist_files
    !
    !---------------------------------------------------------------------------------

       OPEN (unit=21, file=TRIM(input_namelist_file), status='old', &
            action='read', iostat=istat, form='formatted')
            IF (istat /= 0 ) THEN
                CALL SWIM_error ('open', 'model_epa_ps_file_init_no_RF.f90',TRIM(input_namelist_file))
                ierr = istat
                WRITE (*,*) 'cannot open input_namelist_file'
                call exit(1)
            END IF
        ierr = 0

!        read(21, nml=static_state_data)
!        WRITE (*, nml = static_state_data)

        read(21, nml=evolving_model_data)
        WRITE (*, nml = evolving_model_data)

        CLOSE (21)
      
          
    !---------------------------------------------------------------------------------
    !     
    !  Get aux plasma state to copy epa data from.  It comes in from input_state_file
    !
    !---------------------------------------------------------------------------------
        
      CALL ps_get_plasma_state(ierr, TRIM(input_state_file), aux)
      if(ierr .ne. 0) then
          WRITE (*,*) 'call to ps_get_plasma_state for aux state failed'
          CALL exit(1)
      end if

    !--------------------------------------------------------------------------
    !
    !   Copy data in PLASMA and EQ sections of 'aux' state to current (i.e. 'ps') state
    !
    !--------------------------------------------------------------------------
	  
	  CALL ps_cclist_remove("*", cclist, iwarn)
	  CALL ps_cclist_add("PLASMA" ,cclist, iwarn)
	  CALL ps_cclist_add("EQ", cclist, iwarn)
	  CALL ps_cclist_add("GAS", cclist, iwarn)
	  CALL ps_cclist_add("RAD", cclist, iwarn)
	  CALL PS_COPY_DIMS(aux, ps, 0 ,ierr, cclist, 0)
          CALL PS_COPY_PLASMA_STATE(aux, ps, ierr,  cclist, inodims = 1)
      if(ierr .ne. 0) then
         WRITE (*,*) 'call failed to PS_COPY_PLASMA_STATE for aux state to ps state'
         CALL exit(1)
      end if
      WRITE (*,*) 'Updated state from aux plasma state file ', TRIM(input_state_file)

	!-------------------------------------------------------------------------- 
	! Get ICRF machine configuration data from machine description
	! namelist file.  Put into 'output' plasma state.  Note: we are using another of the
	! default plasma state objects 'output'
	!--------------------------------------------------------------------------

	CALL ps_store_plasma_state(ierr, 'initial_output_ps.cdf', output)

	call ps_mdescr_namelist_read(.False., 'mdescr_namelist.dat', ' ',  ' ', output, ierr)
	IF (ierr .ne. 0) THEN
		print*, 'Could not get namelist mdescr'
		call exit(1)
	END IF
	print*,   '   number of icrf sources = ', ps%nicrf_src

	!--------------------------------------------------------------------------
	! Copy Zimp1 and Zimp2 from ps to output to avoid these being clobbered in the
	! plasma state copy from output to ps
	!--------------------------------------------------------------------------
	
	output%Zimp1 = ps%Zimp1
	output%Zimp2 = ps%Zimp2
	
	!--------------------------------------------------------------------------
	! Copy data in RF_IC sections of 'output' state to current (i.e. 'ps') state. 
	!--------------------------------------------------------------------------
	 
	 CALL ps_cclist_remove("*", cclist, iwarn)
	 CALL ps_cclist_add("IC" ,cclist, iwarn)
	 
	 CALL PS_COPY_DIMS(output, ps, 0 ,ierr, cclist, 0)
	 if (ierr .ne. 0) then
		 print*, 'call failed to PS_COPY_DIMS for output state to ps state'
		 call exit(1)
	 end if
	 CALL PS_COPY_PLASMA_STATE(output, ps, ierr, cclist, 1, 0)
	 if (ierr .ne. 0) then
		 print*, 'call failed to PS_COPY_PLASMA_STATE for output state to ps state'
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
    ! Load equilibrium data from input_eqdsk_file if one is provided.
    ! If there is no input_eqdsk_file the ps_update_equilibrium will fail, in which 
    ! case use the equilibrium data from the input_plasma_state.
    !--------------------------------------------------------------------------

	CALL ps_update_equilibrium(ierr, TRIM(input_eqdsk_file) )
	IF (ierr .ne. 0) THEN
		WRITE (*,*) 'model_epa INIT: ps_update_equilibrium call failed'
		WRITE (*,*) 'model_epa INIT: using equilibrium data from', input_state_file
    ELSE
        WRITE (*,*) 'Updated equilibrium from aux eqdsk file ', TRIM(input_eqdsk_file)
	END IF
	
!-------------------------------------------------------------------------- 
! Over-ride ps%eqdsk_file path in input state file with cur_eqdsk_file path
!--------------------------------------------------------------------------

	ps%eqdsk_file = cur_eqdsk_file
	
    !-------------------------------------------------------------------------- 
    ! Over-ride ps%t1 in input state file with time_stamp
    !--------------------------------------------------------------------------

	READ (time_stamp, '(f12.6)') ps%t1

    !-------------------------------------------------------------------------- 
    ! Over-ride initial power_IC if a time dependent model is provided
    !--------------------------------------------------------------------------
        
    IF (TRIM(power_ic_evolution_model) == 'none') THEN
         	
         CONTINUE	! ICRF power doesn't change in time
                  
    ELSE IF (TRIM(power_ic_evolution_model) == 'change_power') THEN
			
		IF (ps%t1 < t_icrf_1) THEN
			ps%power_ic =  P_icrf_0/ps%nicrf_src
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (INIT): ps%power_ic = ', ps%power_ic
		END IF

	ELSE
		WRITE(*,*) 'model_epa_ps_file_init_no_RF (INIT): unrecognized ICRF model ', &
		& power_ic_evolution_model			
		CALL exit(1)
	END IF

    !-------------------------------------------------------------------------- 
    ! Over-ride initial power_nbi if a time dependent model is provided
    !--------------------------------------------------------------------------
        
    IF (TRIM(power_nbi_evolution_model) == 'none') THEN
         	
         CONTINUE	! NB power doesn't change in time
                  
    ELSE IF (TRIM(power_nbi_evolution_model) == 'change_power') THEN
			
		IF (ps%t1 < t_nbi_1) THEN
			ps%power_nbi =  P_nbi_0/ps%nbeam
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (INIT): ps%power_nbi = ', ps%power_nbi
                        WRITE(*,*) 'model_epa_ps_file_init_no_RF (INIT): ps%nbeam =', ps%nbeam

		END IF

	ELSE
		WRITE(*,*) 'model_epa_ps_file_init_no_RF (INIT): unrecognized NBI model ', &
		& power_nbi_evolution_model			
		CALL exit(1)
	END IF

    !-------------------------------------------------------------------------- 
    ! Store initial plasma state
    !--------------------------------------------------------------------------

	CALL    ps_STORE_PLASMA_STATE(ierr, cur_state_file)
	IF (ierr .ne. 0) THEN
		WRITE (*,*) 'model_epa_ps_file_init_no_RF INIT: STORE_PLASMA_STATE call failed'
		CALL exit(1)
	ELSE 
		WRITE (*,*) "model_epa_ps_file_init_no_RF: Stored initial Plasma State" 
	END IF

    !---------------------------------------------------------------------------------
    !  Create scratch file for time dependent model data
    !---------------------------------------------------------------------------------

   OPEN (unit=21, file=TRIM(internal_state_file), status='unknown', &
		action='write', iostat=istat, form='formatted')
		IF (istat /= 0 ) THEN
			CALL SWIM_error ('open', 'model_epa_ps_file_init_no_RF.f90',TRIM(internal_state_file))
			ierr = istat
			WRITE (*,*) 'cannot open input_namelist_file'
			CALL exit(1)
		END IF
	ierr = 0

	write(21, nml=internal_state_data)
	CLOSE (21)
        
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
	if(ierr .ne. 0) then
	    WRITE (*,*)'call failed to ps_get_plasma_state to load profiles '
	    CALL exit(1)
	end if
          
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
			WRITE (*,*) 'cannot open input_namelist_file'
			CALL exit(1)
		END IF
	ierr = 0
	read(21, nml=evolving_model_data)
	CLOSE (21)
        
    !---------------------------------------------------------------------------------
    !  Get scratch data for time dependent model
    !---------------------------------------------------------------------------------

   OPEN (unit=21, file=TRIM(internal_state_file), status='old', &
		action='read', iostat=istat, form='formatted')
		IF (istat /= 0 ) THEN
			CALL SWIM_error ('open', 'model_epa.f90',TRIM(internal_state_file))
			ierr = istat			
			WRITE (*,*) 'cannot open input_namelist_file'
			CALL exit(1)
		END IF
	ierr = 0

	read(21, nml=internal_state_data)
	CLOSE (21)
        
    !--------------------------------------------------------------------------
    !
    !   Evolve state data
    !
    !--------------------------------------------------------------------------
        
        READ (time_stamp, '(f12.6)' ) t
        WRITE (*,*) 'model_epa STEP: time stamp = ', t, '   ps%t1 = ', ps%t1
        
        n_step_calls = n_step_calls + 1 ! increment counter for number of calls

        !---------------------------------------------------------------------------------
        ! Evolve data in the thermal profile state arrays
        !---------------------------------------------------------------------------------
         
        IF (TRIM(Te_evolution_model) == 'none') THEN
         	
         	CONTINUE	! Te profiles don't change in time
                  
        ELSE IF (TRIM(Te_evolution_model) == 'scale_profiles') THEN
         	ps%Ts(:,0) = ps%Ts(:,0)*(1.0_rspec + delta_Te)

        ELSE IF (TRIM(Te_evolution_model) == 'relax_profiles') THEN
        	CALL relax_profile(ps%Ts(:,0), delta_Te, f_Te, x_Te)

		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& Te_evolution_model			
			CALL exit(1)
		END IF
         
        IF (TRIM(Ti_evolution_model) == 'none') THEN
         	
         	CONTINUE	! Ti profiles don't change in time
                  
        ELSE IF (TRIM(Ti_evolution_model) == 'scale_profiles') THEN
         	ps%Ts(:,1:ps%nspec_th) = ps%Ts(:,1:ps%nspec_th)*(1.0_rspec + delta_Ti)
			ps%Ti(:) = ps%Ts(:,1)

        ELSE IF (TRIM(Ti_evolution_model) == 'relax_profiles') THEN
			CALL relax_species_profiles(ps%Ts(:,1:ps%nspec_th), delta_Ti, f_Ti, x_Ti)
			ps%Ti(:) = ps%Ts(:,1) ! Set Ti profile to first thermal ion species
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& Ti_evolution_model			
			CALL exit(1)
		END IF
         
        IF (TRIM(n_evolution_model) == 'none') THEN
         	
         	CONTINUE	! Density profiles don't change in time
                  
        ELSE IF (TRIM(n_evolution_model) == 'scale_profiles') THEN
			ps%ns(:, 0:ps%nspec_th) = ps%ns(:, 0:ps%nspec_th)*(1.0_rspec + delta_n)         

        ELSE IF (TRIM(n_evolution_model) == 'relax_profiles') THEN
			CALL relax_species_profiles(ps%ns, delta_n, f_n, x_n)

		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& n_evolution_model			
			CALL exit(1)
		END IF
         
        IF (TRIM(curr_bootstrap_evolution_model) == 'none') THEN
         	
         	CONTINUE	! curr_bootstrap profiles don't change in time
                  
        ELSE IF (TRIM(curr_bootstrap_evolution_model) == 'scale_profiles') THEN
         	ps%curr_bootstrap = ps%curr_bootstrap*(1.0_rspec + delta_curr_bootstrap)

        ELSE IF (TRIM(curr_bootstrap_evolution_model) == 'relax_profiles') THEN
        	CALL relax_profile(ps%curr_bootstrap, delta_curr_bootstrap, f_curr_bootstrap, &
        		x_curr_bootstrap)
         	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& curr_bootstrap_evolution_model			
			CALL exit(1)
		END IF

        IF (TRIM(curr_ohmic_evolution_model) == 'none') THEN
         	
         	CONTINUE	! curr_ohmic profiles don't change in time
                  
        ELSE IF (TRIM(curr_ohmic_evolution_model) == 'scale_profiles') THEN
         	ps%curr_ohmic = ps%curr_ohmic*(1.0_rspec + delta_curr_ohmic)

        ELSE IF (TRIM(curr_ohmic_evolution_model) == 'relax_profiles') THEN
        	CALL relax_profile(ps%curr_ohmic, delta_curr_ohmic, f_curr_ohmic, x_curr_ohmic)
         	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& curr_ohmic_evolution_model			
			CALL exit(1)
		END IF

       IF (TRIM(pohme_evolution_model) == 'none') THEN
         	
         	CONTINUE	! pohme profiles don't change in time
                  
        ELSE IF (TRIM(pohme_evolution_model) == 'scale_profiles') THEN
         	ps%pohme = ps%pohme*(1.0_rspec + delta_pohme)

        ELSE IF (TRIM(pohme_evolution_model) == 'relax_profiles') THEN
        	CALL relax_profile(ps%pohme, delta_pohme, f_pohme, x_pohme)
         	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& pohme_evolution_model			
			CALL exit(1)
		END IF
         
        IF (TRIM(vsur_evolution_model) == 'none') THEN
         	
         	CONTINUE	! vsur profiles don't change in time
                  
        ELSE IF (TRIM(vsur_evolution_model) == 'scale_profiles') THEN
         	ps%vsur = ps%vsur*(1.0_rspec + delta_vsur)

        ELSE IF (TRIM(vsur_evolution_model) == 'relax_profiles') THEN
        	CALL relax_scalar(ps%vsur, delta_vsur, f_vsur, x_vsur)
         	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& vsur_evolution_model			
			CALL exit(1)
		END IF
         
        IF (TRIM(Zeff_evolution_model) == 'none') THEN
         	
         	CONTINUE	! Zeff profiles don't change in time
                  
        ELSE IF (TRIM(Zeff_evolution_model) == 'scale_profiles') THEN
         	ps%Zeff = ps%Zeff*(1.0_rspec + delta_Zeff)

        ELSE IF (TRIM(Zeff_evolution_model) == 'relax_profiles') THEN
        	CALL relax_profile(ps%Zeff, delta_Zeff, f_Zeff, x_Zeff)
         	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& Zeff_evolution_model			
			CALL exit(1)
		END IF

        !---------------------------------------------------------------------------------
        ! Evolve data in the EQ component arrays
        !---------------------------------------------------------------------------------
         
        IF (TRIM(curt_evolution_model) == 'none') THEN
         	
         	CONTINUE	! curt profiles don't change in time
                  
        ELSE IF (TRIM(curt_evolution_model) == 'scale_profiles') THEN
         	ps%curt = ps%curt*(1.0_rspec + delta_curt)

        ELSE IF (TRIM(curt_evolution_model) == 'relax_profiles') THEN
        	CALL relax_profile(ps%curt, delta_curt, f_curt, x_curt)
         	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& curt_evolution_model			
			CALL exit(1)
		END IF
         
        IF (TRIM(P_eq_evolution_model) == 'none') THEN
         	
         	CONTINUE	! P_eq profiles don't change in time
                  
        ELSE IF (TRIM(P_eq_evolution_model) == 'scale_profiles') THEN
         	ps%P_eq = ps%P_eq*(1.0_rspec + delta_P_eq)

        ELSE IF (TRIM(P_eq_evolution_model) == 'relax_profiles') THEN
        	CALL relax_profile(ps%P_eq, delta_P_eq, f_P_eq, x_P_eq)
         	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& P_eq_evolution_model			
			CALL exit(1)
		END IF
         
        IF (TRIM(q_eq_evolution_model) == 'none') THEN
         	
         	CONTINUE	! q_eq profiles don't change in time
                  
        ELSE IF (TRIM(q_eq_evolution_model) == 'scale_profiles') THEN
         	ps%q_eq = ps%q_eq*(1.0_rspec + delta_q_eq)

        ELSE IF (TRIM(q_eq_evolution_model) == 'relax_profiles') THEN
        	CALL relax_profile(ps%q_eq, delta_q_eq, f_q_eq, x_q_eq)
         	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& q_eq_evolution_model			
			CALL exit(1)
		END IF
         
        IF (TRIM(elong_evolution_model) == 'none') THEN
         	
         	CONTINUE	! elong profiles don't change in time
                  
        ELSE IF (TRIM(elong_evolution_model) == 'scale_profiles') THEN
         	ps%elong = ps%elong*(1.0_rspec + delta_elong)

        ELSE IF (TRIM(elong_evolution_model) == 'relax_profiles') THEN
        	CALL relax_profile(ps%elong, delta_elong, f_elong, x_elong)
         	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& elong_evolution_model			
			CALL exit(1)
		END IF
         
        IF (TRIM(triang_evolution_model) == 'none') THEN
         	
         	CONTINUE	! triang profiles don't change in time
                  
        ELSE IF (TRIM(triang_evolution_model) == 'scale_profiles') THEN
         	ps%triang = ps%triang*(1.0_rspec + delta_triang)

        ELSE IF (TRIM(triang_evolution_model) == 'relax_profiles') THEN
        	CALL relax_profile(ps%triang, delta_triang, f_triang, x_triang)
         	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized profile model ', &
			& triang_evolution_model			
			CALL exit(1)
		END IF

    !--------------------------------------------------------------------------
    !    Change total ICRF power for next time step
    !--------------------------------------------------------------------------
        
        IF (TRIM(power_ic_evolution_model) == 'none') THEN
         	
         	CONTINUE	! ICRF power doesn't change in time
                  
        ELSE IF (TRIM(power_ic_evolution_model) == 'change_power') THEN
				
			IF (ps%t1 < t_icrf_1) THEN
				ps%power_ic =  P_icrf_0
			ELSE
				ps%power_ic =  P_icrf_1
				WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): ps%power_ic = ', ps%power_ic
			END IF
	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized ICRF model ', &
			& power_ic_evolution_model			
			CALL exit(1)
		END IF

    !--------------------------------------------------------------------------
    !    Change total NBI power for next time step
    !--------------------------------------------------------------------------
        
        IF (TRIM(power_nbi_evolution_model) == 'none') THEN
         	
         	CONTINUE	! NBI power doesn't change in time
                  
        ELSE IF (TRIM(power_nbi_evolution_model) == 'change_power') THEN
				
			IF (ps%t1 < t_nbi_1) THEN
				ps%power_nbi =  P_nbi_0/ps%nicrf_src
			ELSE
				ps%power_nbi =  P_nbi_1/ps%nbeam
				WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): ps%power_nbi = ', ps%power_nbi
                                WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): ps%nbeam =', ps%nbeam

			END IF
	
		ELSE
			WRITE(*,*) 'model_epa_ps_file_init_no_RF (STEP): unrecognized NBI model ', &
			& power_nbi_evolution_model			
			CALL exit(1)
		END IF

    !--------------------------------------------------------------------------    !
    ! Store the data in partial plasma_state file
    !--------------------------------------------------------------------------

		CALL PS_WRITE_UPDATE_FILE('EPA_'//cur_state_file, ierr)
        if(ierr .ne. 0) then
            WRITE(*,*) 'model_epa STEP: failed to write partial plasma state'
            CALL exit(1)
        end if
		WRITE (*,*) "Stored Partial EPA Plasma State"    
 
    !--------------------------------------------------------------------------
    ! Store the data in full plasma_state file
    !--------------------------------------------------------------------------
                
	CALL    ps_STORE_PLASMA_STATE(ierr, cur_state_file)
	WRITE (*,*) "model_epa STEP: Stored Plasma State"    
    
        
    !---------------------------------------------------------------------------------
    !  Write scratch data for time dependent model
    !---------------------------------------------------------------------------------

       OPEN (unit=21, file=TRIM(internal_state_file), status='old', &
            action='write', iostat=istat, form='formatted')
            IF (istat /= 0 ) THEN
                CALL SWIM_error ('open', 'model_epa.f90',TRIM(internal_state_file))
                ierr = istat
                WRITE (*,*) 'cannot open internal_state_file for output'
                CALL exit(1)
            END IF
        ierr = 0

        write(21, nml=internal_state_data)
        CLOSE (21)
        
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
    SUBROUTINE relax_profile(profile, delta, f, x)

      !  routine to model exponential relaxation of a profile to a multiple = f of
      !  its initial value.  At each call profile moves a fraction delta from it's
      !  present input value toward f*profile(at t = 0).  x is a running value that
      !  moves from 1 to f as things advance

      REAL(KIND=rspec), intent(in) :: delta  ! fraction to step
      REAL(KIND=rspec), intent(in) :: f  ! asymptotic multiple of profile
      REAL(KIND=rspec), intent(inout) :: x  ! keeps track of advancement
      REAL(KIND=rspec), dimension(:), intent(inout) :: profile ! real vector to relax
      REAL(KIND=rspec) :: delta0
      
      delta0 = delta*(f/x -1)
      x= (1.0_rspec + delta0)*x
      profile = (1.0_rspec + delta0)*profile

    END SUBROUTINE relax_profile

    !------------------------------------
    SUBROUTINE relax_species_profiles(profile, delta, f, x)

      !  routine to model exponential relaxation of a profile to a multiple = f of
      !  its initial value.  At each call profile moves a fraction delta from it's
      !  present input value toward f*profile(at t = 0).  x is a running value that
      !  moves from 1 to f as things advance

      REAL(KIND=rspec), intent(in) :: delta  ! fraction to step
      REAL(KIND=rspec), intent(in) :: f  ! asymptotic multiple of profile
      REAL(KIND=rspec), intent(inout) :: x  ! keeps track of advancement
      REAL(KIND=rspec), dimension(:,:), intent(inout) :: profile ! real vector to relax
      REAL(KIND=rspec) :: delta0
      
      delta0 = delta*(f/x -1)
      x= (1.0_rspec + delta0)*x
      profile = (1.0_rspec + delta0)*profile

    END SUBROUTINE relax_species_profiles

    !------------------------------------
    SUBROUTINE relax_scalar(s, delta, f, x)

      !  routine to model exponential relaxation of a profile to a multiple = f of
      !  its initial value.  At each call profile moves a fraction delta from it's
      !  present input value toward f*profile(at t = 0).  x is a running value that
      !  moves from 1 to f as things advance

      REAL(KIND=rspec), intent(in) :: delta  ! fraction to step
      REAL(KIND=rspec), intent(in) :: f  ! asymptotic multiple of profile
      REAL(KIND=rspec), intent(inout) :: x  ! keeps track of advancement
      REAL(KIND=rspec), intent(inout) :: s ! real vector to relax
      REAL(KIND=rspec) :: delta0
      
      delta0 = delta*(f/x -1)
      x= (1.0_rspec + delta0)*x
      s = (1.0_rspec + delta0)*s

    END SUBROUTINE relax_scalar

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
         CALL exit(1)
      ENDIF
    END SUBROUTINE ckerr


END PROGRAM model_epa_ps_file_init_no_RF

