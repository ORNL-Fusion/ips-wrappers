PROGRAM load_mdescr_sconfig

!--------------------------------------------------------------------------
! Version 0.0 11/2/2012 (Batchelor) 
!
! Helper program for init phase of IPS ICRF components.  It reads a plasma state file
! ps which has EPA data already initialized.  Then it reads namelist files containing
! ICRF machine description and shot configuration namelists using the plasma
! state functions ps_mdescr_read() and ps_sconfig_read() and stores this data in 
! temporary plasma state file , aux.  It copies the new ICRF data to the current plasma
! state and saves it.  The saved state file will require further initialization to
! allocate and fill the profile arrays as is normally done in the component init
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
!	2) Shot configuration namelist file: sconfig_namelist.dat.  sconfig files can be
!	written by plasma state function ps_sconfig_write(). Sample sconfig files are produced
!	by the plasma_state_test routine.  The ICRF relevant items in namelist sconfig are:
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
!	3)	The code takes one command-line argument, the name of the current plasma state
!		to be initialized.  The names of the namelist files are assumed to be
!		mdescr_namelist.dat and sconfig_namelist.dat.
!
! One other point: The EPA component is responsible for setting the number of ICRF
! sources, nicrf_src and allocating the arrays it dimensions.  This is because EPA
! controls the power on the different sources.  Note that nicrf_src also 
! appears in the mdescr namelist.  Although it appears implicitly as the number
! of elements in the namelist that are dimensioned nicrf_src.  For example icrf_src_name,
! ant_model and dx_fshield. The user must make sure this value is consistent
! with what is set by EPA.  I know, data should only be set in one place, but it
! would be complicated to fix, since the EPA may not be initialized from the same,
! or indeed any, mdescr namelist.  Just be careful at run time, ok.
!
! For that matter many of the arrays in the sconfig namelist are also dimensioned
! nicrf_src, for example freq_ic.  So the dimensons in the two namelist files must be
! consistent.
!
!--------------------------------------------------------------------------


    USE plasma_state_mod
    
    USE swim_global_data_mod, only : &
            & rspec, ispec, &               ! int: kind specification for real and integer
            & SWIM_name, SWIM_filename, &   ! derived data types: containing one character
                                                                        ! string
            & SWIM_error                    ! subroutine: a simple error handling routine
    
    IMPLICIT none
    
    INTEGER :: iarg, ierr = 0, iwarn = 0
    
    INTEGER :: cclist(ps_ccount)   ! component activation list 

    CHARACTER (len=256) :: cur_state_file
    
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
				
    
	!------------------------------------------------------------------------------------
	!   Get command line argument
	!------------------------------------------------------------------------------------

      call get_arg_count(iarg)
      if (iarg .ne. 1) then
         print*, 'model_RF_IC usage: '
         print*, ' command line arg = cur_state_file'
         call exit(1)
      end if
      
      call getarg(1,cur_state_file)
      
     WRITE (*,*)
     WRITE (*,*) 'load_mdescr_sconfig'      
     print*, 'cur_state_file = ', trim(cur_state_file)
      
	!------------------------------------------------------------------------------------
	!  Get current plasma state 
	!------------------------------------------------------------------------------------
			
    call ps_get_plasma_state(ierr, trim(cur_state_file))
    if (ierr .ne. 0) then
       print*, 'model_RF_IC:failed to get_plasma_state'
       call exit(1)
    end if


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
		 print*, 'model_RF_IC_2.f90: I give up'
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
    ! Store updated plasma state
    !--------------------------------------------------------------------------

        CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)
        IF (ierr .ne. 0) THEN
			print*,  'load_mdescr_sconfig: PS_STORE_PLASMA_STATE failed'
		call exit(1)
		ELSE
        	WRITE (*,*) "load_mdescr_sconfig: Stored updated Plasma State"    
        END IF


END PROGRAM load_mdescr_sconfig
