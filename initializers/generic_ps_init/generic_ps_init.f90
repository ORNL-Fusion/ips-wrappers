PROGRAM generic_ps_init

!--------------------------------------------------------------------------
!
! Fortran code called by generic_ps_init.py.  The  Swiss army knife of Plasma State initializers.
! 
! This version combines several previous initializer routines and extends them.  There are
! X modes of initialization to be specified by the config file variable INIT_MODE
! 
! INIT_MODE = minimal
! This fortran code only generates a completely empty current state file.
! The metadata that was previously inserted into CURRENT_STATE by minimal_state_init.f90  i.e. 
! time variables - ps%t0, ps%t1, ps%tinit, and ps%tfinal 
! simulation identifiers - ps%tokamak_id, ps%shot_number, ps%run_id.  
! ps%Global_label -> run_id_tokamak_id_shot_number.
! is now done in the python, generic_ps_init.py.  
! 
! INIT_MODE = existing_ps_file
! The copy of the input plasma state to CURRENT_PLASMA_STATE is done in the python.  However if the config variable
! 'GENERATE_EQDSK_FROM_STATE' is 'True', the Plasma State subroutine ps_wr_eqdsk() is used to generate it.
! 
! 
! INIT_MODE = mdescr
! Machine description data is read from a plasma state machine description file <tokamak>.mdescr 
! using the plasma state function ps_mdescr_read().  

! INIT_MODE = mixed (yet to be implemented)
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
    
    INTEGER :: istat, ierr = 0
    INTEGER :: iarg
    
	CHARACTER (len=256) :: cur_state_file, cur_eqdsk_file
	CHARACTER (len=256) :: 	mdescr_file
	CHARACTER (len=256) :: 	sconfig_file = ' '
	CHARACTER (len=256) :: 	input_eqdsk_file = ' '
	CHARACTER(len=32) :: init_mode	
	CHARACTER(len=32) :: generate_eqdsk = 'False'



!--------------------------------------------------------------------------
!
!   Internal data: None
!
!--------------------------------------------------------------------------

	REAL(KIND=rspec) :: t, tinit, tfinal
	CHARACTER(LEN=*), PARAMETER :: ps_init_nml_file = 'generic_ps_init.nml'

!------------------------------------------------------------------------------------
!     
!   Namelists
!
!------------------------------------------------------------------------------------

    namelist /ps_init_nml/ &
          init_mode, generate_eqdsk, cur_state_file, cur_eqdsk_file, &
          mdescr_file, input_eqdsk_file, sconfig_file
          
           
	WRITE (*,*)
	WRITE (*,*) 'generic_ps_init'      

!---------------------------------------------------------------------------------
!     
!  Get init configuration data from ps_init_nml_file
!
!---------------------------------------------------------------------------------

    OPEN (unit=21, file = 'generic_ps_init.nml', status = 'old',   &
         form = 'formatted', iostat = ierr)
    IF (ierr .ne. 0) STOP 'cannot open EPA_model.nml'


!        OPEN (unit=21, file=TRIM(ps_init_nml_file), status='unknown', &
!             action='read', iostat=istat, form='formatted')
!             IF (istat /= 0 ) THEN
!                 CALL SWIM_error ('open', 'generic_ps_init.f90',ps_init_nml_file)
!                 ierr = istat
!                 WRITE (*,*) 'generic_ps_init.f90: Cannot open ', TRIM(ps_init_nml_file)
!                 stop 'generic_ps_init.f90: Cannot open ps_init_nml_file'
!             END IF
!         ierr = 0

        read(21, nml=ps_init_nml)
        CLOSE (21)
        WRITE (*, nml = ps_init_nml)

!------------------------------------------------------------------------------------
!     
!   Do initalizations from machine description file
!
!------------------------------------------------------------------------------------

	IF (TRIM(init_mode) == 'mdescr') THEN
		call ps_mdescr_namelist_read(.False., trim(mdescr_file), ' ',  &
				TRIM(input_eqdsk_file), ps, ierr)
		IF (ierr .ne. 0) THEN
			print*, 'Could not get namelist mdescr'
			call exit(1)
		END IF
	END IF

!------------------------------------------------------------------------------------
!     
!   Load shot configuration data from sconfig file
!
!------------------------------------------------------------------------------------

	IF (TRIM(sconfig_file) /= ' ') THEN
		call ps_sconfig_namelist_read(.False., TRIM(sconfig_file), ' ',  ' ', ps, ierr)
		IF (ierr .ne. 0) THEN
			print*, 'Could not get namelist sconfig'
			call exit(1)
		END IF
	END IF
	
!------------------------------------------------------------------------------------
!     
!   Extract eqdsk file from plasma state
!
!------------------------------------------------------------------------------------
	
	IF (TRIM(generate_eqdsk) == 'True') THEN
	!  Get current plasma state 
			
		call ps_get_plasma_state(ierr, trim(cur_state_file))
		if(ierr .ne. 0) then
		   print*, 'model_EPA:failed to get_plasma_state'
		   stop 1
		end if

		CALL ps_wr_geqdsk(ierr, cur_eqdsk_file)
		IF (ierr .ne. 0) THEN
			print*, 'Could not get generate eqdsk file from plasma state'
			call exit(1)
		END IF
	END IF
		

!------------------------------------------------------------------------------------
!     
!   Store initial plasma state.  If init_mode == 'minimal' the plasma statei file is
!   initialized but otherwise completely empty.
!
!------------------------------------------------------------------------------------

	CALL PS_STORE_PLASMA_STATE(ierr, trim(cur_state_file))
	
	WRITE (*,*) "generic_ps_init.f90: Stored initial Plasma State"    		

END PROGRAM generic_ps_init

