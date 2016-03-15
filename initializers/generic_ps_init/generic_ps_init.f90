PROGRAM generic_ps_init

!--------------------------------------------------------------------------
!
! Fortran code called by generic_ps_init.py.  The  Swiss army knife of Plasma State initializers.
! 
! This version combines several previous initializer routines and extends them.  There are
! X modes of initialization to be specified by the config file variable INIT_MODE
! 
! N.B. The metadata that was previously inserted into CURRENT_STATE by minimal_state_init.f90  i.e. 
! time variables - ps%t0, ps%t1, ps%tinit, and ps%tfinal 
! simulation identifiers - ps%tokamak_id, ps%shot_number, ps%run_id.  
! ps%Global_label -> run_id_tokamak_id_shot_number.
! is now done in the python, generic_ps_init.py.  If INIT_MODE = minimal, this fortran code is not run at all.
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
    
	CHARACTER (len=256) :: cur_state_file
	CHARACTER*32 :: tokamak_id, shot_number, run_id
    CHARACTER(len=32) :: time_stamp, char_tinit, char_tfinal


!--------------------------------------------------------------------------
!
!   Internal data: None
!
!--------------------------------------------------------------------------

	REAL(KIND=rspec) :: t, tinit, tfinal
	PARAMETER, CHARACTER* :: ps_init_nml_file = 'generic_ps_init.nml'

!------------------------------------------------------------------------------------
!     
!   Namelists
!
!------------------------------------------------------------------------------------

    namelist /genric_ps_init/ &
          init_mode, cur_state_file, mdescr_file
           
	WRITE (*,*)
	WRITE (*,*) 'generic_ps_init'      

!---------------------------------------------------------------------------------
!     
!  Get init configuration data from ps_init_nml_file
!
!---------------------------------------------------------------------------------

       OPEN (unit=21, file=TRIM(input_namelist_file), status='unknown', &
            action='write', iostat=istat, form='formatted')
            IF (istat /= 0 ) THEN
                CALL SWIM_error ('open', 'generic_ps_init.f90',ps_init_nml_file)
                ierr = istat
                stop 'cannot open ps_init_nml_file'
            END IF
        ierr = 0

        read(21, nml=genric_ps_init)
        CLOSE (21)
        WRITE (*, nml = genric_ps_init)

!------------------------------------------------------------------------------------
!     
!   Do initalizations from input plasma state file
!
!------------------------------------------------------------------------------------

	IF TRIM(init_mode) == 'existing_ps_file'
    	call ps_mdescr_read(trim(mdescr_file), ierr, state=ps)

!------------------------------------------------------------------------------------
!     
!   Do initalizations from machine description file
!
!------------------------------------------------------------------------------------

	IF TRIM(init_mode) == 'mdescr'
    	call ps_mdescr_read(trim(mdescr_file), ierr, state=ps)

!------------------------------------------------------------------------------------
!     
!   Store minimal initial plasma state
!
!------------------------------------------------------------------------------------

	CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)
	
	WRITE (*,*) "generic_ps_init.f90: Stored initial Plasma State"    

	ELSE
		WRITE (*,*) 'Unknown initialization moe = ', init_mode

END PROGRAM generic_ps_init

