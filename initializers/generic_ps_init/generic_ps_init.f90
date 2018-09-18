PROGRAM generic_ps_init

! version 1.0  (4/4/2016) Batchelor

!--------------------------------------------------------------------------
!
! Fortran code called by generic_ps_init.py.  The  Swiss army knife of Plasma State initializers.
! 
! This version combines several previous initializer routines and extends them.  There are
! 3 modes of initialization which must be specified by the config file variable INIT_MODE
! 
! INIT_MODE = minimal
! This is exactly the same as the previous minimal_state_init.py. It produces a CURRENT_STATE 
! that is empty except for some metadata:
! time variables - ps%t0, ps%t1, ps%tinit, and ps%tfinal 
! simulation identifiers - ps%tokamak_id, ps%shot_number, ps%run_id.  
! ps%Global_label is set to run_id_tokamak_id_shot_number.
! This data is set for all initialization modes, but for 'minimal' this is all the data
! included.
! 
! INIT_MODE = existing_ps_file
! This copies an existing input plasma state file and optionally an existing eqdsk file to
! CURRENT_STATE and CURRENT_EQDSK.  If the config parameter GENERATE_EQDSK is set to 'True'
! the CURRENT_EQDSK file is generated from equilibrium data in the INPUT_STATE_FILE.
! The INPUT_STATE_FILE and INPUT_EQDSK_FILE must be specified in the config file.
! 
! INIT_MODE = mdescr
! This initializes all machine description data from a plasma state machine description 
! file, e.g. <tokamak>.mdescr, as specified by config parameter MDESCR_FILE. In addition
! if a shot configuration file config parameter, SCONFIG_FILE, is specified, the shot config
! data is also loaded into CURRENT_STATE.  Machine description and shot configuration files
! are namelist files that can be read and loaded using Plasma State subroutines ps_mdescr_read()
! and ps_scongif_read().  Note:  machine description and shot configuration do not define
! the MHD equilibrium, so the equilibrium must be specified during further component 
! initializations
! 
! INIT_MODE = mixed (yet to be implemented)
! 
! Except for possibly mode = existing_ps_file, all modes call on the fortran helper code 
! generic_ps_file_init.f90 to interact with the Plasma State. The fortran code is also used
! in existing_ps_file mode to extract the CURRENT_EQDSK when GENERATE_EQDSK = true.
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
    INTEGER :: iarg, status = 0
    
    CHARACTER (len=256) :: cur_state_file, cur_eqdsk_file
    CHARACTER (len=256) ::  input_eqdsk_file = ' '
    CHARACTER (len=32) ::   mdescr_file = ' '
    CHARACTER (len=32) ::   sconfig_file = ' '
    CHARACTER(len=32) :: init_mode  
    CHARACTER(len=32) :: generate_eqdsk = 'False'

    LOGICAL :: file_exists


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
    IF (ierr .ne. 0) THEN
		CALL SWIM_error ('open', 'generic_ps_init.f90',ps_init_nml_file)
		WRITE (*,*) 'generic_ps_init.f90: Cannot open ', TRIM(ps_init_nml_file)
		call exit(1)
	END IF

	read(21, nml=ps_init_nml)
	CLOSE (21)
	WRITE (*, nml = ps_init_nml)

!------------------------------------------------------------------------------------
!     
!   Do initalizations from machine description file
!
!------------------------------------------------------------------------------------

    IF (TRIM(init_mode) == 'mdescr') THEN
        inquire(file=trim(mdescr_file), exist=file_exists)
        if(.not.file_exists)then
            write(*,*)'generic_ps_init : ERROR - mdescr_file not found'  
            write(*,*) trim(mdescr_file)
            call exit(1)
        endif
        write(*,*) 'generic_ps_init: mdescr_file = ', trim(mdescr_file)
        call ps_mdescr_read(trim(mdescr_file), ierr, state=ps)
    END IF

!------------------------------------------------------------------------------------
!     
!   Load shot configuration data from sconfig file
!
!------------------------------------------------------------------------------------

    IF (TRIM(sconfig_file) /= ' ') THEN
        inquire(file=trim(sconfig_file), exist=file_exists)
        if(.not.file_exists)then
            write(*,*)'generic_ps_init : ERROR - sconfig_file not found'  
            write(*,*) trim(sconfig_file)
            status = 1
            call exit(status)
        endif
        write(*,*) 'generic_ps_init: sconfig_file = ', trim(sconfig_file)
        call ps_sconfig_read(trim(sconfig_file), ierr, state=ps)
    END IF

!------------------------------------------------------------------------------------
!     
!   If init_mode = 'mixed' load current state file as aux and copy state data to ps
!
!------------------------------------------------------------------------------------

    IF (TRIM(init_mode) == 'mixed') THEN
          CALL ps_get_plasma_state(ierr, TRIM(cur_state_file), aux)
          if (ierr .ne. 0) then
              print*, 'call failed to ps_get_plasma_state for aux state'
              call exit(1)
          end if
          
         CALL PS_COPY_PLASMA_STATE(aux, ps, ierr, inodims = 1)
         if (ierr .ne. 0) then
             print*, 'call failed to PS_COPY_PLASMA_STATE for aux state to ps state'
             call exit(1)
         end if
    	


!--------------------------------------------------------------------------
! 
! Load equilibrium data from input_eqdsk_file if one is provided.
!
!--------------------------------------------------------------------------

    IF (TRIM(input_eqdsk_file) /= ' ') THEN
        inquire(file=trim(input_eqdsk_file), exist=file_exists)
        if(.not.file_exists)then
            write(*,*)'generic_ps_init : ERROR - input_eqdsk_file not found'  
            write(*,*) trim(input_eqdsk_file)
            status = 1
            call exit(status)
        endif
        write(*,*) 'generic_ps_init: input_eqdsk_file = ', trim(input_eqdsk_file)
		CALL ps_update_equilibrium(ierr, TRIM(input_eqdsk_file) )
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
           print*, 'generic_ps_init: failed to get_plasma_state'
           call exit(1)
        end if

        CALL ps_wr_geqdsk(ierr, cur_eqdsk_file)
        IF (ierr .ne. 0) THEN
            print*, 'Could not get generate eqdsk file from plasma state'
            call exit(1)
        END IF
    END IF
        

!------------------------------------------------------------------------------------
!     
!   Store initial plasma state.  If init_mode == 'minimal' the plasma state file is
!   initialized but otherwise completely empty.
!
!------------------------------------------------------------------------------------

    CALL PS_STORE_PLASMA_STATE(ierr, trim(cur_state_file))    
    WRITE (*,*) "generic_ps_init.f90: Stored initial Plasma State"          

END PROGRAM generic_ps_init

