PROGRAM minimal_state_init

! version 2.0 2/2/2010 (Batchelor)
! version 1.1 2/17/2009 (Batchelor

!--------------------------------------------------------------------------
!
! Simulation init component that produces a plasma state that is empty except for:
! time variables - t0, t1, tinit, and tfinal, 
! and
! simulation identifiers - tokamak_id, shot_number, run_id.  
! Also the plasma state variable Global_label is set to run_id_tokamak_id_shot_number.
!
!   Note: There is a potential problem with Global_label.  
!   Global_lable is dimensioned CHARACTER(*80) in the state, whereas the other identifiers
!   are CHARACTER(*32).  So if they were all maximal length, Global_label would turn out 
!   to be 3*32 + 2 = 98.  (The extra 2 are the joining underscores). I hope the value 
!   would just be truncated.
!
! This code is launched by the component script minimal_state_init.py
! 
! This version takes 7 command line arguments: the path to the current plasma state, 
! tokamk_id, shot_number, run_id, tinit, tfinal of the time loop, and time_stamp.
!
! N.B. Both ps%t0 and ps%t1 are set to the value time_stamp.  ps%tinit and ps%tfinal
!      are generated by the Python component from the TIME_LOOP variable in the
!      simulation config file.  Note that the initial ps%t0 can be different from 
!      ps%tinit, as might be needed in a restart.

!! N.B. The other plasma state files that were produced by the previous versions of
!      this fortran code are now produced in the minimal_state_init.py. These files
!      include prior_state file and next file as well as the dummy files: cur_cql_file
!      cur_eqdsk_file, cur_dql_file, and cur_jsdsk_file.
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

!------------------------------------------------------------------------------------
!     
!   Get command line arguments
!
!------------------------------------------------------------------------------------

      call get_arg_count(iarg)
    if(iarg .ne. 7) then
        print*, 'minimal_state_init: '
	 	print*, ' command line args = cur_state_file, tokamak_id, &
				& shot_number, run_id, tinit, tfinal, time_stamp'
	 	stop 'incorrect command line arguments'
    end if
      
      call getarg(1,cur_state_file)
      call getarg(2,tokamak_id)
      call getarg(3,shot_number)
      call getarg(4,run_id)
      call getarg(5,char_tinit)
      call getarg(6,char_tfinal)
      call getarg(7,time_stamp)

      
	WRITE (*,*)
	WRITE (*,*) 'minimal_state_init'      
      print*, 'cur_state_file = ', trim(cur_state_file)      
      print*, 'tokamak = ', trim(tokamak_id)
      print*, 'shot number = ', trim(shot_number)
      print*, 'run_id = ', trim(run_id)
      print*, 'tinit = ', trim(char_tinit)
      print*, 'tfinal = ', trim(char_tfinal)
      print*, 'time_stamp = ', trim(time_stamp)



!------------------------------------------------------------------------------------
!     
!   Put run idnetifiers into plasma state obtained from USE plasma_state_mod
!
!------------------------------------------------------------------------------------

	ps%tokamak_id = trim( tokamak_id(1:32) )
	ps%RunID = trim( run_id(1:32) )
	
	read (shot_number, *) ps%shot_number
	
	ps%Global_label = trim(run_id)//'_'//trim(tokamak_id)//'_'//trim(shot_number)
 
!-----------------------------------
! Insert time at beginning and end of time step (both are time_stamp on initialization)
!-----------------------------------
	READ (time_stamp, '(f12.6)') t
	WRITE (*,*) 'minimal_state_init INIT: time stamp t = ', t
	 
	ps%t0 = t
	ps%t1 = ps%t0

	READ (char_tinit, '(f12.6)') tinit
	READ (char_tfinal, '(f12.6)') tfinal
	WRITE (*,*) 'minimal_state_init INIT: tinit = ', tinit, '   tfinal = ', tfinal
	 
	ps%tinit = tinit
	ps%tfinal = tfinal

!------------------------------------------------------------------------------------
!     
!   Store minimal initial plasma state
!
!------------------------------------------------------------------------------------

	CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)
	
	WRITE (*,*) "minimal_state_init.f90: Stored initial Plasma State"    

END PROGRAM minimal_state_init

