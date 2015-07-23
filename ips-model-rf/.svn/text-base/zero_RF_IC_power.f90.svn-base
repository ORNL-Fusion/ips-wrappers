PROGRAM zero_RF_IC_power

! version 0.0 5/3/2010 (Batchelor)

!--------------------------------------------------------------------------
!
!   Simple code to set all ICRF source profiles in plasma state to zero and
!   then write a partial plasma state.  This called from RF_IC python component
!   to produce a valid update partial plasma state without running an ICRF code.
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
!   Takes one command line argument -> current plasma state file
!	Writes a partial state -> RF_IC_ + current_state_file name
!
!--------------------------------------------------------------------------


    USE plasma_state_mod
    
    USE swim_global_data_mod, only : &
            & rspec, ispec, &               ! int: kind specification for real and integer
            & SWIM_name, SWIM_filename, &   ! derived data types: containing one character
                                                                        ! string
            & SWIM_error                    ! subroutine: a simple error handling routine
    
    IMPLICIT none
    
    INTEGER :: ierr = 0, iwarn = 0, iarg
    
    INTEGER :: cclist(ps_ccount)   ! component activation list 

    CHARACTER (len=256) :: cur_state_file

	!------------------------------------------------------------------------------------
	!   Get command line arguments
	!------------------------------------------------------------------------------------

      call get_arg_count(iarg)
      if(iarg .ne. 1) then
         print*, 'zero_RF_IC_power usage: '
         print*, ' command line args = cur_state_file'
         stop 'incorrect command line arguments'
      end if
      
      call getarg(1,cur_state_file)
            
	!------------------------------------------------------------------------------------
	!  Get current plasma state 
	!------------------------------------------------------------------------------------
			
    call ps_get_plasma_state(ierr, trim(cur_state_file))
    if(ierr .ne. 0) then
       print*, 'zero_RF_IC_power:failed to get_plasma_state'
       stop 1
    end if

		
    !-------------------------------------------------------------------------- 
    !  Set all outputs to zero
    ! --------------------------------------------------------------------------

		ps%picrf_srcs = 0.        
		ps%picrf_totals = 0.        
		ps%picth = 0.
		ps%curich = 0.
	
    !--------------------------------------------------------------------------    !
    ! Store the data in plasma_state file
    !--------------------------------------------------------------------------

	CALL PS_WRITE_UPDATE_FILE('RF_IC_'//cur_state_file, ierr)

	WRITE (*,*) "zero_RF_IC_power: Stored Partial RF Plasma State with zero outputs"    

END PROGRAM zero_RF_IC_power

