PROGRAM write_medescr_sconfig

! Version 1.0 4/8/2016 (Batchelor)
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

    INTEGER :: ierr = 0, iarg = 0
 
    CHARACTER (len=256) :: state_file        

       
!------------------------------------------------------------------------------------
!     
!   Namelists
!
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
!     
!   Get command line arguments
!
!------------------------------------------------------------------------------------

    iarg = command_argument_count()
    if(iarg .ne. 1) then
        print*, 'write_medescr_sconfig: '
        print*, ' command line args = state_file'
        print*, 'incorrect command line arguments'
        call exit(1)
    end if
      
      call get_command_argument(1,state_file)
      
    WRITE (*,*)
    WRITE (*,*) 'write_medescr_sconfig'      
      print*, 'state_file = ', trim(state_file)
                
!---------------------------------------------------------------------------------
!     
!  Get plasma state 
!
!---------------------------------------------------------------------------------
        
      call   ps_get_plasma_state(ierr, trim(state_file))
      if(ierr .ne. 0) then 
          WRITE (*,*) 'write_medescr_sconfig.f90: call ps_get_plasma_state failed'
          call exit(1)
      end if
                
!---------------------------------------------------------------------------------
!     
!  Write mdescr and sconfig namelists 
!
!---------------------------------------------------------------------------------
          
	  write(iout,*) '   mdescr namelist write test...'
	  call ps_mdescr_write('mdescr_sample_namelist.dat',ierr)
	  call ckerr('ps_mdescr_write')

	  write(iout,*) '   sconfig namelist read test...'
	  call ps_sconfig_write('sconfig_sample_namelist.dat',ierr)
	  call ckerr('ps_sconfig_write')

!---------------------------------------------------------------------------------
!     
! read these namelists into a separate state object
!
!---------------------------------------------------------------------------------

	  write(iout,*) '   mdescr namelist read test...'
	  call ps_mdescr_read('mdescr_namelist.dat',ierr, state=aux)
	  call ckerr('ps_mdescr_write')

	  write(iout,*) '   sconfig namelist read test...'
	  call ps_sconfig_read('sconfig_namelist.dat',ierr, state=aux)
	  call ckerr('ps_sconfig_write')

	  ! store to file
	  CALL ps_store_plasma_state(ierr, state=aux)

	  write(iout,*) ' namelist I/O tests completed.'

    !------------------------------------
    SUBROUTINE ckerr(sbrtn)
      character*(*), intent(in) :: sbrtn

      IF(ierr.NE.0) then
         write(6,*) ' ?write_medescr_sconfig: error in call: '//trim(sbrtn)
         CALL exit(1)
      ENDIF
    END SUBROUTINE ckerr

END PROGRAM write_medescr_sconfig

