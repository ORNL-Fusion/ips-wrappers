PROGRAM set_ps_beam_power

! version 1.0 5/23/2013 (Batchelor) 

!--------------------------------------------------------------------------
!
!   Simple code to set NB power for the beam sources in plasma state to values specified
!   in command line arguments.  To be launched by nubeam component.  The purpose is to
!   allow beam powers to be changed individually from within the component.
!
!       The code requires 2 or more command-line arguments
!       1) path to the current plasma state file
!       2) string containing power level of beam_1 -> str_power_beam(1)
!          and optionallly
!       3) string containing power level of beam_2 ... -> str_power_beam(2), ...
!   
!
!       Don Batchelor
!       ORNL
!       Oak Ridge, TN 37831
!       batchelordb@ornl.gov
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
INTEGER :: i

CHARACTER (len=256) :: cur_state_file

!--------------------------------------------------------------------------
!   Internal data
!--------------------------------------------------------------------------

CHARACTER(len=12), ALLOCATABLE :: beam_name(:)
CHARACTER(len=32), ALLOCATABLE :: str_power_beam(:)
INTEGER :: nbeams
REAL(KIND=rspec), ALLOCATABLE :: power_beam(:)

WRITE (*,*)
WRITE (*,*) 'set_ps_beam_power'      

!------------------------------------------------------------------------------------
!   Get command line arguments
!------------------------------------------------------------------------------------

call get_arg_count(iarg)
if(iarg .lt. 2) then
	PRINT*, 'set_ps_beam_power usage: '
	PRINT*, ' command line args = cur_state_file and power level for one or more beams'
	stop 'incorrect command line arguments'
end if

call getarg(1,cur_state_file)
PRINT*, 'cur_state_file = ', trim(cur_state_file)
nbeams = iarg-1

!------------------------------------------------------------------------------------
!  Get current plasma state 
!------------------------------------------------------------------------------------
		
call ps_get_plasma_state(ierr, trim(cur_state_file))
if(ierr .ne. 0) then
   PRINT*, 'set_ps_beam_power:failed to get_plasma_state'
   stop 1
end if
	

!---------------------------------------------------------------------------------
!  Check that number of beams in command line matches number in plasma state
!---------------------------------------------------------------------------------
			 
IF ( allocated(ps%nbi_src_name) ) THEN
	 PRINT*, 'NB source arrays allocated in input Plasma State'             
	 IF ( ps%nbeam .ne. nbeams) THEN
	 	PRINT*, 'number of nbi sources in command line arg inconsistent with plasma state'
	 	STOP
	 END IF		
	 PRINT*,   '   number of nbi sources = ', ps%nbeam

	 allocate(beam_name(nbeams),stat=istat)
	 DO i = 1, SIZE(ps%nbi_src_name)
	 	beam_name(i) = trim(ps%nbi_src_name(i))
	 END DO
	 WRITE(*, fmt = '(a, 4(3x,a))')   '    names of sources = ', beam_name
	 PRINT*,   '   initial nbi source powers = ', ps%power_nbi
	 
ELSE		  
	 PRINT*, 'NB sources NOT allocated in initial Plasma State'
	 STOP				
END IF

!------------------------------------------------------------------------------------
!  Allocate and get new beam powers from command line args 
!------------------------------------------------------------------------------------

allocate(str_power_beam(nbeams),stat=istat)
allocate(power_beam(nbeams),stat=istat)

DO i = 1, nbeams
	call getarg(i+1,str_power_beam(i))
	read(str_power_beam(i), *) power_beam(i)  ! Convert powers from string to real
END DO

      

!-------------------------------------------------------------------------- 
! modify plasma state
!--------------------------------------------------------------------------

DO i = 1, nbeams
	ps%power_nbi(i) = power_beam(i)
END DO
	 
PRINT*,   '   new nbi source powers = ', ps%power_nbi

!--------------------------------------------------------------------------    !
! Store the data in plasma_state file
!--------------------------------------------------------------------------

CALL PS_STORE_PLASMA_STATE(ierr, cur_state_file)

WRITE (*,*) "set_ps_beam_power: Stored Plasma State"    

END PROGRAM set_ps_beam_power
