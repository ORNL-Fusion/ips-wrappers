Module toric_utils_mod

private
public:: getlun, assert

contains

SUBROUTINE getlun (ilun,ierr)
  !
  !-----------------------------------------------------------------------
  !
  ! ****** Return an unused logical unit identifier.
  !
  !-----------------------------------------------------------------------
  !
  ! ****** Upon successful completion (IERR=0), the first
  ! ****** unused logical unit number between MINLUN and
  ! ****** MAXLUN, inclusive, is returned in variable ILUN.
  ! ****** If all units between these limits are busy,
  ! ****** IERR=1 is returned.
  !
  !-----------------------------------------------------------------------
  !
  INTEGER, INTENT(OUT) :: ilun, ierr
  !
  !-----------------------------------------------------------------------
  !
  ! ****** Range of valid units.
  !
  INTEGER, PARAMETER :: minlun=30, maxlun=99
  LOGICAL :: busy
  !
  !-----------------------------------------------------------------------
  !
  ierr=0
  !
  ! ****** Find an unused unit number.
  !
  DO i=minlun,maxlun
     INQUIRE (unit=i,opened=busy)
     IF (.NOT.busy) THEN
        ilun=1
        RETURN
     END IF
  END DO
  !
  ! ****** Fall through here if all units are busy.
  !
  ierr=1
  RETURN

END subroutine getlun

subroutine assert( lcond, mesg, ivalue )
  logical, intent(in) :: lcond
  character*(*),intent(in) :: mesg
  integer, intent(in) :: ivalue

  if (.not.lcond) then
     write(0,*) mesg,ivalue
     stop '** assertion error **'
  endif
  return
end subroutine assert

end Module toric_utils_mod
