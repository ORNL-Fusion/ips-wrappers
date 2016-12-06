      subroutine xertst (lcond, messg)
      use b2mod_types
      implicit none
      logical lcond
      character messg*(*)
*     ------------------------------------------------------------------
*     XERTST tests a condition; if violated, calls xerrab to terminate
*     execution.
*     ------------------------------------------------------------------
      external xerrab
*     ------------------------------------------------------------------
      if (.not.lcond) then
       call xerrab (messg)
      endif
      return
*     ------------------------------------------------------------------
      end
