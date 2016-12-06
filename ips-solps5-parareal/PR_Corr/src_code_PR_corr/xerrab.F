      subroutine xerrab (messg)
      use b2mod_types
      implicit none
      character messg*(*)
*     ------------------------------------------------------------------
*     XERRAB writes an error message and the call chain, and terminates
*     execution of the program.
*     The output unit number that is used by xerrab may be controlled
*     by an earlier call to xerset; see the description of that routine.
*     The call chain is maintained by the group of routines (entry
*     points) in subsys. It is not required that these routines have
*     been invoked at any stage; if they have not been invoked, no call
*     chain will be written.
*     ------------------------------------------------------------------
      integer nprt
      logical lwf, lbot
      character fform*12, srname*8
      logical streql
      external xerget, streql, subxnd
*     ------------------------------------------------------------------
*     Note: Any routine that is called by xerrab must not call xerrab,
*     directly or indirectly. Routines xerget and streql are trivial
*     and should have no need for error handling, but routine subxnd
*     is more elaborate. It must not call xerrab.
*     ------------------------------------------------------------------
*   ..determine the output file
      call xerget (nprt)
      inquire (nprt,form=fform,err=1)
      lwf = streql(fform,'formatted')
      goto 2
    1 lwf = .false.
    2 continue
*     (lwf equivales "unit nprt is connected to a formatted file")
*   ..write error message
      if (lwf) then
       write (nprt,'(1x,a/1x,a)')
     &   '*** XERRAB: program will stop. ***', messg
      else
       write (*,'(1x,a/1x,a)')
     &   '*** XERRAB: program will stop. ***', messg
      endif
*   ..write the call chain
      call subxnd (lbot, srname)
      if (lwf) then
       if (lbot) then
        write (nprt,'(1x,a)') 'Call chain is not available.'
       else
        write (nprt,'(1x,a)') 'Call chain follows.'
       endif
      else
       if (lbot) then
        write (*,'(1x,a)') 'Call chain is not available.'
       else
        write (*,'(1x,a)') 'Call chain follows.'
       endif
      endif
    3 if (.not.lbot) then
       if (lwf) then
        write (nprt,'(2x,a)') srname
       else
        write (*,'(2x,a)') srname
       endif
       call subxnd (lbot, srname)
       goto 3
      endif
      stop 'xerrab--error'
*     ------------------------------------------------------------------
      end subroutine xerrab
