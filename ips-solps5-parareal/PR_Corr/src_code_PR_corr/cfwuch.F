      subroutine cfwuch (nput, n, chfun, id)
      use b2mod_types
      implicit none
      integer nput, n
      character chfun*(*), id*(*)
*     (Please think of chfun as character*n. The declaration must be
*     character*(*) to allow the case n.eq.0.)
*     ------------------------------------------------------------------
*     CFWUCH writes a character array to an un*formatted file.
*     ------------------------------------------------------------------
      character chf*12, idcod*8, idtyp*8, id0*32
      intrinsic len
      logical streql
      external xertst, xerrab, streql
*     ------------------------------------------------------------------
*   ..preliminaries
      call xertst (n.le.len(chfun), 'cfwuch--faulty length for chfun')
      idcod = '*cf:'
      idtyp = 'char'
      id0 = id
      inquire (nput,form=chf,err=91)
*   ..branch on form
      if (streql(chf,'formatted')) then
*    ..formatted write
       write (nput,'(2a8,i12,4x,a32)') idcod, idtyp, n, id0
       if (0.lt.n) then
        write (nput,'(1x,a)') chfun(1:n)
       endif
      else if (streql(chf,'unformatted')) then
*    ..unformatted write
       write (nput) idcod, idtyp, n, id0
       if (0.lt.n) then
        write (nput) chfun(1:n)
       endif
      else
       call xerrab ('cfwuch--invalid file form-property: '//chf)
      endif
      return
*     ------------------------------------------------------------------
   91 call xerrab ('cfwuch--error return from inquire statement')
*     ------------------------------------------------------------------
      end
