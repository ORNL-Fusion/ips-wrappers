      subroutine cfwure (nput, n, refun, id)
      use b2mod_types
      implicit none
      integer nput, n
      real (kind=R8) ::
     *   refun(0:n-1)
      character id*(*)
*     ------------------------------------------------------------------
*     CFWURE writes a real*8 array to an un*formatted file.
*!!   Note: If formatted output is selected, real*8 values must be
*     representable using three digits in the exponent.
*     ------------------------------------------------------------------
      integer ifmt, i
      parameter (ifmt=1)
      character fmt(0:5)*32, chf*12, idcod*8, idtyp*8, id0*32
      save fmt
      logical streql
      external xerrab, streql
      data fmt /'(1p,5(1x,e25.17e3))', '(1p,6(1x,e21.13e3))',
     &  '(1p,7(1x,e17.9e3))', '(1p,8(1x,e15.7e3))',
     &  '(1p,9(1x,e13.5e3))', '(1p,10(1x,e12.4e3))'/
*     ------------------------------------------------------------------
*!!   Note: set ifmt to select the desired format.
*     ------------------------------------------------------------------
*   ..preliminaries
      idcod = '*cf:'
      idtyp = 'real'
      id0 = id
      inquire (nput,form=chf,err=91)
*   ..branch on form
      if (streql(chf,'formatted')) then
*    ..formatted write
       write (nput,'(2a8,i12,4x,a32)') idcod, idtyp, n, id0
       if (0.lt.n) then
        write (nput,fmt(ifmt)) (refun(i),i=0,n-1)
       endif
      else if (streql(chf,'unformatted')) then
*    ..unformatted write
       write (nput) idcod, idtyp, n, id0
       if (0.lt.n) then
        write (nput) (refun(i),i=0,n-1)
       endif
      else
       call xerrab ('cfwure--invalid file form-property: '//chf)
      endif
      return
*     ------------------------------------------------------------------
   91 call xerrab ('cfwure--error return from inquire statement')
*     ------------------------------------------------------------------
      end
