      subroutine cfrure (nget, n, refun, id)
      use b2mod_types
      implicit none
      integer nget, n
      real (kind=R8) ::   refun(0:n-1)
      character id*(*)
!*     ------------------------------------------------------------------
!*     CFRURE reads a real*8 array from an un*formatted file.
!*     ------------------------------------------------------------------
      integer n0, i
      character chf*12, idcod*8, idtyp*8, id0*32, id1*32
      logical streql
      external xertst, xerrab, streql
!*     ------------------------------------------------------------------
      id0 = id
!      inquire (nget,form=chf,err=91)
!*    ..formatted read
!      read (nget,'(2a8,i12,4x,a32)',err=92) idcod, idtyp, n0, id1
       if (0.lt.n) then
!     read (nget,*,err=93) (refun(i),i=0,n-1)
       endif
!*    ..unformatted read
      return
!*     ------------------------------------------------------------------
!   91 call xerrab ('cfrure--error return from inquire statement')
!   92 call xerrab ('cfrure--read error at head, formatted, id='//id0)
!   93 call xerrab ('cfrure--read error at data, formatted, id='//id0)
!   94 call xerrab ('cfrure--read error at head, unformatted, id='//id0)
!   95 call xerrab ('cfrure--read error at data, unformatted, id='//id0)
!*     ------------------------------------------------------------------
   92 write(6,*) 'cfrure--read error at head, formatted'
   93  write(6,*) 'cfrure--read error at data, formatted'  
   91  write(6,*) 'cfrure--error return from inquire statement'
      end
