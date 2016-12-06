      subroutine cfruch (nget, n, chfun, id)
      use b2mod_types
      implicit none
      integer nget, n
      character chfun*(*), id*(*)
*     (Please think of chfun as character*n. The declaration must be
*     character*(*) to allow the case n.eq.0.)
*     ------------------------------------------------------------------
*     CFRUCH reads a character array from an un*formatted file.
*     ------------------------------------------------------------------
      integer n0
      character chf*12, idcod*8, idtyp*8, id0*32, id1*32
      intrinsic len
      logical streql
      external xertst, xerrab, streql
*     ------------------------------------------------------------------
*   ..preliminaries
      call xertst (n.le.len(chfun), 'cfruch--faulty length for chfun')
      id0 = id
      inquire (nget,form=chf,err=91)
*   ..branch on form
      if (streql(chf,'formatted')) then
*    ..formatted read
       read (nget,'(2a8,i12,4x,a32)',err=92) idcod, idtyp, n0, id1
       call xertst (streql('*cf:',idcod),
     &   'cfruch--not a header line, expected *cf:, found '//idcod)
       call xertst (streql(id0,id1),
     &   'cfruch--wrong label: expected '//id0//', found '//id1)
       call xertst (streql('char',idtyp),
     &   'cfruch--wrong type: expected char, found '//idtyp)
       call xertst (n.eq.n0, 'cfruch--dimension mismatch at '//id0)
       if (0.lt.n) then
        read (nget,'(1x,a)',err=93) chfun(1:n)
       endif
      else if (streql(chf,'unformatted')) then
*    ..unformatted read
       read (nget,err=94) idcod, idtyp, n0, id1
       call xertst (streql('*cf:',idcod),
     &   'cfruch--not a header line, expected *cf:, found '//idcod)
       call xertst (streql(id0,id1),
     &   'cfruch--wrong label: expected '//id0//', found '//id1)
       call xertst (streql('char',idtyp),
     &   'cfruch--wrong type: expected char, found '//idtyp)
       call xertst (n.eq.n0, 'cfruch--dimension mismatch at '//id0)
       if (0.lt.n) then
        read (nget,err=95) chfun(1:n)
       endif
      else
       call xerrab ('cfruch--invalid file form-property: '//chf)
      endif
      return
*     ------------------------------------------------------------------
   91 call xerrab ('cfruch--error return from inquire statement')
   92 call xerrab ('cfruch--read error at head, formatted, id='//id0)
   93 call xerrab ('cfruch--read error at data, formatted, id='//id0)
   94 call xerrab ('cfruch--read error at head, unformatted, id='//id0)
   95 call xerrab ('cfruch--read error at data, unformatted, id='//id0)
*     ------------------------------------------------------------------
      end
