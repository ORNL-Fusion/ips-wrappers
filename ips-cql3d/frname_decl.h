c     frname_decl.h
c
c.......................................................................
c     This file will hold all declarations of the namelist variables
c     for the FREYA module given in frname.h, type and dimensions.
c.......................................................................

      common /numbrs_/   nprim,nimp

      character*8 ashape,bshape
      common /nub_/      anglev(kb),angleh(kb),naptr,ialign20
     &  ,ashape(nap,kb),aheigh(nap,kb),awidth(nap,kb),bcur(kb)
     &  ,alen(nap,kb),bleni(kb),blenp(kb)
     &  ,bhofset(kb),bvofset(kb),bptor(kb)
     &  ,bshape(kb),bheigh(kb),bwidth(kb)
     &  ,bhfoc(kb),bvfoc(kb)
     &  ,bhdiv(kb),bvdiv(kb)
     &  ,beamon,btime,ebkev(kb)
     &  ,fbcur(ke,kb)
     &  ,mf,nbeams,npart,npskip
     &  ,nsourc,ranseed,relnub
     &  ,sfrac1(kb),rpivot(kb),zpivot(kb)
     &  ,timbplt(5)
     &  ,fionx,hdepsmth


c................................................................

      character*8 frplt,frmod,multiply

      common /nub2_/ 
     &  ibcur,ibcx,iborb,ibslow,inubpat,npat(2),itrapfi,itrapech
c     ONETWO DIVERGENCE
     &  ,smooth,multiply,bmsprd,frmod,frplt,nfrplt,multiplyn


c................................................................

      common /nub3_/      iexcit,ilorent,mstate,kdene
     &  ,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh
     &  ,ngl
     &  ,izstrp(kimp)

c................................................................

      common /io_/  xdebug(20),nouthx

c................................................................

      character*8 nameb
      common /ions_/ nameb
