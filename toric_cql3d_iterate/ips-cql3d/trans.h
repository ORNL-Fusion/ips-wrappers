c     trans.h

c..............................................................
c     If the transport model is utilized, define the alpha, beta
c     gamma and delta for the tridiagonal sweep option next.
c     h_r defined =H*rho=del_V/del_rho/(4*pi**2*R0)
c     For new soln_method="it3drv" option:
c     BH080307:
c     In alprp,betrp,gamrp, a zmaxpsi(l) factor is added
c     BH080410: 
c     vptb/zmaxpsi factor at appropriate l absorbed into each
c     alprp,betrp,gamrp.
c     BH171226:
c     For difus_io(k)= "drrin" and "drrdrin":
c     Added scale factor drrt(k) for d_rr and drt(k) for d_r.  These
c     scale factors are defaulted to 1d0, and set with function
c     difus_io_scale(,).  The d_rr (and d_r for difus_io(k).eq."drrdrin"
c     are not changed by the scale factor, rather the scale factors
c     are applied as appropriate.
c..............................................................

      ztr(i,l)=cosovb(idx(i,l),l)/dvol(lrindx(l))*4.*pi**2*radmaj
      ztrp(i,l)=cosovb(idx(i,l),l)/dvol(lrindx(l))*4.*pi**2*radmaj*
     1          zmaxpsi(lrindx(l))
      ytr(i,l)=h_r(lrindx(l))*drrt(k)*d_rr(idx(i,l),j,k,indxlr(l))/
     1  drp5(lrindx(l))*bovcos(idx(i,l),l)
c BH080410:  Changing sign of xtr (compatible with usual positive
c BH080410:  radial drift in the positive radial direction.
c            Need to check have changed sign for 
c            soln_method=direct, transp=enabled case.
c     The drt(k) scale factor is only applied for difus_io(k).eq.
c     "drrdrin"; else d_r is calculated from the scaled d_rr.
      xtr(i,l)=cvmgt(-h_r(lrindx(l))*drt(k)*d_r(idx(i,l),j,k,indxlr(l))*
     1  bovcos(idx(i,l),l),-h_r(lrindx(l))*d_r(idx(i,l),j,k,indxlr(l))*
     1  bovcos(idx(i,l),l),difus_io(k).eq."drrdrin")

      alpr(i,l)=ztr(i,l)*(ytr(i,l)+xtr(i,l)*(1.-dl(idx(i,l),j,k,l)))
      alprp(i,l)=ztrp(i,l)*(ytr(i,l)+xtr(i,l)*(1.-dl(idx(i,l),j,k,l)))*
     1  vptb_(idx(i,l+1),lrindx(l+1))*zmaxpsii(lrindx(l+1))
      betr(i,l)=ztr(i,l)*(ytr(i,l)-xtr(i,l)*dl(idx(i,l),j,k,l)+
     1  ytr(i,l-1)+xtr(i,l-1)*(1.-dl(idx(i,l-1),j,k,l-1)))+1./dttr
c BH071109:  Removing 1./dttr for soln_method="it3drv"
      betrp(i,l)=ztrp(i,l)*(ytr(i,l)-xtr(i,l)*dl(idx(i,l),j,k,l)+
     1  ytr(i,l-1)+xtr(i,l-1)*(1.-dl(idx(i,l-1),j,k,l-1)))*
     1  vptb_(idx(i,l),lrindx(l))*zmaxpsii(lrindx(l))
      gamr(i,l)=ztr(i,l)*(ytr(i,l-1)-xtr(i,l-1)*dl(idx(i,l-1),j,k,l-1))
      gamrp(i,l)=ztrp(i,l)*(ytr(i,l-1)-xtr(i,l-1)*
     1  dl(idx(i,l-1),j,k,l-1))*
     1  vptb_(idx(i,l-1),lrindx(l-1))*zmaxpsii(lrindx(l-1))
      delr(i,l)=frn(idx(i,l),j,k,l)/dttr*vptb_(idx(i,l),lrindx(l))*
     *  zmaxpsii(lrindx(l)) + velsou(idx(i,l),j,k,l)*
     *  zmaxpsii(lrindx(l))

      f1l(i,j,k,l)=frn(idx(i,l),j,k,l)*vptb_(idx(i,l),lrindx(l))*
     *  zmaxpsii(lrindx(l))*dl(idx(i,l),j,k,l)+
     +  frn(idx(i,l+1),j,k,l+1)*vptb_(idx(i,l+1),lrindx(l)+1)
     *  *zmaxpsii(lrindx(l)+1)*(1.-dl(idx(i,l),j,k,l))

      sfu(i,j,k,l)=drrt(k)*d_rr(idx(i,l),j,k,indxlr(l))
     1  *(frn(idx(i,l+1),j,k,l+1)
     1  *vptb_(idx(i,l+1),lrindx(l)+1)*zmaxpsii(lrindx(l)+1)-
     -  frn(idx(i,l),j,k,l)*
     *  vptb_(idx(i,l),lrindx(l))*zmaxpsii(lrindx(l)))/drp5(lrindx(l))+
     +  cvmgt(drt(k)*d_r(idx(i,l),j,k,indxlr(l))*f1l(i,j,k,l),
     ,  d_r(idx(i,l),j,k,indxlr(l))*f1l(i,j,k,l),
     ,  difus_io(k).eq."drrdrin")

cBH080425:  Changed sign on d_r, consistent above BH080410
      sfup(i,j,k,l)=drrt(k)*d_rr(idx(i,l),j,k,indxlr(l))*
     1  (frn(idx(i,l+1),j,k,l+1)*vptb_(idx(i,l+1),lrindx(l)+1)*
     1  zmaxpsii(lrindx(l)+1)-frn(idx(i,l),j,k,l)*
     1  vptb_(idx(i,l),lrindx(l))*zmaxpsii(lrindx(l)))/drp5(lrindx(l))-
     1  cvmgt(drt(k)*d_r(idx(i,l),j,k,indxlr(l))*f1l(i,j,k,l),
     1  d_r(idx(i,l),j,k,indxlr(l))*f1l(i,j,k,l),
     ,  difus_io(k).eq."drrdrin")
c*************************************************************

