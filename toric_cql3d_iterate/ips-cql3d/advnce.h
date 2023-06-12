c     advnce.h
c***********************************************************************
c***********************************************************************

c..................................................................
c     advnce contains the statement functions utilized by splitting
c     or implicit time advancement routines and by their
c     diagnostic routines.
c..................................................................

      parameter(idima=2000)
c
      common/sge/
     1  nelt

      common/v2ge/
     1  avar(idima),
     1  ia(idima),
     1  ja(idima)


c..................................................................
c     Define some integration coefficients.
c..................................................................


      qz(j)=1./cint2(j)
      ry(i,j)=dx(j)*twopi/(cynt2(i,l_)*cint2(j))
      cl(i,j)=.25*vptb(itl,lr_)/vptb(i,lr_)*dc(i,j)
      r2y(j)=ry(itl,j)*.5
c%OS  dithta(i,j,l)=di(i,j,k,l)
      dithta(i,j,l_)=0.5-0.5*(1/(i+1))+0.5*(1/(iy_(l_)+1-i))

c###########################################################
c     Statement functions for explicit time advancement (splitting)
c     follow.
c###########################################################

c..................................................................
c     Using Chang-Cooper weights dj(i,j,k,l_) and di(i,j,k,l_) define
c     the weighted averages of the distribution function:
c     f1* - before split;    f2* - after split
c..................................................................

c      f1j(i,j)=temp1(i,j+1)*(1.-dj(i,j,k,l_))+temp1(i,j)*dj(i,j,k,l_)
c      f2j(i,j)=temp2(i,j+1)*(1.-dj(i,j,k,l_))+temp2(i,j)*dj(i,j,k,l_)
c      f1i(i,j)=temp1(i+1,j)*(1.-di(i,j,k,l_))+temp1(i,j)*di(i,j,k,l_)
c      f2i(i,j)=temp2(i+1,j)*(1.-di(i,j,k,l_))+temp2(i,j)*di(i,j,k,l_)
      ! YuP-101228: same as above, but re-arranged to have one '*'
      f1j(i,j)=temp1(i,j+1) + (temp1(i,j)-temp1(i,j+1))*dj(i,j,k,l_)
      f2j(i,j)=temp2(i,j+1) + (temp2(i,j)-temp2(i,j+1))*dj(i,j,k,l_)
      f1i(i,j)=temp1(i+1,j) + (temp1(i,j)-temp1(i+1,j))*di(i,j,k,l_)
      f2i(i,j)=temp2(i+1,j) + (temp2(i,j)-temp2(i+1,j))*di(i,j,k,l_)

c..................................................................
c     The first half of the splitting scheme consists of forward and
c     backward sweeps in the "x" or velocity direction. The forward
c     sweep requires the coefficients alpx betx gamx and the r.h.s.
c     delx of the equation:
c
c     -alpx(i,j)*f(i,j+1,l_) +betx(i,j)*f(i,j,l_) -gamx(i,j)*f(i,j-1,l_)
c
c     =delx(i,j)
c
c     Boundary conditions at x=0  (x=xmax)  that is j=1 (j=jx)
c     automatically force gamx(i,1) (alpx(i,jx)) to be zero.
c..................................................................

      alpx(i,j) = (da(i,j)*(1.-dj(i,j,k,l_))
     1  +db(i,j)*exp5(j)) * qz(j)

      betx(i,j) = (db(i,j)*exp5(j)
     1  +db(i,j-1)*exm5(j)
     1  -da(i,j)*dj(i,j,k,l_)
     1  +da(i,j-1)* (1.-dj(i,j-1,k,l_))) * qz(j)
     1  -vptb(i,lr_)*cah(i,j)
     1  +vptb(i,lr_)*rbgn

      gamx(i,j) = (-da(i,j-1)*dj(i,j-1,k,l_)
     1  +db(i,j-1)*exm5(j)) * qz(j)

c..................................................................
c     The quantity delx is complicated by the averaging done at the
c     pass/trapped boundary (i=itl).
c
c     Define the appropriate average at the p/t boundary
c..................................................................


      cdf(j) = (cl(itl-1,j)*(f1j(itl,j)-f1j(itl-1,j))
     1  *eyp5(itl-1,l_)
     1  +2.*cl(itl+1,j)*(f1j(itl+1,j)-f1j(itl,j))*eyp5(itl,l_)
     1  +cl(itu+1,j)*(f1j(itu+1,j)-f1j(itu,j))*eyp5(itu,l_))

      delx(i,j) = cvmgt(
     1  (dc(i,j)*(f1j(i+1,j)-f1j(i-1,j))*0.5*dyi(i,l_)
     1  -dc(i,j-1)*(f1j(i+1,j-1)-f1j(i-1,j-1))*0.5*dyi(i,l_)) *qz(j),
     1  (cdf(j)-cdf(j-1)) * qz(j)*ident(i),
     1  iota(i).ne.itl .and. iota(i).ne.itu)
     1  +vptb(i,lr_)*.5*so(i,j)
     1  +vptb(i,lr_)*temp1(i,j)*rbgn

c..................................................................
c     The sweep (y) requires the coefficients alpy, bety, gamy, and the
c     r.h.s. dely to the equation:
c
c     -alpy(i,j)*f(i+1,j,l_) +bety(i,j)*f(i,j,l_) -gamy(i,j)*f(i-1,j,l_)
c
c     = dely(i,j)
c
c     Boundary conditions at y=0 and y=pi automatically force gamy(1,j)
c     and alpy(iy,j) to be equal to zero.
c..................................................................

      alpy(i,j) = (dd(i,j)*(1.-di(i,j,k,l_))
     1  +df(i,j)*eyp5(i,l_)) *ry(i,j)

      bety(i,j) = (-dd(i,j)*di(i,j,k,l_)
     1  +df(i,j)*eyp5(i,l_)
     1  +dd(i-1,j)*(1.-di(i-1,j,k,l_))
     1  +df(i-1,j)*eym5(i,l_)) *ry(i,j)
     1  +vptb(i,lr_)*rbgn

      gamy(i,j) = -ry(i,j)*(dd(i-1,j)*di(i-1,j,k,l_)
     1  -df(i-1,j)*eym5(i,l_))

      dely(i,j) = ry(i,j)*0.5*dxi(j)*(de(i,j)*(f1i(i,ifp(j))-f1i(i,j-1))
     1  -de(i-1,j)*(f1i(i-1,ifp(j))-f1i(i-1,j-1)))
     1  +vptb(i,lr_)*rbgn*temp1(i,j)
     1  +vptb(i,lr_)*.5*so(i,j)

c..................................................................
c     Define the flux related quantities G and H the are used in
c     the r.h.s. of the Fokker-Planck equation.
c..................................................................
c-YuP: Moved gfu to diagentr.f, to avoid cvmgt() construct
c      gfu(i,j)=da(i,j)*f2j(i,j)
c     1  +db(i,j)*(temp2(i,j+1)-temp2(i,j))*exp5(j)
c     1  +cvmgt(  dc(i,j)*(f1j(i+1,j)-f1j(i-1,j))*0.5*dyi(i,l_),
c     1  cdf(j),
c     1  i .ne. itl  .and. i .ne. itu)

      hfu(i,j)=dd(i,j)*f2i(i,j)
     1  +de(i,j)*(f1i(i,ifp(j))-f1i(i,j-1))*0.5*dxi(j)
     1  +df(i,j)*(temp2(i+1,j)-temp2(i,j))*eyp5(i,l_)


c..................................................................
c     End splitting scheme time advancement statement functions...
c..................................................................


c######################################################
c     Begin implicit time advancement statement functions..
c######################################################

c..................................................................
c     advnce contains all of the parameters, arrays, and function
c     statements necessary to create the sparse matrix which represents
c     the implicit set of Fokker-Planck equations. The matrix is inverted
c     via Gaussian elimination (White routines ZSGBFA,ZSGBZL)

c
c     Below are statement functions which are used to determine
c     the matrix to be inverted. The i and j are the indices of the
c     equation considered - i being the theta index, j the velocity
c     index. On the left hand side the "m" stands for minus, the "p"
c     for plus and the "u" for upper pass/trapped mesh point itu and
c     "0" for neutral. The "x" stands for coefficients that hold
c     everywhere but the pass/trapped boundary and the "t" for
c     the equation at the pass/trapped boundary. For example,
c     the equation which represents the mesh point (i,j) involves
c     a sum of products of the distribution function on the left
c     hand side. The quantity which multiplies f(i-1,j,l_) is
c     xm0(i,j).
c
c     In terms of coefficients in Killeen et al.(1986) book
c     within a multiplicative constant:
c     da(i,j)=A_i,j+1/2      db(i,j)=B_i,j+1/2      dc(i,j)=C_i,j+1/2
c     dd(i,j)=D_i+1/2,j      de(i,j)=D+i+1/2,j      df(i,j)=F_i+1/2,j
c     di(i,j)=delta_i+1/2,j  dj(i,j)=delta_i,j+1/2
c
c.......................................................................

      xmm(i,j)=(-qz(j)*dc(i,j-1)*dj(i-1,j-1,k,l_)*dyi(i,l_)
     1  -ry(i,j)*de(i-1,j)*di(i-1,j-1,k,l_)*dxi(j))*.5

      x0m(i,j)=qz(j)*da(i,j-1)*dj(i,j-1,k,l_)+ry(i,j)*(de(i,j)*
     1  di(i,j-1,k,l_)-de(i-1,j)*(1.-di(i-1,j-1,k,l_)))*0.5*dxi(j)
     1  -qz(j)*db(i,j-1)*exm5(j)

      xpm(i,j)=(qz(j)*dc(i,j-1)*dj(i+1,j-1,k,l_)*dyi(i,l_)+
     1  ry(i,j)*de(i,j)*(1.-di(i,j-1,k,l_))*dxi(j))*.5

      xm0(i,j)=qz(j)*(dc(i,j)*dj(i-1,j,k,l_)
     1  -dc(i,j-1)*(1.-dj(i-1,j-1,k,l_)))
     1  *0.5*dyi(i,l_)+ry(i,j)*(dd(i-1,j)*di(i-1,j,k,l_)
     1  -df(i-1,j)*eym5(i,l_))
     1  +cthta(i,j)*dithta(i-1,j,l_)                      !Added since 1992

      x00(i,j)=qz(j)*
     1 (-da(i,j)*dj(i,j,k,l_)+da(i,j-1)*(1.-dj(i,j-1,k,l_))
     1  +db(i,j)*exp5(j)+db(i,j-1)*exm5(j))
     1  +ry(i,j)*(-dd(i,j)*di(i,j,k,l_)
     1  +dd(i-1,j)*(1.-di(i-1,j,k,l_))+df(i,j)*eyp5(i,l_)+df(i-1,j)
     1  *eym5(i,l_))
     1  -vptb(i,lr_)*(cah(i,j)-1./dtreff)
     1  +cthta(i,j)*(1.-dithta(i-1,j,l_)-dithta(i,j,l_))  !Added since 1992

      xp0(i,j)=qz(j)*(-dc(i,j)*dj(i+1,j,k,l_)
     1  +dc(i,j-1)*(1.-dj(i+1,j-1,k,l_)))*0.5*dyi(i,l_)
     1  -ry(i,j)*(dd(i,j)
     1  *(1.-di(i,j,k,l_))+df(i,j)*eyp5(i,l_))
     1  -cthta(i,j)*(1.-dithta(i,j,l_))                   !Added since 1992

      xmp(i,j)=qz(j)*dc(i,j)*(1.-dj(i-1,j,k,l_))*.5*dyi(i,l_)+
     1  ry(i,j)*de(i-1,j)*.5*dxi(j)*di(i-1,j+1,k,l_)

      x0p(i,j)=qz(j)*(-da(i,j)*(1.-dj(i,j,k,l_))-db(i,j)*exp5(j))
     1  +ry(i,j)*(-de(i,j)*di(i,j+1,k,l_)
     1  +de(i-1,j)*(1.-di(i-1,j+1,k,l_)))*0.5*dxi(j)

      xpp(i,j)=-qz(j)*dc(i,j)*(1.-dj(i+1,j,k,l_))*0.5*dyi(i,l_)
     1  -ry(i,j)*de(i,j)*(1.-di(i,j+1,k,l_))*0.5*dxi(j)

c.......................................................................
c     z00 is the right hand side of the equation, and holds the 
c     explicit-in-time rhs of the FP difference equations.
c
c     The terms involving the factors bsl, bsu , x**_ and t0**_ 
c     are related to calculation of the bootstrap effect.
c     We assume virtually that the distribution is skewed asymetrically
c     in the trapped region...that is we assume (virtually) that 
c     f(itl).ne.f(itu) and that the difference is driven by
c     a df/dr term through bsl and bsu. Since this term involves f at
c     different radial positions, it cannot figure into the solution 
c     implicitly, that is, it is differenced explicitly. The resulting 
c     contributions appear below. There will be contributions from 
c     i=itl-1, itu+1, itl and itu only.
c     All contributions are zero elsewhere, and are zero everywhere 
c     if bootcalc= "disabled".   (Refer to Harvey et al, 1993 Sherwood
c     Theory Mtg; E. Westerhof and A.G. Peters, Computer Physics Comm.,
c     Vol. 95, p. 131-138 (1996).)
c.......................................................................


ccc      z00f(i,j)=vptb(i,lr_)*(f_(i,j,k,l_)/dtreff+so(i,j)) +
ccc     +  spasou(i,j,k,l_)

c     itl-1 case:

ccc      z00itl1(i,j)=z00f(i,j)
ccc     1           -xpm(i,j)*bsl(j-1,k,l_)-xp0(i,j)*bsl(j,k,l_)
ccc     2           -xpp(i,j)*bsl(j+1,k,l_)

c     itu+1 case:

ccc      z00itu1(i,j)=z00f(i,j)
ccc     1           -xmm(i,j)*bsu(j-1,k,l_)-xm0(i,j)*bsu(j,k,l_)
ccc     2           -xmp(i,j)*bsu(j+1,k,l_)

c     itl or itu case:

ccc      t0ml_(j)=qz(j)*(
ccc     1 cl(itl-1,j-1)*dj(itl,j-1,k,l_)*eym5(itl,l_))
ccc     1 +r2y(j)*(-de(itl-1,j)*(1.-di(itl-1,j-1,k,l_)))
ccc     1 *0.5*dxi(j)

ccc      t00l_(j)=
ccc     1 +qz(j)*(
ccc     1 -cl(itl-1,j)*dj(itl,j,k,l_)*eym5(itl,l_)
ccc     1 +cl(itl-1,j-1)*(1.-dj(itl,j-1,k,l_))*eym5(itl,l_))
ccc     1 +r2y(j)*(dd(itl-1,j)*(1.-di(itl-1,j,k,l_))
ccc     1 +df(itl-1,j)*eym5(itl,l_))

ccc      t0pl_(j)=qz(j)*(
ccc     1 -cl(itl-1,j)*eym5(itl,l_)*(1.-dj(itl,j,k,l_)))
ccc     1 +r2y(j)*de(itl-1,j)*0.5*dxi(j)*(1.-di(itl-1,j+1,k,l_))


ccc      t0mu_(j)=qz(j)*(
ccc     1 -cl(itu+1,j-1)*dj(itu,j-1,k,l_)*eyp5(itu,l_))
ccc     1 +r2y(j)*(
ccc     1 +de(itu,j)*di(itu,j-1,k,l_))*0.5*dxi(j)

ccc      t00u_(j)=
ccc     1 +qz(j)*(
ccc     1 +cl(itu+1,j)*dj(itu,j,k,l_)*eyp5(itu,l_)
ccc     1 -cl(itu+1,j-1)*(1.-dj(itu,j-1,k,l_))*eyp5(itu,l_))
ccc     1 +r2y(j)*(
ccc     1 -dd(itu,j)
ccc     1 *di(itu,j,k,l_)
ccc     1 +df(itu,j)*eyp5(itu,l_))

ccc      t0pu_(j)=qz(j)*(
ccc     1 +cl(itu+1,j)*(1.-dj(itu,j,k,l_))*eyp5(itu,l_))
ccc     1 +r2y(j)*(-de(itu,j)*di(itu,j+1,k,l_)*0.5*dxi(j))

ccc      z00itl(i,j)=z00f(i,j)
ccc     1           -(t0ml_(j)*bsl(j-1,k,l_)+t00l_(j)*bsl(j,k,l_)
ccc     2           +t0pl_(j)*bsl(j+1,k,l_)+t0mu_(j)*bsu(j-1,k,l_)
ccc     3           +t00u_(j)*bsu(j,k,l_)+t0pu_(j)*bsu(j+1,k,l_))

ccc      z00ff(i,j)=cvmgt(z00itu1(i,j),z00itl(i,j),i.eq.(itu+1))

ccc      z00t(i,j)=cvmgt(z00itl1(i,j),z00ff(i,j),i.eq.(itl-1))

ccc      z00(i,j)=cvmgt(z00t(i,j),z00f(i,j),bootcalc.ne."disabled".and.
ccc     1                                  (i.eq.(itl-1).or.i.eq.itl.or.
ccc     2                                   i.eq.itu.or.i.eq.(itu+1)))

c.......................................................................
c     Pass/Trapped boundary statement functions follow..
c.......................................................................


      tmm(j)=-qz(j)*cl(itl-1,j-1)*dj(itl-1,j-1,k,l_)*eym5(itl,l_)
     1  -r2y(j)*di(itl-1,j-1,k,l_)*de(itl-1,j)*0.5*dxi(j)

      tm0(j)=qz(j)*cl(itl-1,j)*dj(itl-1,j,k,l_)*eym5(itl,l_)
     1  -qz(j)*cl(itl-1,j-1)*(1.-dj(itl-1,j-1,k,l_))*eym5(itl,l_)
     1  +r2y(j)*(dd(itl-1,j)*di(itl-1,j,k,l_)
     1  -df(itl-1,j)*eym5(itl,l_))

      tmp(j)=qz(j)*cl(itl-1,j)*(1.-dj(itl-1,j,k,l_))*eym5(itl,l_)
     1  +r2y(j)*di(itl-1,j+1,k,l_)*de(itl-1,j)*0.5*dxi(j)

      t0m(j)=qz(j)*(da(itl,j-1)*dj(itl,j-1,k,l_)-db(itl,j-1)*exm5(j)+
     1  cl(itl-1,j-1)*dj(itl,j-1,k,l_)*eym5(itl,l_)-2.*cl(itl+1,j-1)*
     1  dj(itl,j-1,k,l_)*eyp5(itl,l_)
     1  -cl(itu+1,j-1)*dj(itu,j-1,k,l_)*eyp5(itu,l_))
     1  +r2y(j)*(-de(itl-1,j)*(1.-di(itl-1,j-1,k,l_))+2.*de(itl,j)
     1  *di(itl,j-1,k,l_)+de(itu,j)*di(itu,j-1,k,l_))*0.5*dxi(j)

      t00(j)=vptb(itl,lr_)/dtreff
     1  +qz(j)*(-da(itl,j)*dj(itl,j,k,l_)+db(itl,j)*
     1  exp5(j)-cl(itl-1,j)*dj(itl,j,k,l_)*eym5(itl,l_)
     1  +2.*cl(itl+1,j)*dj(itl,j,k,l_)
     1  *eyp5(itl,l_)
     1  +cl(itu+1,j)*dj(itu,j,k,l_)*eyp5(itu,l_)
     1  +da(itl,j-1)*(1.-dj(itl,j-1,k,l_))
     1  +db(itl,j-1)*exm5(j)
     1  +cl(itl-1,j-1)*(1.-dj(itl,j-1,k,l_))*eym5(itl,l_)
     1  -2.*cl(itl+1,j-1)*(1.-dj(itl,j-1,k,l_))*eyp5(itl,l_)
     1  -cl(itu+1,j-1)*(1.-dj(itu,j-1,k,l_))*eyp5(itu,l_))
     1  +r2y(j)*(dd(itl-1,j)*(1.-di(itl-1,j,k,l_))
     1  +df(itl-1,j)*eym5(itl,l_)
     1  -2.*dd(itl,j)*di(itl,j,k,l_)
     1  +2.*df(itl,j)*eyp5(itl,l_)-dd(itu,j)
     1  *di(itu,j,k,l_)
     1  +df(itu,j)*eyp5(itu,l_))-vptb(itl,lr_)*cah(itl,j)

      t0p(j)=qz(j)*(-da(itl,j)*(1.-dj(itl,j,k,l_))-db(itl,j)*exp5(j)
     1  -cl(itl-1,j)*eym5(itl,l_)*(1.-dj(itl,j,k,l_))+2.*cl(itl+1,j)
     1  *eyp5(itl,l_)*(1.-dj(itl,j,k,l_))
     1  +cl(itu+1,j)*(1.-dj(itu,j,k,l_))*eyp5(itu,l_))
     1  +r2y(j)*(de(itl-1,j)*0.5*dxi(j)*(1.-di(itl-1,j+1,k,l_))-
     1  2.*de(itl,j)*di(itl,j+1,k,l_)*0.5*dxi(j)
     1  -de(itu,j)*di(itu,j+1,k,l_)*0.5*dxi(j) )

      tpm(j)=2.*qz(j)*cl(itl+1,j-1)*eyp5(itl,l_)*dj(itl+1,j-1,k,l_)
     1  +2.*r2y(j)*de(itl,j)*(1.-di(itl,j-1,k,l_))*0.5*dxi(j)

      tp0(j)=-2.*qz(j)*(cl(itl+1,j)*dj(itl+1,j,k,l_)*eyp5(itl,l_)-
     1  cl(itl+1,j-1)*(1.-dj(itl+1,j-1,k,l_))*eyp5(itl,l_))
     1  -2.*r2y(j)*df(itl,j)*eyp5(itl,l_)
     1  -2.*r2y(j)*dd(itl,j)*(1.-di(itl,j,k,l_))

      tpp(j)=-2*qz(j)*cl(itl+1,j)*eyp5(itl,l_)*(1.-dj(itl+1,j,k,l_))
     1  -2.*r2y(j)*de(itl,j)*0.5*dxi(j)*(1.-di(itl,j+1,k,l_))

      tum(j)=qz(j)*cl(itu+1,j-1)*eyp5(itu,l_)*dj(itu+1,j-1,k,l_)
     1  +r2y(j)*de(itu,j)*0.5*dxi(j)*(1.-di(itu,j-1,k,l_))

      tu0(j)=-qz(j)*cl(itu+1,j)*dj(itu+1,j,k,l_)*eyp5(itu,l_)
     1  +qz(j)*cl(itu+1,j-1)*eyp5(itu,l_)*(1.-dj(itu+1,j-1,k,l_))
     1  -r2y(j)*(dd(itu,j)*(1.-di(itu,j,k,l_))
     1  +df(itu,j)*eyp5(itu,l_))

      tup(j)=-qz(j)*cl(itu+1,j)*eyp5(itu,l_)*(1.-dj(itu+1,j,k,l_))
     1  -r2y(j)*de(itu,j)*(1.-di(itu,j+1,k,l_))*0.5*dxi(j)

c..................................................................
c     Express the Chang-Cooper weighted average f(i,j+1/2,l_): fpj
c..................................................................

      fpj(i,j)=f(i,j+1,k,l_)+ (f(i,j,k,l_)-f(i,j+1,k,l_))*dj(i,j,k,l_)

      fpjp(i,j)=fpj(i+1,j)
     1          +cvmgt(bsl(j,k,l_),zero,(i+1).eq.itl)
      fpj0(i,j)=fpj(i,j)
     1          +cvmgt(bsu(j,k,l_),zero,i.eq.itu) 
                 !YuP-110106: error corrected:  l_
     
c  Note: no need to check bootcalc="disabled" or not,
c  because when bootcalc="disabled",  bsl==0 and bsu==0.     

c..................................................................
c     Express the velocity flux at (i,j+1/2)
c..................................................................
c-YuP: Moved gft,gfi to diagentr.f, to avoid cvmgt() construct
c      gft(j)=
c     1  +cl(itl-1,j)*eyp5(itl-1,l_)*(fpjp(itl-1,j)-fpj(itl-1,j))
c     1  +2.*cl(itl+1,j)*eyp5(itl,l_)*(fpj(itl+1,j)-fpj(itl,j))
c     1  +cl(itu+1,j)*eyp5(itu,l_)*(fpj(itu+1,j)-fpj0(itu,j))
c
c      gfi(i,j)=da(i,j)*fpj(i,j)
c     1  +db(i,j)*exp5(j)*(f(i,j+1,k,l_)-f(i,j,k,l_))
c     1  +cvmgt(dc(i,j)*0.5*dyi(i,l_)*(fpjp(i,j)-fpj0(i-1,j)),
c     1  gft(j),
c     1  (i.ne.itl .and. i.ne.itu) .or. symtrap.ne."enabled")

c..................................................................
c     Express the Chang-Cooper weighted average f(i+1/2,j,l_): fpi
c..................................................................

      fpip(i,j)=f(i+1,j,k,l_)
     1          +cvmgt(bsl(j,k,l_),zero,(i+1).eq.itl)
      fpi0(i,j)=f(i,j,k,l_)
     1          +cvmgt(bsu(j,k,l_),zero,i.eq.itu)

c  Note: no need to check bootcalc="disabled" or not,
c  because when bootcalc="disabled",  bsl==0 and bsu==0.     

      fpi(i,j)=fpip(i,j)*(1.-di(i,j,k,l_)) + fpi0(i,j)*di(i,j,k,l_)

c..................................................................
c     Express the theta flux at (i+1/2,j)
c..................................................................

      hfi(i,j)=dd(i,j)*fpi(i,j)
     1  +de(i,j)*0.5*dxi(j)*(fpi(i,ifp(j))-fpi(i,j-1))
     1  +df(i,j)*eyp5(i,l_)*(fpip(i,j)-fpi0(i,j))

c..................................................................
c     End of statement functions used for implicit time advancement.
c..................................................................
