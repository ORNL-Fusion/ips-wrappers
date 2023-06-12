c
c
c        Notes File.  (Try starting this as of June 21, 1994).
c
c.......................................................................
c     Things to watch for in old cqlinput NAMELIST FILES:
c     (BobH, 990702)
c.......................................................................
c
c
c      -To make namelist files maximally compatible
c       across platforms (Joe Freeman experience),
c       begin and end all namelist sections with "$",
c       (rather than "&" or "/".
c      -Replace all "-quotes with '-quotes.
c
c      -On CRAY J90/f90 (only), need to remove non-namelist data
c       such as "**" between lines, header data before
c       first namelist section, and comments following ";".
c      -Fix eegy(,,) ==> eegy(,,,)
c      -Fix asor(,)  ==> asor(,,)
c      -Add noplots="enabled" to first invocation of "setup"
c       namelist, until all the dead-end references to GRAFLIB
c       routines are completely replaced by PGPLOT
c      -Check nconteq,izeff and multiply are character*8
c       (integer values now put in through nconteqn,multiplyn).
c      -Check fpld(,) is type real
c      -Careful that nrfstep1(3),pwrscale(3),wdscale(3) 
c       not set if nmodsa.le.2,
c       nonrf(2),noffrf(2) not set if ngena=1.
c      -Set bsign=-1. in eqsetup namelist, if the toroidal field 
c       in the eqdsk is negative.
c      -Need to set gamaset to nonzero value for ion simulation.
c       (Upgrade this?).
c
c
c.......................................................................
c     Modifications to 32-bit code to give code running on Cray J90
c     64-bit environment:                       990701 (BobH)
c.......................................................................
c
c      -Changes are indicated in the code by comment lines.
c       c_cray, as compared to c_pc. In further detail:
c      -Change machinea=2==>1 in param.h, for cray
c      -Reverse use of dgbtrf,dgbtrs to sgbtrf,sgbtrs (can use
c       replace_sgbtrf script).
c      -malloc,free ==> hpalloc,hpdeallc
c      -Reverse use of drand(iflag) and ranf()
c      -Reverse use of call system and call ishell
c      -Need to set environmenatal variables 
c       PGPLOT_DIR=(pgplot directory) ,PGPLOT_DEV=/XWINDOW, e.g.
c      -makefile:remove references to r8lsode and urfpackm.
c      -cqlinput namelist file cannot have header of characters
c       between the namelist sections (such as "**").
c      -use netcdf.inc compatible with installed netCDF.
c      -note on netcdf: the NCDOUBLE designation for storage
c       of netcdf variables is netCDF based, rather than on
c       the architecture of the computer.
c       Thus, REAL*8 is NCDOUBLE, regardless.
c
c
c.......................................................................
c
c  f_code, indicating distribution functions calculated in the
c     code, is equal to f*vnorm**3,
c     where  f=distribution function in particles per cm**3, and
c              per velocity**3 (cgs).  For the relatavistic case
c              velocity is actually (momentum/rest_mass).
c            vnorm=maximum velocity of the grid (maximum 
c              momentum/rest_mass, for the relativistic case).
c
c
c
c  c..................................................................
c  c     generate radial (rho) mesh. rya will be the normalized mesh.
c  c     rrz will be the intermediate (mid) mesh points. Later rz
c  c     will be defined to be the non-normalized actual radial
c  c     toroidal flux mesh.
c  c..................................................................
c  c
c  c
c  c.......................................................................
c  c     determine array rovera
c  c.......................................................................
c
c        do 30 ll=1,lrzmax
c          rovera(ll)=rya(ll)
c          if (0.lt.rovera(ll) .and. rovera(ll).lt.1.e-8)
c       +    rovera(ll)=1.e-8
c   30   continue
c
c  c..................................................................
c  c     fill in input arrays between rya=0. and rya=1.
c  c     This determines density, temperature, source, etc profiles as
c  c     functions of the normalized radial coordinate rho. If density
c  c     and temperatures are to be specified as functions of bbpsi
c  c     poloidal flux coordinate these arrays will be overwritten
c  c     in subroutine eqcoord
c  c     If eqsource="tsc" (running with tsc):
c  c     iprone=iprote=iproti="disabled" (set in tdtscinp).
c  c..................................................................
c  c
c  c
c  c..................................................................
c  c..................................................................
c  c       cosz(1:iy,1:lzmax,1:lrzmax)  is cosine of pitch angle
c  c         corresponding to y(i,lmdlpln_), as given by cnsts of motion.
c  c..................................................................
c          if (ilcqlp .eq. 0) then
c            do 140 i=1,iiyh
c              cosz(i,l,lr_)=cvmgt((1.-bbpsi(l,lr_)*(sinn(i,lmdpln_)**2)),
c       1        1.e-90,1.-bbpsi(l,lr_)*(sinn(i,lmdpln_)**2) .gt. 1.e-11)
c              cosz(i,l,lr_)=sqrt(cosz(i,l,lr_))
c   140      continue
c          else
c            do 141 i=1,iiyh
c              cosz(i,l,lr_)=coss(i,indxls(l))
c   141      continue
c          endif
c  c..................................................................
c  c     imax(l,lr_) is the highest i such that a particle with a midplane
c  c     pitch angle y(i,lmdpln_) attains and passes z(l,lr_) before turning.
c  c     exception: imax(lz,lr_)=itl.
c  c..................................................................
c          do 150 i=1,iiyh
c            iii=iiy+1-i
c   150    cosz(iii,l,lr_)=-cosz(i,l,lr_)
c
c          imax(l,lr_)=1
c          if (ilcqlp .eq. 0) then
c            do 160 i=2,iyh_(lmdpln_)
c              if (cosz(i,l,lr_) .le. 1.e-10) go to 170
c   160      continue
c            imax(l,lr_)=iyh_(lmdpln_)
c            go to 180
c   170      imax(l,lr_)=i-1
c   180      continue
c          else
c            imax(l,lr_)=iiyh
c          endif
c
c  c..................................................................
c  c     sinz and tanz are defined analogously to cosz (above)
c  c..................................................................
c  c
c  c
c  c
c  c
c  c
c  c..................................................................
c  c     lmax(i,lr_) is the highest l such theta particle with midplane
c  c     pitch angle y(i,l_) attains and passes z(l,lr_). exception lmax(itl
c  c..................................................................
c
c  c
c  c
c  c
c  c
c  c
c  c
c  c..................................................................
c  c     zboun(i,lr_) is the point along a orbit at which a trapped
c  c     particle bounces. Passing particles are assigned zmax(lr_).
c  c     Particles at the passed trapped boundary are temporarily
c  c     assigned zstar. See subroutine baviorbt.
c  c..................................................................
c  c
c  c
c  c
c
c
c
c  There are (at least) 2 meshes along the magnetic field line:
c  (1)Subroutine eqfndpsi works with 2D arrays with arguments
c     l=1:lorbit(lr_),lr_=1,lrzmax.  lorbit(lr_) takes values 
c     up to lfielda (=250, presently), and gives a fine mesh of
c     points along a flux surface:
c        es(1:lorbit(lr_),lr_) is distance along the field line,
c          measured from the outer minimum B-field point (i.e., the
c          equatorial plane.
c        solr(1:lorbit(lr_),lr_) corresonding major radii
c        solz(             ,   )              verical height
c        eqbpol(               )   B_poloidal (T)
c        bpsi(                 )   B(es)/B(0.)
c        thtpol(               )   atan(solz/(rmag-solr))
c        eqdell(               )   delta(poloidal distance between
c                                  mesh points).
c  (2)Subroutine micxiniz works with a courser mesh along the field
c     line (requiring less storage and computations), but derived
c     from the above mesh:
c     z(1:lz,1:lrzmax) distance along field line,
c                      where  namelist lz.le.lza=80 presently.
c                      Spacing determined by tfacz.
c     pol(1:lz,lrzmax) from interpolation of thtpol at z( , )
c                           (psimodel.eq."spline")
c     bbpsi(           )                     bpsi at z( , )
c     solrz(  ,      )                       solr at z( , )
c     
