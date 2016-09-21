c     name_decl.h
c
c.......................................................................
c     This file will hold all declarations of the namelist variables
c     (given in name.h), type and dimensions.
c.......................................................................



c.......................................................................
c     nml variables that take on the values assigned to parameters.
c.......................................................................
      common /params/
     1  iy,jx,
     1  lfield,lz,lrz,
     1  lrzmax,ls,lsmax,
     1  mx,
     1  nbctime,negyrg,ngen,nmax

c.......................................................................
c     SCALAR INPUT FOR NAMELIST SETUP...
c.......................................................................

      character*8
     1  chang,
     1  eqmod,eleccomp,
     1  iactst,ineg,idrop,idskf,idskrf,ichkpnt,implct,
     1  lbdry0,locquas,lrzdiff,lsdiff,taunew,
     1  machine,meshy,manymat,
     1  netcdfvecal,netcdfvecc,netcdfvece,netcdfvecrf,netcdfvecs,
     1  noplots,
     1  psimodel,pltpowe,pltend,pltinput,pltlim,
     1  pltrst,plturfb,pltvflu,pltra,pltfvs,pltd,pltprpp,pltfofv,pltlos,
     1  pltdn,pltvecal,pltvecc,pltvecrf,pltvece,pltstrm,pltflux,
     1  pltsig,
     1  profpsi,
     1  qsineut,trapmod,scatmod,
     1  relativ,sigmamod,soln_method,
     1  symtrap,syncrad,bremsrad,
     1  tandem,gamafac,
     1  yreset,
     1  nmlstout

      character*256 netcdfnm
      character*256 mnemonic
      character*8 iuser,izeff,netcdfshort

      integer colmodl


      common /readscal/
     1  btor,bth,
     1  contrmin,constr,chang,colmodl,
     1  deltabdb,droptol,dtr,dtr0,
     1  esink,ephicc,esfac,eoved,enorm,enorme,enormi,eleccomp,
     1  eqmod,
     1  gsla,gslb,gamaset,gamafac,
     1  iactst,ineg,idrop,idskf,idskrf,ichkpnt,implct,
     1  isigtst,isigsgv1,isigsgv2,
     1  iuser,izeff,
     1  kenorm,lfil,
     1  lbdry0,locquas,lrzdiff,lsdiff,taunew,
     1  nnspec

      complex*16 vlfemin,vlfeplus

      common /readscal/
     1  machine,mnemonic,meshy,manymat,netcdfnm,netcdfshort,
     1  netcdfvecal,netcdfvecc,netcdfvece,netcdfvecrf,netcdfvecs,
     1  nstop,ncoef,nchec,ncont,nrstrt,nstps,nfpld,
     1  noncntrl,nonel,noffel,nonvphi,noffvphi,nonloss,noffloss,
     1  numby,lnwidth,noplots,nmlstout,
     1  psimodel,pltpowe,pltend,pltinput,pltlim,pltlimm,
     1  pltrst,plturfb,pltvflu,pltra,pltfvs,pltd,pltprpp,pltfofv,pltlos,
     1  pltdn,pltvecal,pltvecc,pltvecrf,pltvece,
     1  pltstrm,pltflux,pltmag,pltsig,profpsi,
     1  xprpmax,
     1  qsineut,trapmod,trapredc,scatmod,scatfrac,
     1  radmaj,radmin,rmirror,relativ,
     1  sigvi,sigmamod,sigvcx,soln_method,
     1  symtrap,syncrad,bremsrad,brfac,brfac1,brfacgm3,
     1  tfac,tfacz,tandem,temp_den,
     1  veclnth,
     1  vnorm

      real*8 mpwr,npwr,megy,negy,mtorloss,ntorloss,
     1  mpwrzeff,npwrzeff,mpwrvphi,npwrvphi,
     1  mpwrelec,npwrelec,mpwrxj,npwrxj

      common /readscal/
     1  xfac,xpctlwr,xpctmdl,xlwr,xmdl,xsinkm,
     1  yreset,ylower,yupper,
     1  mpwrzeff,npwrzeff,mpwrvphi,npwrvphi,mpwrxj,npwrxj,
     1  mpwrelec,npwrelec,
     1  elecscal,enescal,zeffscal,vphiscal,tescal,tiscal

c..................................................................
c     VECTOR DIMENSIONED NAMELIST SETUP COMMON BLOCK.....
c..................................................................

      character*8
     1  ibox,
     1  kpress,kfield,
     1  lbdry,lossmode,
     1  regy,
     1  torloss,
     1  difus_type

      character*256  lossfile

      common /readvec/
     1  bnumb(ntotala),
     1  fmass(ntotala),
     1  ibox(3),ioutput(2),isigmas(6),pltflux1(7),
     1  kpress(ntotala),kfield(ntotala),
     1  mpwr(0:ntotala),megy(ngena),mtorloss(ngena),
     1  npwr(0:ntotala),negy(ngena),ntorloss(ngena),
     1  lbdry(ngena),
     1  enloss(ngena),lossfile(ngena),
     1  lossmode(ngena),regy(ngena),gamegy(ngena),
     1  torloss(ngena),paregy(ngena),peregy(ngena),pegy(ngena),
     1  zeffin(0:njenea),vphiplin(0:njenea),
     1  ennl(npaproca),ennb(npaproca),ennscal(npaproca),
     1  bctime(nbctimea),dtr1(ndtr1a),nondtr1(ndtr1a),
     1  zeffc(nbctimea),zeffb(nbctimea),
     1  elecc(nbctimea),elecb(nbctimea),
     1  vphic(nbctimea),vphib(nbctimea),
     1  xjc(nbctimea),xjb(nbctimea),
     1  totcrt(nbctimea),
     1  nplot(nplota),
     1  difus_type(ngena)

c..................................................................
c     TWO DIMENSIONAL NAMELIST SETUP COMMON BLOCK.....
c..................................................................

      character*8
     1  kspeci

      common /readarr/
     1  tauloss(3,ngena),
     1  kspeci(2,ntotala),
     1  fpld(10,ngena),
     1  redenc(nbctimea,ntotala),redenb(nbctimea,ntotala),
     1  tempc(nbctimea,ntotala),tempb(nbctimea,ntotala)



c**********************************************************************
c     Variables in common block diskx have as their last dimension lrza.
c     Thus they are dimensioned with respect to the radial coordinate.
c     If lrza=1  ==> single flux surface CQL run.
c     Note that input variable lrz can be less than lrza. For some of
c     the larger arrays we allocate space using lrz rather than
c     dimensioning with lrza to save some space.
c
c     VECTORS
c
c**********************************************************************
c
      common /diskx/
     1  elecfld(0:lrza),
     1  tbnd(lrorsa),
     1  rovera(lrza),
     1  zmax(0:lrza),
     1  lrindx(0:lrorsa),lsindx(0:lrorsa)

c..................................................................
c     2-D ARRAYS
c..................................................................

      common /diskx/
     1  reden(ntotala,0:lrza),
     1  temp(ntotala,0:lrza)

      common /diskx/
     1  tauegy(ngena,0:lrza),eparc(ngena,0:lrza),
     1  eperc(ngena,0:lrza),
     1  isoucof,faccof

c..................................................................
c     3-D ARRAYS
c..................................................................

c..................................................................
c     4-D ARRAYS
c..................................................................

      common /diskx/
     1  eegy(negyrga,2,ngena,lrza),
     1  jegy(negyrga,2,ngena,lrza)

c*****************************************************************
c     BEGIN arrays for analytic ion source (sou..) routines
c*****************************************************************

      real*8 mpwrsou,npwrsou

      common /diskx/
     1  mpwrsou(0:ngena),npwrsou(0:ngena)

      common /diskx/
     1  asor(ngena,nsoa,lrza)

      common /params/
     1  nso

      character*8
     1  pltso,
     1  soucoord,
     1  knockon,komodel,
     1  flemodel

      common /readscal/
     1  nsou,
     1  pltso,
     1  soucoord,
     1  knockon,komodel,nkorfn,nonko,noffko,soffvte,soffpr,
     1  flemodel,jfl,xlfac,xlpctlwr,xlpctmdl,xllwr,xlmdl

      common /readarr/
     1  nonso(ngena,nsoa),noffso(ngena,nsoa),
     1  sellm1(ngena,nsoa),sellm2(ngena,nsoa),seppm1(ngena,nsoa),
     1  seppm2(ngena,nsoa),sem1(ngena,nsoa),sem2(ngena,nsoa),
     1  sthm1(ngena,nsoa),scm2(ngena,nsoa),szm1(ngena,nsoa),
     1  szm2(ngena,nsoa)


c*****************************************************************
c     BEGIN arrays for rf package..(rf...,vlh[B,...,vlf...) routines
c*****************************************************************
      character*8
     1  urfmod,
     1  vlfmod,vlfbes,vlfnpvar,vlhmod,vprprop,vlhplse,vlhprprp,
     1  rfread,rdcmod,rdc_clipping
      character*256 rffile

      common /params/
     1  nrf
      common/readscal/ urfmod,
     1  vlhmod,vlhmodes,vprprop,vdalp,vlh_karney,
     1  vlhplse,vlhpon,vlhpoff,
     1  vlfmod,vlfmodes,
     1  vlfbes,vlfnpvar,rdc_upar_sign,
     1  rfread,rdcmod,rdc_clipping

      common /readvec/
     1  nonrf(ngena),noffrf(ngena),
     1  dlndau(nmodsa),vparmin(nmodsa),vparmax(nmodsa),
     1  vlhprprp(nmodsa),vprpmin(nmodsa),vprpmax(nmodsa),
     1  vlhpolmn(nmodsa),vlhpolmx(nmodsa),
     1  vlffreq(nmodsa),vlfharms(nmodsa),vlfharm1(nmodsa),
     1  vlfnp(nmodsa),vlfdnp(nmodsa),vlfddnp(nmodsa),
     1  vlfeplus(nmodsa),vlfemin(nmodsa),
     1  vlfpol(nmodsa),vlfdpol(nmodsa),vlfddpol(nmodsa),
     1  vlfnperp(nmodsa),vlfdnorm(nmodsa),
     1  vlfparmn(nmodsa),vlfparmx(nmodsa),
     1  vlfprpmn(nmodsa),vlfprpmx(nmodsa),
     1  rffile(nmodsa)



c..................................................................
c     arrays in input
c..................................................................

      common/arr3d/
     1  acoefne(4),acoefte(4),ryain(njenea),elecin(njenea),
     1  enein(njenea,ntotala),tein(njenea),tiin(njenea),
     1  ennin(njenea,npaproca),irzplt(lrorsa),
     1  rya(0:lrza+1),
     1  thet1(nva),thet2(nva),thet1_npa(nva),thet2_npa(nva),
     1  nplt3d(nplota)

      common/arr3d/
     1  sellm1z(ngena,nsoa,0:lrza),sellm2z(ngena,nsoa,0:lrza),
     1  seppm1z(ngena,nsoa,0:lrza),sem1z(ngena,nsoa,0:lrza),
     1  sem2z(ngena,nsoa,0:lrza),sthm1z(ngena,nsoa,0:lrza),
     1  scm2z(ngena,nsoa,0:lrza),szm1z(ngena,nsoa,0:lrza),
     1  seppm2z(ngena,nsoa,0:lrza),
     1  szm2z(ngena,nsoa,0:lrza),asorz(ngena,nsoa,0:lrza)


c****************************************************************
c     BEGIN arrays for 3-d (td..) driver.
c****************************************************************


c..................................................................
c     scalars in input
c..................................................................

      character*8
     1  bootst,bootcalc,bootupdt,
     1  iprone,iprote,iproti,iprozeff,iprovphi,iproelec,ipronn,
     1  tmdmeth,partner,pinch,plt3d,pltvs,radcoord,
     1  relaxtsp,rzset,
     1  ndeltarho,softxry,npa_diag,atten_npa,
     1  transp,adimeth,
     1  efswtch,efswtchn,efiter,efflag,npa_process

      common /s3d/
     1  difusr,advectr,
     1  enmin,enmax,bootst,bootcalc,bootupdt,bootsign,pinch,
     1  fds,
     1  iprone,iprote,iproti,iprozeff,iprovphi,iproelec,ipronn,
     1  tmdmeth,kfrsou,
     1  mmsv,msxr,njene,njte,njti,nonboot,jhirsh,
     1  nrskip,nen,nv,nen_npa,nv_npa,npaproc,
     1  nr_delta,nz_delta,nt_delta,
     1  plt3d,pltvs,partner,radcoord,
     1  rfacz,rzset,roveram,relaxden,relaxtsp,
     1  ndeltarho,softxry,npa_diag,
     1  transp,
     1  enmin_npa,enmax_npa,fds_npa,
     1  atten_npa,adimeth,nonadi,
     1  efswtch,efswtchn,efiter,efflag,
     1  curr_edge,efrelax,efrelax1,currerr

      common/ar3d/
     1  difus_rshape(8),difus_vshape(4),difin(njenea),
     1  rd(nva),thetd(nva),x_sxr(nva),z_sxr(nva),
     1  rd_npa(nva),thetd_npa(nva),x_npa(nva),z_npa(nva),
     1  npa_process(npaproca)



c******************************************************************
c     BEGIN arrays for EQUILIBRIUM MODEL (eq..) (NON-CIRCULAR CROSS
c     SECTIONS).
c******************************************************************

      character*8 nconteq

      common/params/ nconteq,nconteqn

      character*8
     1  eqsym,eqdskalt,eqsource,eqmodel,
     1  fpsimodl

c     ONETWO uses character*60 for eqdskin.
      character*256 eqdskin

      common/readscal/
     1  atol,
     1  ellptcty,eqmodel,eqpower,eqsource,eqdskin,bsign,
     1  eqsym,eqdskalt,
     1  fpsimodl,
     1  methflag,
     1  povdelp,
     1  rtol,rmag,rbox,rboxdst,
     1  zbox


c*********************************************************************
c     BEGIN arrays for LOWER HYBRID FAST WAVE and ECH Module.
c*********************************************************************


      common/params/
     1  nurftime

      character*8
     1  urfdmp,iurfcoll,iurfl,
     1  call_lh,
     1  call_ech,
     1  call_fw,
     1  ech,fw,lh,
     1  scaleurf,urfrstrt,urfwrray,
     1  rftype

      common/readscal/
     1  call_lh,call_ech,call_fw,ieqbrurf,urfncoef,
     1  ech,
     1  fw,
     1  lh,
     1  nmods,
     1  nbssltbl,nondamp,nrfstep2,
     1  nrfpwr,nrfitr1,nrfitr2,nrfitr3,
     1  scaleurf,urfrstrt,urfwrray,
     1  urfdmp,urfmult

      common/readvec/
     1  pwrscale(nmodsa),wdscale(nmodsa),nrfstep1(nmodsa),
     1  nharms(nmodsa),nharm1(nmodsa),nrfspecies(nmodsa),
     1  iurfcoll(nmodsa),iurfl(nmodsa),rftype(nmodsa)

      common/readvec/
     1  pwrscale1(nbctimea),urftime(nbctimea)

c-----------------------------------------------------------------------
c     BEGIN variables for WP... modules for CQLP case
c-----------------------------------------------------------------------

      character*8  
     1  special_calls,cqlpmod,
     1  oldiag,
     1  sbdry,scheck,eseswtch,
     1  updown

c      logical
      character*8
     1  nlrestrt,nlwritf

      logical
     1  nlotp1,nlotp2,nlotp3,nlotp4

      common /readscal/
     1  special_calls,cqlpmod,
     1  epsthet,
     1  elpar0,
     1  lmidpln,lmidvel,laddbnd,
     1  nchgdy,ngauss,nlagran,
     1  nonavgf,nofavgf,
     1  nontran, nofftran, nonelpr, noffelpr,
     1  nlrestrt,nlwritf,nummods,numixts,
     1  oldiag,
     1  sbdry,scheck,eseswtch,
     1  updown

      common /readvec/
     1  denpar(ntotala,0:lza+1),
     1  nkconro(ntotala),nlotp1(noutpta),nlotp2(noutpta),
     1  nlotp3(noutpta),nlotp4(noutpta),
     1  temppar(ntotala,0:lza+1)



c.......................................................................
c     Setup block for finite orbit width (FOW) calculations
c.......................................................................

      character*8  fow
      character*16 outorb 
      character*38 file_fow_plt ! for saving data on orbit to a file
      
      common/fow_control/fow,outorb,file_fow_plt, nmu,npfi,
     +                  nsteps_orb,nptsorb,i_orb_width,iorb2,
     +  j0_ini,j0_end,inc_j0, i0_ini,i0_end,inc_i0,
     +  j2_ini,j2_end,inc_j2, i2_ini,i2_end,inc_i2
      
      ! fow= 'enabled' or 'disabled' 
      
      ! outorb  ! 'detailed' or 'Not-detailed'
                ! (saving/not-saving data to a file for plotting)  

      ! nmu     ! grid sizes for ad.ivariant mu 
      ! npfi    ! and canonical momentum Pfi; 
                ! to setup COM->R lookup table.

      ! nsteps_orb ! Max.number of time steps for orbit integration.
                ! Also used to trace Pfi=const levels for COM->R table
                ! in order to find intersections with mu=const levels.
      
      ! nptsorb ! Number of points on a complete orbit 
                ! (ityp=0 "main" orbit)
                ! from which ityp=1 "secondary" orbits are launched.
                ! ityp=1 orbit is stopped when it reaches the midplane.
                ! (Note: secondary orbits are not traced usually, 
                ! see below, iorb2=0)

      ! i_orb_width ! 1 -> Normal finite-orbit-width calculations. 
                    ! 0 -> V_drift_perp is set to 0 (ZOW approximation)
   	                     
      ! iorb2  ! set to 1 to perform Runge-Kutta integration for tracing
               ! SECONDARY orbits to midplane; 0 - no RK tracing.
               ! This option (1) can be used for plotting orbits
               ! (together with outorb="detailed"),
               ! otherwise not needed.
