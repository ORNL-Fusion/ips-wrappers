c     name.h

c..................................................................
c     NAMELIST (SETUP0) DECLARATION FOR INPUT
c..................................................................

      namelist/setup0/mnemonic,ioutput,iuser,ibox,noplots,lnwidth,
     1  nmlstout,special_calls,cqlpmod,lrz,lrzdiff,lrzmax,lrindx,
     1  ls,lsmax,lsdiff,lsindx,nlrestrt,nlwritf

c..................................................................
c     NAMELIST (SETUP) DECLARATION FOR INPUT
c..................................................................

      namelist/setup/
     1  acoefne,acoefte,
     1  ampfmod,nampfmax,nonampf,ampferr,
     1  bnumb,btor,bth,bootst,bootcalc,bootupdt,bootsign,nonboot,jhirsh,
     1  contrmin,constr,chang,colmodl,
     1  deltabdb,denpar,droptol,dtr,dtr1,
     1  eegy,eparc,eperc,simpbfac,epar,eper,
     1  elecfld,elpar0,enorm,enorme,enormi,eleccomp,
     1  elecin,elecin_t,elecscal,enein,enein_t,ennin_t,
     +  enescal,enloss,epsthet,
     1  enmin,enmax,ennb,ennin,ennl,ennscal,enmin_npa,enmax_npa,
     1  eseswtch,xsink,esink,ephicc,esfac,eoved,
     1  fds,fds_npa,fmass,f4d_out,
     1  tavg,tavg1,tavg2,
     1  gsla,gslb,gamaset,gamafac,gamegy,
     1  iactst,ineg,idskf,idskrf,ichkpnt,implct,
     1  iprone,iprote,iproti,iprozeff,iprovphi,iproelec,ipronn,iprocur,
     1  tmdmeth,isigmas,isigtst,isigsgv1,isigsgv2,
     1  pltflux1,
     1  irzplt,izeff,ioutime,
     1  iy,
     1  jx,
     1  kenorm,lfil,kfrsou,kpress,kfield,kspeci,fpld,
     1  lmidpln,locquas,lbdry,lbdry0,lossfile,lossmode,lmidvel,laddbnd,
     1  lz,
     1  machine,meshy,manymat,netcdfnm,netcdfshort,
     1  netcdfvecal,netcdfvecc,netcdfvece,netcdfvecrf,netcdfvecs,
     1  nnspec,mpwr,megy,mtorloss,mmsv,msxr,mx,
     1  nchgdy,ngauss,nlagran,
     1  nlotp1,nlotp2,nlotp3,nlotp4,
     1  nmax,ngen,nkconro,nplt3d,nrskip,nen,nv,nen_npa,nv_npa,
     1  npaproc,npa_process,
     1  nr_delta,nz_delta,nt_delta,
     1  nr_f4d,nz_f4d,nv_f4d,nt_f4d,
     1  npwr,negy,ntorloss,njene,njte,njti,
     1  nstop,nondtr1,nplot,ncoef,nchec,ncont,nrstrt,nstps,nfpld,
     1  noncntrl,nonel,noffel,nonvphi,noffvphi,nonavgf,nofavgf,
     1  nonloss,noffloss,nummods,numixts,
     1  numby,negyrg,
     1  oldiag,
     1  plt3d,pltvs,partner,paregy,peregy,pegy,
     1  zeffin,zeffin_t,zeffscal,vphiplin,vphiplin_t,vphiscal,
     1  pltdn,pltvecal,pltvecc,pltvecrf,pltvece,
     1  pltstrm,pltflux,pltmag,pltsig,pltlim,pltlimm,
     1  pltrst,plturfb,pltvflu,pltra,pltfvs,pltd,pltprpp,pltfofv,pltlos,
     1  pltrdc,profpsi,
     1  psimodel,pltpowe,pltend,pltinput,pltview,
     1  qsineut,trapmod,trapredc,scatmod,scatfrac,
     1  ryain,radmaj,radmin,rmirror,relativ,
     1  reden,regy,rfacz,rzset,rd,roveram,
     1  rovera,rya,radcoord,
     1  sbdry,scheck,ndeltarho,softxry,npa_diag,symtrap,syncrad,
     1  bremsrad,brfac,brfac1,brfacgm3,sigmamod,sigvcx,sigvi,
     1  soln_method,tauegy,taunew,tein,tein_t,tescal,tiin,tiin_t,tiscal,
     1  tauloss,temp,temppar,
     1  tfac,tfacz,tbnd,tandem,
     1  thetd,torloss,thet1,thet2,x_sxr,z_sxr,
     1  rd_npa,thetd_npa,x_npa,z_npa,thet1_npa,thet2_npa,atten_npa,
     1  updown,
     1  veclnth,vnorm,
     1  xfac,xpctlwr,xpctmdl,xlwr,xmdl,xsinkm,
     1  xprpmax,ipxy,jpxy,
     1  yreset,ylower,yupper,
     1  mpwrzeff,npwrzeff,mpwrvphi,npwrvphi,mpwrxj,npwrxj,
     1  mpwrelec,npwrelec,
     1  redenc,redenb,temp_den,tempc,tempb,zeffc,zeffb,elecc,elecb,
     1  vphic,vphib,xjc,xjb,xjin_t,totcrt,efswtch,efswtchn,
     1  efiter,efflag,curr_edge,efrelax,efrelax1,currerr,
     1  bctime,nbctime,
     1  zmax,
     1  fow,outorb,nmu,npfi,nsteps_orb,nptsorb,i_orb_width,iorb2,
     1  j0_ini,j0_end,inc_j0, i0_ini,i0_end,inc_i0,
     1  j2_ini,j2_end,inc_j2, i2_ini,i2_end,inc_i2

c..................................................................
c     Namelist (trsetup) for "tr" transport routines
c..................................................................

      namelist/trsetup/
     1  advectr,
     1  difus_type,difusr,difus_rshape,difus_vshape,difin,
     1  difus_io,difus_io_file,
     1  difus_io_drrscale,difus_io_drscale,difus_io_t,
     1  pinch,
     1  relaxden,relaxtsp,
     1  transp,adimeth,nonadi,
     1  nontran,nofftran,nonelpr,noffelpr,ndifus_io_t


c..................................................................
c     Namelist (sousetup) for "sou" simple source routines.
c..................................................................

      namelist/sousetup/
     1  asorz,asor,flemodel,
     1  nonso,noffso,nso,nsou,
     1  pltso,mpwrsou,npwrsou,
     1  scm2z,szm1z,scm2,sellm1,sellm2,seppm1,
     1  sellm1z,sellm2z,seppm2,sem1,sem2,
     1  seppm1z,sem1z,sem2z,sthm1z,
     1  seppm2z,soucoord,knockon,komodel,nkorfn,nonko,noffko,soffvte,
     1  soffpr,isoucof,faccof,jfl,xlfac,xlpctlwr,xlpctmdl,xllwr,xlmdl,
     1  szm2z,sthm1,szm1,szm2

c..................................................................
c     Namelist (eqsetup) for "eq"(equilibrium geometry) routines.
c..................................................................

      namelist/eqsetup/
     1  atol,
     1  ellptcty,eqmodel,eqpower,eqsource,eqdskin,bsign,eqmod,
     1  eqsym,eqdskalt,
     1  fpsimodl,
     1  lfield,
     1  methflag,
     1  nconteq,nconteqn,
     1  povdelp,
     1  rtol,rmag,rbox,rboxdst,
     1  zbox

c..................................................................
c      Namelist (rfsetup) for "urf", "vlh", "vlf", rdcmod  modules.
c..................................................................

      namelist/rfsetup/
     1  call_lh,call_ech,call_fw,ieqbrurf,urfncoef,
     1  dlndau,
     1  lh,
     1  ech, 
     1  fw,
     1  rftype,
     1  rffile,rdcfile,rfread,
     1  nharms,nharm1,nrfspecies,iurfcoll,iurfl,
     1  nbssltbl,nondamp,nrfstep1,nrfstep2,
     1  nrfpwr,nrfitr1,nrfitr2,nrfitr3,
     1  nonrf,noffrf,nrf,
     1  scaleurf,pwrscale,wdscale,urfrstrt,urfwrray,
     1  nurftime,urftime,pwrscale1,
     1  urfdmp,urfmult,urfmod,
     1  vlhmod,vlhmodes,vparmin,vparmax,vprprop,vdalp,vlh_karney,
     1  vlhprprp,vlhplse,vlhpon,vlhpoff,vprpmin,vprpmax,vlhpolmn,
     1  vlhpolmx,vlfmod,vlfmodes,vlffreq,vlfnp,vlfdnp,vlfddnp,
     1  vlfeplus,vlfemin,vlfpol,vlfdpol,vlfddpol,vlfbes,vlfnpvar,
     1  vlfharms,vlfharm1,vlfnperp,vlfdnorm,
     1  vlfparmn,vlfparmx,vlfprpmn,vlfprpmx,rdc_upar_sign,nrdc,
     1  rdcmod,rdc_clipping,nrdcspecies,rdcscale,rdc_netcdf


     
   	               
