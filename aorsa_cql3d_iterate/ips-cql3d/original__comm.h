c     comm.h
c
c.......................................................................
c     There are a number of arrays intended for temporary storage
c     in a subroutine, or for simple passes of data to
c     associated subroutines:
c
c     tem1(iyjx2a) --> tem6(iyjx2a)
c     temp1(0:iyp1a,0,jxp1a) --> temp6(0:iyp1a,0:jxp1a), 
c      "equivalenced" to tem1(iyjxa) --> tem6(iyjxa).
c     item1(iyjx2a) --> item6(iyjx2a)
c      "equivalenced" to tem1(iyjx2a) --> tem6(iyjx2a).
c     iyjx2a=(iya+2)*(jxa+2)
c     MOREOVER: We assume temp[1-6] are in contiguous storage,
c               so we can reference the whole six arrays through tem1.
c     tam1(jxa) --> tam30(jxa)
c     temc1(iya) --> temc4(iya)
c     tz1(lza) --> tz2(lza)
c     tr(0:lrza)
c     tr1(0:lrza) --> tr5(0:lrza)
c     itemc1(iya) --> itemc2(iya)
c     urftmp(nrayelta*nraya*2)
c.......................................................................
c
c
c.......................................................................
c     Add in type,size,common declarations for namelist variables.
c     This needs to precede the rest of comm.h declarations so
c     that name_decl.h can be used by itself in aindflt.f.
c     aindflt.f will also be used for setting defaults in the
c     SWIM project Integrated Plasma Simulation (IPS) modules.
c.......................................................................

      include 'name_decl.h'

c.......................................................................
c     variables that take on the values assigned to parameters.
c.......................................................................

      common /params/
     1  idim,iyjx,
     1  lrors,
     1  mxp1,mbet,
     1  nonch,niong,nionm,ntotal


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
     1  bmod0(lrza),btor0(lrza),bthr(lrza),btoru(lrza),bthr0(lrza),
     1  consn(lrorsa),consn0(lrorsa),currt(lrza),currxj0(0:lrza),
     1  currtp(lrza),currmtp(lrorsa),currmt(lrorsa),currxj(0:lrza),
     1  currpar(lrza),curreq(lrza),
     1  zreshin(lrza),zreskim(lrza),rovsc(lrza),rovsc_hi(lrza),
     1  curtor(lrza),curpol(lrza),ccurtor(0:lrza),ccurpol(0:lrza),
     1  deltapsi(lrza),
     1  eps(lrza),etll0(lrza),
     1  jchang(lrorsa),
     1  itl_(lrorsa),itu_(lrorsa),
     1  iy_(lrorsa),iyh_(lrorsa),iyjx_(lrorsa),
     1  indxlr(0:lrorsa),indxls(0:lrorsa),
     1  lorbit(lrza),lmdpln(0:lrza),
     1  n_(lrorsa),nch(lrorsa),
     1  psimx(lrza),pibzmax(lrza),psidz(lrza),
     1  qsafety(lrza),
     1  r0geom(lrza), r0drdz(0:lrza),rgeom(lrza), zgeom(lrza),
     1  rovs(lrza),rovsn(lrza),rovsloc(lrorsa)
      common /diskx/
     1  sptzr(lrorsa),sgaint1(lrorsa),starnue(lrorsa),
     1  thb(lrorsa),tauee(lrza),taueeh(lrorsa),time_(lrorsa),
     1  twoint(lrorsa),
     1  vthe(lrza),
     1  xlbnd(lrza),xlndn0(lrza),
     1  zmaxpsi(0:lrza),zmaxi(0:lrza),
     1  zmaxpsii(0:lrza),zeff(lrza),zeff4(lrza),enn(lrza),
     1  vphipl(lrza),
     1  srckotot(lrza),elecr(lrza),
     1  tauii(lrza),tau_neo(lrza),drr_gs(lrza),
     1  rhol(lrza),rhol_pol(lrza)
      common /diskx/
     1  taubi(lrza),tau_neo_b(lrza),drr_gs_b(lrza),
     1  rhol_b(lrza),rhol_pol_b(lrza)

c..................................................................
c     2-D ARRAYS
c..................................................................

      common /diskx/
     1  energy(ntotala,lrza),
     1  vth(ntotala,lrza)


      common /diskx/
     1  currm(ngena,lrorsa),curr(ngena,lrza),
     1  hnis(ngena,lrza),
     1  ratio(ngena,lrza),
     1  wpar(ngena,lrza),wperp(ngena,lrza),
     1  xlndn00(ngena,lrza),xlncur(ngena,lrza),
     1  xlndn(ngena,lrza),energym(ngena,lrorsa),
     1  xlndnr(ntotala,lrza),energyr(ntotala,lrza),
     1  currr(ngena,lrza),xlndnv(ntotala,lrza),
     1  energyv(ntotala,lrza),currv_(ngena,lrza),eoe0(ngena,lrza),
     1  ucrit(ngena,lrza),jxcrit(ngena,lrza),
     1  denra(ngena,lrza),curra(ngena,lrza),
     1  fdenra(ngena,lrza),fcurra(ngena,lrza)
      common /diskx/
     1  pm(0:mxp1a,lrorsa)
      common /diskx/
     1  alm(0:mbeta,lrza),
     1  betta(0:mbeta,lrza),
     1  entrintr(ngena,-1:15)

c..................................................................
c     3-D ARRAYS
c..................................................................

      common /diskx/
     1  sgaint(8,ngena,lrorsa)
      common /diskx/
     1  entr(ngena,-1:15,lrorsa)

c..................................................................
c     4-D ARRAYS
c..................................................................

      common /diskx/
     1  pentr(noncha,ngena,-1:15,lrorsa)



c..................................................................
c     Next follow variables and arrays not dimensioned by lrza.
c..................................................................

c..................................................................
c     SCALARS.....
c..................................................................

      character*8
     1  outfile,
     1  prefix,
     1  frmodp

      character*512 t_

      common
     1  clight,charge,restmkev,
     1  cnorm,cnorm2,cnorm3,clite2,cnorm2i,cnormi,
     1  ergtkev,eps0,elect,eratio,eovsz,eovedd,
     1  em6,em8,em10,em12,em14,em37,em40,ep37,em90,ep90,
     1  ep100,em100,em120,em300,four,fourpi,fc,
     1  gszb,gsb,gsem,gsep,gszm,gszp,gseb,gszm2,gszp2,gszb2,gsla2,
     1  gslb2,gsnm,half,
     1  iyh,istate,itl,itu,
     1  impcoef,imprf,impelec,impcah,irstart,impadi,
     1  jxm1,jlwr,jmdl,
     1  kgnande,kelecg,kelecm,kelec,kionn,kiong(ngena),kiongg(ngena),
     1  kionm(nmaxa),
     1  l_,lr_,indxlr_,indxls_,lmdpln_,ls_,
     1  miyjx,
     1  frmodp                  !This variable is reset to frmod, 
                                !namelist declared in frname.h and
                                !passed to the comm.h related subroutines
                                !as an argument of subroutine frnfreya
      common
     1  ipad1,r0geomp, rgeomp, zgeomp,
     1  niyjx,nccoef,ncentr,navnc,ncplt,n,
     1  nflux,nframe,npltmx,
     1  outfile,ipad2,one,one_,
     1  prefix,
     1  pi,pio180,psir,proton,
     1  rgyro,resist,resistn,rovsf,rmp0,rovscf,rbgn,
     1  stopm,
     1  tei,timet,t_,three,third,sevenhun,two,twopi,
     1  vparl,vpovth,
     1  xconn,xmax,ipxy,jpxy,
     1  zero,zstar,
     1  iplot,nirzplt

c..................................................................
c     VECTORS.....
c..................................................................

      character*8
     1  text

      common
     1  ipad3,text(4),
     1  iota(0:750),realiota(0:750),
     1  i1p(2),itab(3),ipad4,tab(3),
     1  imsh(20),
     1  frbuf(1024),
     1  ws2d(2500),
     1  work(tlfld1a)
      common
     1  ekev(ngena),
     1  engain(ngena),
     1  tnorm(ngena)
      common
     1  fions(ntotala),
     1  gama0(ntotala),
     1  tz1(lza),tz2(lza),
     1  ix1(0:mxa),ix2(0:mxa),ix3(0:mxa),ix4(0:mxa),
     1  ix5(0:mxa),ix6(0:mxa),ix7(0:mxa),ix8(0:mxa),
     1  tom1(0:mxp1a),tom2(0:mxp1a),tom3(0:mxp1a),tom4(0:mxp1a)
      common
     1  wkbc(3*nbctimea),iopbc(2)

      character*8
     1  mplot

       common
     1  mplot(lrorsa),
     1  tplot(nplota)

c..................................................................
c     TWO-DIMENSIONAL ARRAYS.....
c..................................................................

      common
     1  xlndnz(ngenpa,negyrga)

      common
     1  gamt(ntotala,ntotala),gama(ntotala,ntotala),
     1  satioz2(ntotala,ntotala),satiom(ntotala,ntotala),
     1  cog(0:mxa,15),
     1  sgain(8,ngena),
     1  choose(0:mxa+2,0:mxa+2),fctrl(0:mxa+2)

c..................................................................
c     Allocatable arrays allocated in subroutine ainalloc
c..................................................................

CPTR>>>REPLACE PTR-AINALLOC
      common /dptr/ dumptr
      common /dptr/ fptr,f_ptr,fxspptr,goneptr,sourcptr,eyp5ptr,psiivptr
      common /dptr/ i0tranptr
      common /dptr/ calptr,cblptr,cclptr,cdlptr,celptr,cflptr
      common /dptr/ ealptr,eblptr,scalptr
      common /dptr/ dzptr,cetptr,cexptr,diptr,djptr,dym5ptr,dyp5ptr
      common /dptr/ idxptr,yptr,dyptr,yptbptr,cossptr,cynt2ptr,batotptr
      common /dptr/ lmaxptr,vpintptr
      common /dptr/ psibaptr,psisqptr,psicuptr,psiquptr
      common /dptr/ bavpdptr,psiirptr
      common /dptr/ sinnptr,ymidptr,tauptr,imaxptr,eym5ptr,bavdnptr
      common /dptr/ tannptr,vptbptr,zbounptr,polptr,solrzptr,solzzptr
      common /dptr/ thtabptr,zptr,psiptr
      common /dptr/ consnptr,ptimptr,rovspptr,sptzrptr,vderbptr,vpovptr
      common /dptr/ restpptr,restnpptr,pefldptr
      common /dptr/ esptr,bpsiptr
      common /dptr/ d2bpsptr,d2srzptr,d2szzptr,d2thtptr,d2esptr
      common /dptr/ thtpoptr,esfiptr,espsiptr,psiesptr,psifiptr
      common /dptr/ souppptr
      common /dptr/ coszptr,dtauptr,sinzptr,tanzptr,yzptr,totptr
      common /dptr/ waaptr,wbbptr
      common /dptr/ vfluxptr
      common /dptr/ favptr
      common /dptr/ sincoptr,syncaptr,syncdptr,tauloptr
      common /dptr/ efldnptr
      common /dptr/ egyloptr,denszptr,ssptr
      common /dptr/ dcoflegptr,dpcosptr,ssyptr,ssyyptr,ssyiptr,ssyyyptr
      common /dptr/ pcurrptr,pdensptr,pengyptr
      common /dptr/ pdraptr,pcraptr,pfdraptr,pfcraptr,pucriptr,peoe0ptr
      common /dptr/ psrcptr,peoedptr
      common /dptr/ cint2ptr
      common /dptr/ dxptr
      common /dptr/ ifpptr
      common /dptr/ sgptr,sgxptr,sgxxptr
      common /dptr/ shptr,shxptr,shxxptr,shxxxptr
      common /dptr/ tam1ptr,tam2ptr,tam3ptr,tam4ptr,tam5ptr
      common /dptr/ tam6ptr,tam7ptr,tam8ptr,tam9ptr,tam10ptr
      common /dptr/ tam11ptr,tam12ptr,tam13ptr
      common /dptr/ tam14ptr,tam15ptr,tam16ptr,tam17ptr
      common /dptr/ tam18ptr,tam19ptr,tam20ptr,tam21ptr,tam22ptr
      common /dptr/ tam23ptr,tam24ptr,tam25ptr,tam26ptr,tam27ptr
      common /dptr/ tam28ptr,tam29ptr,tam30ptr
      common /dptr/ xptr,xmidpptr,xiptr,xsqptr,x3iptr,x2iptr
      common /dptr/ xcuptr,xcentptr,xcensptr,uocptr,enerkptr
      common /dptr/ gammaptr,gamsqptr,gamcuptr,gammiptr,gamm2ptr
      common /dptr/ gamm1ptr,gamnptr,alpnptr,asnhptr,tcsgmptr
      common /dptr/ gamfacptr
      common /dptr/ dxm5ptr,dxp5ptr
      common /dptr/ exm5ptr,exp5ptr
      common /dptr/ identptr
      common /dptr/ temc1ptr,temc2ptr,temc3ptr,temc4ptr
      common /dptr/ itmc1ptr,itmc2ptr
      common /dptr/ l_lowptr,lptptr,munptr
      common /dptr/ fllptr,xparptr,rheadptr,dfvleptr,dfvliptr
      common /dptr/ flptr,fl1ptr,fl2ptr
      common /dptr/ xlptr,jmaxxptr,xlmptr,dxlptr,pparsptr,pprpsptr
      common /dptr/ faciptr,pareaptr,wtfl0ptr,wtflmptr,jflbnptr
      common /dptr/ xperpptr
      common /dptr/ xmptr
      common /dptr/ daptr,dbptr,dcptr,dbbptr
      common /dptr/ ddptr,deptr,dfptr,dffptr
      common /dptr/ cahptr,cthtaptr,gonptr
      common /dptr/ soptr, spasptr, velsptr, vels2ptr
      common /dptr/ temp1ptr,temp2ptr,temp3ptr,temp4ptr
      common /dptr/ temp5ptr,temp6ptr
      common /dptr/ xheadptr,yheadptr,xtailptr,ytailptr,fpnptr
      common /dpre/ xlljiptr,xppjiptr
      common /dptr/ egptr,fgptr
      common /dptr/ tem1ptr,tem2ptr,tem3ptr,tem4ptr,tem5ptr,tem6ptr
      common /dptr/item1ptr,item2ptr,item3ptr,item4ptr,item5ptr,item6ptr
      common /dptr/ captr,cbptr,ccptr,cdptr,ceptr,cfptr
      common /dptr/ currvptr,curr_ptr
      common /dptr/ pwrrfptr,pwrr_ptr
      common /dptr/ talptr,tblptr,tflptr
      common /dptr/ plegptr
      common /dptr/ fetaptr,fetbptr
      common /dptr/ wfluxptr
      common /dptr/ rhsptr
      common /dptr/ sovtptr
      common /dptr/ sigsxptr
      common /dptr/ tamt1ptr,tamt2ptr

      pointer(dumptr,dum(1))
      pointer(fptr,f(0:iyp1a,0:jxp1a,ngena,lrors))
      pointer(fxspptr,fxsp(0:iyp1a,0:jxp1a,ngena,lrors))
      pointer(f_ptr,f_(0:iyp1a,0:jxp1a,ngena,lrors))
      pointer(spasptr,spasou(0:iyp1a,0:jxp1a,ngena,lrors))

      pointer(velsptr,velsou(0:iyp1a,0:jxp1a,ngena,0:lrors+1))
      pointer(vels2ptr,velsou2(0:iyp1a,0:jxp1a,ngena,0:lrors+1))

      pointer(sourcptr,source(0:iyp1a,0:jxp1a,ngena,lrz))
      pointer(goneptr,gone(0:iyp1a,0:jxp1a,ngena,lrz))
      pointer(egyloptr,egylosa(0:iyp1a,0:jxp1a,ngena,lrz))

      pointer(i0tranptr,i0tran(i0param+1,lz,lrz))

      pointer(calptr,cal(iya,jxa,ngena,lrors))
      pointer(cblptr,cbl(iya,jxa,ngena,lrors))
      pointer(cclptr,ccl(iya,jxa,ngena,lrors))
      pointer(cdlptr,cdl(iya,jxa,ngena,lrors))
      pointer(celptr,cel(iya,jxa,ngena,lrors))
      pointer(cflptr,cfl(iya,jxa,ngena,lrors))
      pointer(ealptr,eal(iya,jxa,ngena,2,lrors))
      pointer(eblptr,ebl(iya,jxa,ngena,2,lrors))

      pointer(scalptr,scal(iyjxnga,lrors))

      pointer(cetptr,cet(iya,jxa,lrors))
      pointer(cexptr,cex(iya,jxa,lrors))

      pointer(syncaptr,synca(iya,jxa,lrz))
      pointer(syncdptr,syncd(iya,jxa,lrz))
      pointer(tauloptr,taulos(iya,jxa,lrz))

      pointer(efldnptr,elecfldn(nefitera,noncha,lrz))

      pointer(diptr,di(0:iya,0:jxp1a,lrors))

      pointer(djptr,dj(0:iyp1a,0:jxa,lrors))

      pointer(dym5ptr,dym5(0:iya,lrors))
      pointer(dyp5ptr,dyp5(0:iya,lrors))
      pointer(eyp5ptr,eyp5(0:iya,lrors))
      pointer(eym5ptr,eym5(0:iya,lrors))

      pointer(yptr,y(iya,lrors))
      pointer(dyptr,dy(iya,lrors))
      pointer(yptbptr,yptb(iya,lrors))
      pointer(cossptr,coss(iya,lrors))
      pointer(cynt2ptr,cynt2(iya,lrors))

      pointer(batotptr,batot(iya,lrzmax))
      pointer(lmaxptr,lmax(iya,lrzmax))
      pointer(vpintptr,vpint(iya,lrzmax))
      pointer(psiivptr,psiiv(iya,lrzmax))
      pointer(psibaptr,psiba(iya,lrzmax))
      pointer(psisqptr,psisq(iya,lrzmax))
      pointer(psicuptr,psicu(iya,lrzmax))
      pointer(psiquptr,psiqu(iya,lrzmax))
      pointer(bavpdptr,bavpd(iya,lrzmax))
      pointer(bavdnptr,bavdn(iya,lrzmax))
      pointer(psiirptr,psiir(iya,lrzmax))
      pointer(vderbptr,vderb(iya,lrzmax))

      pointer(sinnptr,sinn(iya,lrors))
      pointer(tannptr,tann(iya,lrors))
      pointer(ymidptr,ymid(iya,lrors))

      pointer(tauptr,tau(iya,lrzmax))
      pointer(vptbptr,vptb(iya,lrzmax))
      pointer(zbounptr,zboun(iya,lrzmax))

      pointer(idxptr,idx(iya,0:lrors))

      pointer(imaxptr,imax(lza,lrzmax))
      pointer(dzptr,dz(lza,lrzmax))
      pointer(polptr,pol(lza,lrzmax))
      pointer(solrzptr,solrz(lza,lrzmax))
      pointer(solzzptr,solzz(lza,lrzmax))
      pointer(thtabptr,thtab(lza,lrzmax))
      pointer(zptr,z(lza,lrzmax))
      pointer(psiptr,psi(lza,lrzmax))

      pointer(consnptr,consnp(noncha,lrors))
      pointer(ptimptr,ptime(noncha,lrors))
      pointer(sptzrptr,sptzrp(noncha,lrors))
      pointer(pefldptr,pefld(noncha,lrors))

      pointer(rovspptr,rovsp(noncha,lrorsa))

      pointer(restpptr,restp(noncha,lrzmax))
      pointer(restnpptr,restnp(noncha,lrzmax))
      pointer(vpovptr,vpov(noncha,lrzmax))

      pointer(esptr,es(lfielda,lrzmax))
      pointer(bpsiptr,bpsi(lfielda,lrzmax))
      pointer(d2bpsptr,d2bpsi(lfielda,lrzmax))
      pointer(d2srzptr,d2solrz(lfielda,lrzmax))
      pointer(d2szzptr,d2solzz(lfielda,lrzmax))
      pointer(d2thtptr,d2thtpol(lfielda,lrzmax))
      pointer(d2esptr,d2es(lfielda,lrzmax))
      pointer(thtpoptr,thtpol(lfielda,lrzmax))

      pointer(esfiptr,esfi(incza,lrzmax))
      pointer(psiesptr,psiesfi(incza,lrzmax))

      pointer(psifiptr,psifi(inczpa,lrzmax))
      pointer(espsiptr,espsifi(inczpa,lrzmax))

      pointer(souppptr,soupp(jxa,lrzmax))

      pointer(waaptr,waa(lza,0:mxa,lrzmax))
      pointer(wbbptr,wbb(lza,0:mxa,lrzmax))

      pointer(coszptr,cosz(iya,lza,lrzmax))
      pointer(dtauptr,dtau(iya,lza,lrzmax))
      pointer(sinzptr,sinz(iya,lza,lrzmax))
      pointer(tanzptr,tanz(iya,lza,lrzmax))
      pointer(yzptr,yz(iya,lza,lrzmax))
      pointer(totptr,tot(iya,lza,lrzmax))

      pointer(vfluxptr,vflux(jxa,ngena,lrors))

      pointer(favptr,f_aveth(jxa,ngena,lrors,5))

      pointer(sincoptr,sincosba(iya,ngena,lrzmax))

      pointer(denszptr,densz(lza,ngenpa,negyrga,lrzmax))

      pointer(ssptr,ss(iya,lza,0:mxp1a,lrzmax))

      pointer(dcoflegptr,dcofleg(iya,lza,0:mxa,lrz))

      pointer(dpcosptr,dpcosz(iya,lza,0:mxa,lrzmax))
      pointer(ssyptr,ssy(iya,lza,0:mxa,lrzmax))
      pointer(ssyyptr,ssyy(iya,lza,0:mxa,lrzmax))
      pointer(ssyiptr,ssyi(iya,lza,0:mxa,lrzmax))
      pointer(ssyyyptr,ssyyy(iya,lza,0:mxa,lrzmax))

      pointer(pcurrptr,pcurr(noncha,ngena,lrorsa))
      pointer(pdensptr,pdens(noncha,ngena,lrorsa))
      pointer(pengyptr,pengy(noncha,ngena,lrorsa))

      pointer(pdraptr,pdenra(noncha,lrz))
      pointer(pcraptr,pcurra(noncha,lrz))
      pointer(pfdraptr,pfdenra(noncha,lrz))
      pointer(pfcraptr,pfcurra(noncha,lrz))
      pointer(pucriptr,pucrit(noncha,lrz))
      pointer(peoe0ptr,peoe0(noncha,lrz))
      pointer(psrcptr,psrc(noncha,lrz))
      pointer(peoedptr,peoed(noncha,lrz))

      pointer(cint2ptr,cint2(jxa))
      pointer(dxptr,dx(jxa))
      pointer(ifpptr,ifp(jxa))
      pointer(sgptr,sg(jxa))
      pointer(sgxptr,sgx(jxa))
      pointer(sgxxptr,sgxx(jxa))
      pointer(shptr,sh(jxa))
      pointer(shxptr,shx(jxa))
      pointer(shxxptr,shxx(jxa))
      pointer(shxxxptr,shxxx(jxa))
      pointer(tam1ptr,tam1(jxa))
      pointer(tam2ptr,tam2(jxa))
      pointer(tam3ptr,tam3(jxa))
      pointer(tam4ptr,tam4(jxa))
      pointer(tam5ptr,tam5(jxa))
      pointer(tam6ptr,tam6(jxa))
      pointer(tam7ptr,tam7(jxa))
      pointer(tam8ptr,tam8(jxa))
      pointer(tam9ptr,tam9(jxa))
      pointer(tam10ptr,tam10(jxa))
      pointer(tam11ptr,tam11(jxa))
      pointer(tam12ptr,tam12(jxa))
      pointer(tam13ptr,tam13(jxa))
      pointer(tam14ptr,tam14(jxa))
      pointer(tam15ptr,tam15(jxa))
      pointer(tam16ptr,tam16(jxa))
      pointer(tam17ptr,tam17(jxa))
      pointer(tam18ptr,tam18(jxa))
      pointer(tam19ptr,tam19(jxa))
      pointer(tam20ptr,tam20(jxa))
      pointer(tam21ptr,tam21(jxa))
      pointer(tam22ptr,tam22(jxa))
      pointer(tam23ptr,tam23(jxa))
      pointer(tam24ptr,tam24(jxa))
      pointer(tam25ptr,tam25(jxa))
      pointer(tam26ptr,tam26(jxa))
      pointer(tam27ptr,tam27(jxa))
      pointer(tam28ptr,tam28(jxa))
      pointer(tam29ptr,tam29(jxa))
      pointer(tam30ptr,tam30(jxa))
      pointer(xptr,x(jxa))
      pointer(xmidpptr,xmidpt(jxa))
      pointer(xiptr,xi(jxa))
      pointer(xsqptr,xsq(jxa))
      pointer(x3iptr,x3i(jxa))
      pointer(x2iptr,x2i(jxa))
      pointer(xcuptr,xcu(jxa))
      pointer(xcentptr,xcenter(jxa))
      pointer(xcensptr,xcensq(jxa))
      pointer(uocptr,uoc(jxa))
      pointer(enerkptr,enerkev(jxa))
      pointer(gammaptr,gamma(jxa))
      pointer(gamsqptr,gamsqr(jxa))
      pointer(gamcuptr,gamcub(jxa))
      pointer(gammiptr,gammi(jxa))
      pointer(gamm2ptr,gamm2i(jxa))
      pointer(gamm1ptr,gamm1(jxa))
      pointer(tcsgmptr,tcsgm1(jxa))
      pointer(gamfacptr,gamefac(jxa))
      
      pointer(dxm5ptr,dxm5(jxp1a))
      pointer(exm5ptr,exm5(jxp1a))
c     dxp5ptr=dxm5ptr
      pointer(dxp5ptr,dxp5(0:jxa))
c     exp5ptr=exm5ptr
      pointer(exp5ptr,exp5(0:jxa))

      pointer(identptr,ident(iya))
      pointer(temc1ptr,temc1(iya))
      pointer(temc2ptr,temc2(iya))
      pointer(temc3ptr,temc3(iya))
      pointer(temc4ptr,temc4(iya))
      pointer(itmc1ptr,itemc1(iya))
      pointer(itmc2ptr,itemc2(iya))
      pointer(l_lowptr,l_lower(iya))
      pointer(lptptr,lpt(iya))
      pointer(munptr,mun(iya))

      pointer(fllptr,fll(jpxya))
      pointer(xparptr,xpar(jpxya))
      pointer(rheadptr,rheads(jpxya))
      pointer(dfvleptr,dfvlle(jpxya))
      pointer(dfvliptr,dfvlli(jpxya))

      pointer(xperpptr,xperp(ipxya))

      pointer(xlptr,xl(1:jfla))
      pointer(jmaxxptr,jmaxxl(1:jfla))
      pointer(xlmptr,xlm(0:jfla))
      pointer(dxlptr,dxl(0:jfla))
      pointer(flptr,fl(0:jfla))
      pointer(fl1ptr,fl1(0:jfla))
      pointer(fl2ptr,fl2(0:jfla))
      pointer(pparsptr,ppars(jxa,jfl))
      pointer(pprpsptr,pprps(jxa,jfl))
      pointer(faciptr,faci(jxa,jfl))
      pointer(pareaptr,pparea(jxa,jfl))
      pointer(wtfl0ptr,wtfl0(iya,jxa,lz))
      pointer(wtflmptr,wtflm(iya,jxa,lz))
      pointer(jflbnptr,jflbin(iya,jxa,lz))


      pointer(xmptr,xm(jxa,-5-mxa:4+mxa))

c     note: "ci" are arrays equivalenced to "di", i=a,..,f
c           for purposes of shifting the columns by one.
      pointer(daptr,da(iya,0:jxa))
      pointer(dbptr,db(iya,0:jxa))
      pointer(dcptr,dc(iya,0:jxa))
      pointer(dbbptr,dbb(iya,0:jxa))

      pointer(ddptr,dd(0:iya,jxa))
      pointer(deptr,de(0:iya,jxa))
      pointer(dfptr,df(0:iya,jxa))
      pointer(dffptr,dff(0:iya,jxa))

      pointer(cahptr,cah(iya,jxa))
      pointer(cthtaptr,cthta(iya,jxa))

      pointer(gonptr,gon(0:iyp1a,0:jxp1a))
      pointer(soptr,so(0:iyp1a,0:jxp1a))

c     NOTE:     We assume temp[1-6] are in contiguous storage,
c               so we can reference the whole six arrays through tem1.

c     NOTE: tem1,item1,xhead,tamt1,eg are arrays equivalenced to temp1
      pointer(temp1ptr,temp1(0:iyp1a,0:jxp1a))
c     note: tem2,item2,xtail,fg,fpn   are arrays equivalenced to temp2
      pointer(temp2ptr,temp2(0:iyp1a,0:jxp1a))
c     note: tem3,item3,ytail          are arrays equivalenced to temp3
      pointer(temp3ptr,temp3(0:iyp1a,0:jxp1a))
c     note: tem4,item4,yhead,tamt2    are arrays equivalenced to temp4
      pointer(temp4ptr,temp4(0:iyp1a,0:jxp1a))
c     note: tem5,item5                are arrays equivalenced to temp5
      pointer(temp5ptr,temp5(0:iyp1a,0:jxp1a))
c     note: tem6,item6                are arrays equivalenced to temp6
      pointer(temp6ptr,temp6(0:iyp1a,0:jxp1a))

c     xheadptr=temp1ptr
      pointer(xheadptr,xhead(jpxya,ipxya))
c     xtailptr=temp2ptr
      pointer(xtailptr,xtail(jpxya,ipxya))
c     ytailptr=temp3ptr
      pointer(ytailptr,ytail(jpxya,ipxya))
c     yheadptr=temp4ptr
      pointer(yheadptr,yhead(jpxya,ipxya))
c     (changed to separate allocation: bh960724) fpnptr=temp2ptr
c     bh970331: Changed back to fpnptr=temp2ptr
      pointer(fpnptr,fpn(jpxya,ipxya))
c     xlljiptr=ddptr
      pointer(xlljiptr,xllji(jpxya,ipxya))
c     xppjiptr=dfptr
      pointer(xppjiptr,xppji(jpxya,ipxya))

c     egptr=temp1ptr
      pointer(egptr,eg(iya,jxa))
c     fgptr=temp2ptr
      pointer(fgptr,fg(iya,jxa))

c     tem1ptr=temp1ptr
      pointer(tem1ptr,tem1(iyjx2a))
c     tem2ptr=temp2ptr
      pointer(tem2ptr,tem2(iyjx2a))
c     tem3ptr=temp3ptr
      pointer(tem3ptr,tem3(iyjx2a))
c     tem4ptr=temp4ptr
      pointer(tem4ptr,tem4(iyjx2a))
c     tem5ptr=temp5ptr
      pointer(tem5ptr,tem5(iyjx2a))
c     tem6ptr=temp6ptr
      pointer(tem6ptr,tem6(iyjx2a))

c     item1ptr=temp1ptr
      pointer(item1ptr,item1(iyjx2a))
c     item2ptr=temp2ptr
      pointer(item2ptr,item2(iyjx2a))
c     item3ptr=temp3ptr
      pointer(item3ptr,item3(iyjx2a))
c     item4ptr=temp4ptr
      pointer(item4ptr,item4(iyjx2a))
c     item5ptr=temp5ptr
      pointer(item5ptr,item5(iyjx2a))
c     item6ptr=temp6ptr
      pointer(item6ptr,item6(iyjx2a))

c     captr=daptr
      pointer(captr,ca(iya,jxa))
c     cbptr=dbptr
      pointer(cbptr,cb(iya,jxa))
c     ccptr=dcptr
      pointer(ccptr,cc(iya,jxa))
c     cdptr=ddptr
      pointer(cdptr,cd(iya,jxa))
c     ceptr=deptr
      pointer(ceptr,ce(iya,jxa))
c     cfptr=dfptr
      pointer(cfptr,cf(iya,jxa))

      pointer(currvptr,currv(jx,ngen,lrors))
      pointer(curr_ptr,currvs(jxa,ngena))
      pointer(pwrrfptr,pwrrf(jxa,ngena))
      pointer(talptr,tal(jxa,ngena))
      pointer(tblptr,tbl(jxa,ngena))
      pointer(tflptr,tfl(jxa,ngena))

      pointer(pwrr_ptr,pwrrfs(0:jxa,ngena))

      pointer(plegptr,pleg(0:mxa,iya))

      pointer(fetaptr,feta(jxa,0:mxa))
      pointer(fetbptr,fetb(jxa,0:mxa))

      pointer(wfluxptr,wflux(jxa,ngena,0:11))

      pointer(rhsptr,rhs(iyjxa))

      pointer(sovtptr,sovt(jxa,ngena,nsoa,lrzmax))

      pointer(sigsxptr,sigsxr(jxa,0:msxra,nena,2))

c.......................................................................
c*****arrays related to relativ=fully option
c.......................................................................
      pointer(gamnptr,gamman(jxa,-2:mxa+2))
      pointer(alpnptr,alphan(jxa,-2:mxa+2))

      pointer(asnhptr,asnha(jxa))

c     tamt1ptr=temp1ptr
      pointer(tamt1ptr,tamt1(2,jxa,-1:mxa+1,0:mxa+2))
c     tamt2ptr=temp4ptr
      pointer(tamt2ptr,tamt2(2,jxa,-1:mxa+1,0:mxa+2))
CPTR<<<END PTR-AINALLOC

c*****************************************************************
c     BEGIN arrays for analytic ion source (sou..) routines
c*****************************************************************

      common /diskx/
     1  bdre(lrza),bdrep(lrza),
     1  sorpwt(lrza),sorpwti(0:lrza),
     1  xlncurt(lrza)

      common /diskx/
     1  sorpw(0:ngena,lrza)

      common /diskx/
     1  cosm1(ngena,nsoa,lrza),cosm2(ngena,nsoa,lrza),
     1  sxllm1(ngena,nsoa,lrza),sxllm2(ngena,nsoa,lrza),
     1  sxppm1(ngena,nsoa,lrza),
     1  sxppm2(ngena,nsoa,lrza),
     1  xem1(ngena,nsoa,lrza),xem2(ngena,nsoa,lrza),
     1  zm1(ngena,nsoa,lrza),zm2(ngena,nsoa,lrza)

      common /diskx/
     1  sounor(ngena,nsoa,lza,lrza)

      common /scalar/
     1  isounor



c*****************************************************************
c     BEGIN arrays for rf package..(rf...,vlh[B,...,vlf...) routines
c*****************************************************************


      common/qlcoef/cqlb(iya,jxa,nmodsa),cqlc(iya,jxa,nmodsa)
      common/qlcoef/cqle(iya,jxa,nmodsa),cqlf(iya,jxa,nmodsa)

      common/qlcoef/ bqlm(iya,jxa)



c****************************************************************
c     BEGIN arrays for 3-d (td..) driver.
c****************************************************************


c..................................................................
c     scalars used in CQL3D
c..................................................................

      real*8 li
      common/sc3d/
     1  dttr,dtreff,currtza,currtpza,conserv,
     1  fom,fomp,fompla,fomtot,flxout,
     1  eden,edenlavg,etemp,ethtemp,edntmp,pden,pdntmp,psynct,
     1  li,
     1  ncount,iplt3d,
cBH070408(Not used excet in fr routines):    1  smooth,
     1  toteqd,cursign,totcurza,total,total0,totcurtt,curxjtot

c..................................................................
c     all arrays used only in CQL3D
c..................................................................

      real*8 jparb,jparbt,jparbp,mun
      common/ar3d/ rrz(0:lrza),
     1  tr(0:lrza),tr1(0:lrza),tr2(0:lrza),
     1  tr3(0:lrza),tr4(0:lrza),tr5(0:lrza),drp5(0:lrza),
     1  dpsi(0:lrza),dpsidrho(lrza),
     1  iytr(lrorsa),h_r(0:lrza),
     1  area(0:lrza),equilpsp(0:lrza),equilpsi(0:lrza),
     1  areamid(0:lrza),volmid(0:lrza),
     1  psivalm(lrza),rpconz(lrza),rmconz(lrza),rpmconz(0:lrza),
     1  bpolsqaz(0:lrza),aspin(lrza),trapfrac(lrza),
     1  currz(ngena,lrza),currtpz(lrza),currtz(lrza),
     1  currza(ngena),currtzi(0:lrza),currtpzi(0:lrza),
     1  currmtz(lrorsa),currmtpz(lrorsa),
     1  totcurzi(0:lrza),totcurz(lrza),
     1  fpsiz2(0:lrza),fpsiz(lrza),ffpsiz(lrza),
     1  jparb(lrza),jparbt(lrza),jparbp(lrza),
     1  prestp(lrza),prest(lrza),d2prest(lrza),d2fpsiz(lrza),
     1  d2ffpsiz(lrza),
     1  bmdplne(lrza),d2bmdpl(lrza),
     1  gkpwrz(ngena,lrza),
     1  rfpwrz(ngena,lrza)
      common/ar3d/
     1  dvol(lrza),darea(lrza),
     1  psyncz(lrza),pegyz(ngena,lrza),pplossz(ngena,lrza),
     1  wparzt(ngena),wperpzt(ngena),pegyt(ngena),pplosst(ngena),
     1  rfpwrt(ngena),
     1  gkpwrt(ngena),energyt(ntotala),
     1  rz(0:lrza),
     1  vfluxz(lrorsa),vol(0:lrza),
     1  onovrpz(lrza,2),
     1  constp(noncha,lrorsa),
     1  tplt3d(nplota),
     1  bscurma(2,2),bscurm(0:lrza,2,2),bscurmi(0:lrza,2,2)

      common/sc3d/
     1  sorpwtza

      common /csxr/ sxry(lrza,4),sang(lrza,4),spol(lrza,4),
     1  ibin(lrza,4),eflux(nena,nva),efluxt(nva),alphad(3),xs_(3),
     1  enk(nena),en_(nena),jval_(nena),
     1  inegsxr(nva),lensxr(nva)

      common /csigma/ mtab,msig,jxis,elmin,delegy,iind(1:jxa),
     1  imaxwln(2,4),igenrl(2,4),tamm1(0:mmsva),
     1  sigm(4,lrorsa),sigf(4,lrorsa),sigmt(4),sigft(4),
     1  sigmtt(noncha,4),sigftt(noncha,4),
     1  fuspwrv(4,lrorsa),fuspwrvt(4),fuspwrm(4,lrorsa),fuspwrmt(4)

c..............................................................
c     Set up pointers for sigma-v
c..............................................................

CPTR>>>REPLACE PTR-SIGMA
      common /dptr/ sgdumptr
      common /dptr/ csvptr,svtabptr

      pointer(sgdumptr,sgdum(1))
      pointer(csvptr,csv(jxis,0:mmsv,msig))
      pointer(svtabptr,svtab(mtab))
CPTR<<<END PTR-SIGMA

c..............................................................
c     Set up pointers to allocatable arrays for transport model.
c     Space allocated in subroutine tdtraloc
c..............................................................


CPTR>>>REPLACE PTR-TDTRALOC
      common /dptr/ tddumptr
      common /dptr/ frn_2ptr,frn_1ptr,frnptr,fvn_1ptr,fvnptr
      common /dptr/ d_rrptr,d_rptr,dlptr
      common /dptr/ f_lmptr,f_lpptr,f_upptr
      common /dptr/ f_vrptr
      common /dptr/ vpin_ptr,cynt_ptr,cosovptr,bovcoptr,vptb_ptr
      common /dptr/ advptr,dentaptr
      common /dptr/ eg_ptr,fg_ptr

      pointer(tddumptr,tddum(1))
      pointer(frn_2ptr,frn_2(0:iyp1a,0:jxp1a,ngena,0:lrz))
      pointer(frn_1ptr,frn_1(0:iyp1a,0:jxp1a,ngena,0:lrz))
      pointer(frnptr,frn(0:iyp1a,0:jxp1a,ngena,0:lrz))
      pointer(fvn_1ptr,fvn_1(0:iyp1a,0:jxp1a,ngena,0:lrz))
      pointer(fvnptr,fvn(0:iyp1a,0:jxp1a,ngena,0:lrz))
      pointer(dlptr,dl(0:iyp1a,0:jxp1a,ngena,0:lrz))
      pointer(d_rrptr,d_rr(0:iyp1a,0:jxp1a,ngena,0:lrz))
      pointer(d_rptr,d_r(0:iyp1a,0:jxp1a,ngena,0:lrz))

      pointer(f_lmptr,f_lm(jxa,ngena,lrz))
      pointer(f_lpptr,f_lp(jxa,ngena,lrz))
      pointer(f_upptr,f_up(jxa,ngena,lrz))

      pointer(f_vrptr,f_vtor(jxa,ngena,18,lrz))

      pointer(cynt_ptr,cynt2_(iya,lrz))

      pointer(vpin_ptr,vpint_(iya,lrzmax))

      pointer(vptb_ptr,vptb_(iya,0:lrzmax))

      pointer(cosovptr,cosovb(0:iya,0:lrz))
      pointer(bovcoptr,bovcos(0:iya,0:lrz))

      pointer(advptr,adv(ngena,lrzmax))

      pointer(dentaptr,dentarget(lrz))

      pointer(eg_ptr,eg_(jxa,2,lrz))
      pointer(fg_ptr,fg_(jxa,2,lrz))
CPTR<<<END PTR-TDTRALOC


c******************************************************************
c     BEGIN arrays for EQUILIBRIUM MODEL (eq..) (NON-CIRCULAR CROSS
c     SECTIONS).
c******************************************************************


      common /diskx/
     1  areacon(lrza),
     1  bmidplne(lrza),bpolsqa(lrza),
     1  eqdells(lrza),epsicon(lrza),erhocon(lrza),
     1  fpsi(lrza),flxavgd(lrza),
     1  psiovr(lrza),psiavg(2,lrza),onovpsir3(lrza),
     1  rpcon(lrza),rmcon(lrza),
     1  volcon(lrza),
     1  fppsi(lrza),pppsi(lrza)

      common/params/ nnz,nnr,nj12


      character*8
     1  eqorb,eqcall

      common
     1  bpolsqlm,
     1  eqorb,eqcall,
     1  imag,
     1  nmag,nrc,nzc,nfp,nnv,
     1  psimag,psilim,
     1  rmaxcon,rmincon,rhomax,
     1  zmag,zmaxcon,zmincon,zshift

      common
     1  ibd(4),
     1  eqpsi(nconteqa),eqvol(nconteqa),eqfopsi(nconteqa),
     1  q_(nconteqa),eqrho(nconteqa),eqarea(nconteqa),
     1  eqrpcon(nconteqa),eqrmcon(nconteqa),
     1  d2fpsiar(nnra),
     1  ez(nnza),dummyaz(nnza),
     1  fpsiar(nnra),ffpar(nnra),d2ffpar(nnra),qar(nnra),d2qar(nnra),
     1  prar(nnra),d2prar(nnra),ppar(nnra),d2ppar(nnra),psiar(nnra),
     1  er(nnra),dummyar(nnra),
     1  wkepsi(nrz3p1a),
     1  tlorb1(lfielda),tlorb2(lfielda)

      common
     1  epsi(nnra,nnza),epsirr(nnra,nnza),
     1  epsizz(nnra,nnza),epsirz(nnra,nnza),
     1  dummypsi(nnra,nnza),eqovrp(nconteqa,2)

      common/output/ lorbit_,ialign14,rmcon_,rpcon_,bthr_,btoru_,
     1  eqdells_,fpsi_,fppsi_,zmax_,btor0_,bthr0_,bmidplne_,
     1  solr_(lfielda),solz_(lfielda),es_(lfielda),
     1  eqbpol_(lfielda),bpsi_(lfielda),thtpol_(lfielda),
     1  eqdell_(lfielda)

c..................................................................
c     Allocatable arrays allocated in subroutine eqalloc
c..................................................................

      common/diskx/ onovrp(2,lrza)

CPTR>>>REPLACE PTR-EQALLOC
      common/dptr/ eqdumptr
      common/dptr/ eqdelptr,solrptr,eqbpoptr,solzptr

      pointer(eqdumptr,eqdum(1))

      pointer(eqdelptr,eqdell(lfielda,lrzmax))
      pointer(eqbpoptr,eqbpol(lfielda,lrzmax))
      pointer(solrptr,solr(lfielda,lrzmax))
      pointer(solzptr,solz(lfielda,lrzmax))
CPTR<<<END PTR-EQALLOC


c*********************************************************************
c     BEGIN arrays for LOWER HYBRID FAST WAVE and ECH Module.
c*********************************************************************


      common/params/
     1  jjx

      real*8 jbm1,jb0,jbp1

CPTR>>>REPLACE PTR-JB
      pointer(jbm1ptr,jbm1(nbssltbl,mrfn))
      pointer(jb0ptr,jb0(nbssltbl,mrfn))
      pointer(jbp1ptr,jbp1(nbssltbl,mrfn))
CPTR<<<END PTR-JB

      common
     1  argmax,
     1  nurf,
CPTR>>>DELETE PTR-JB
     1  ipad5,jbm1ptr,jb0ptr,jbp1ptr,
CPTR<<<END PTR-JB
     1  lenj0,
     1  mrf,mrfn,
     1  irftype,
     1  powray,
     1  vnorm2,vnorm3,vnorm4

      complex*16 cosz1,sinz1,sinz2
      common
     1  besl(nharm2a),
     1  bsslstp(nmodsa),
     1  ncontrib(lrza),
     1  powrf(lrza,nmodsa),powrfc(lrza,nmodsa),
     1  powrfl(lrza,nmodsa),
     1  powurf(0:nmodsa),powurfc(0:nmodsa),powurfl(0:nmodsa),
     1  powurfi(0:lrza,0:nmodsa),powrft(lrza)

      common
     1  freqcy(nmodsa),omega(nmodsa),nharm(nmodsa),nray(nmodsa),
     1  bsign1(nmodsa),krfn(nmodsa),irfn(nmodsa)


c..................................................................
c     Allocatable arrays allocated in subroutine urfalloc
c..................................................................

CPTR>>>DELETE PTR-URFALLOC
      common /dptr/ urfbptr,urfcptr,urfeptr,urffptr,cosmzptr,g_ptr
      common /dptr/ urfduptr
      common /dptr/ alfiptr,alfaptr
      common /dptr/ alfagptr,argmtptr
      common /dptr/ ilim1ptr,ilim2ptr,ilimdptr,ili_ptr
      common /dptr/ ifct1ptr,ifct2ptr,if1ptr,if2ptr
      common /dptr/ sxptr
      common /dptr/ xmdxptr
      common /dptr/ cosz1ptr,sinz1ptr,sinz2ptr
      common /dptr/ thtf1ptr,thtf2ptr
      common /dptr/ ilin1ptr,ilin2ptr
      common /dptr/ urftmpptr,urfpwptr,urfpcptr,urfplptr
      common /dptr/ jminptr,jmaxptr
      common /dptr/ llocptr,llrayptr
      common /dptr/ psiloptr
      common /dptr/ scal_ptr
      common /dptr/ cwexdptr
      common /dptr/ cweydptr,cwezdptr
      common /dptr/ delpwptr
      common /dptr/ fluxnptr
      common /dptr/ seikoptr,spsiptr
      common /dptr/ sdpwrptr,sbtotptr
      common /dptr/ seneptr
      common /dptr/ salp1ptr,salp2ptr
      common /dptr/ wsptr,wrptr,wzptr
      common /dptr/ wnparptr,wdnpaptr
      common /dptr/ wnperptr,wphiptr
      common /dptr/ ilowptr,iupptr
      common /dptr/ nrayeptr,jslosptr,nurefptr,keiksptr
      common /dptr/ jpesptr,jpisptr,istatptr
      common /dptr/ iprmtptr,jhlfsptr
      common /dptr/ sxxrtptr,skpsiptr,skthptr,skphiptr
      common /dptr/ lrayeptr
      common /dptr/ delp0ptr
      common /dptr/ nray0ptr
CPTR<<<END PTR-URFALLOC

      complex*16 cwexde,cweyde,cwezde

CPTR>>>REPLACE PTR-URFALLOC
      pointer(urfduptr,urfdum(1))

      pointer(urfbptr,urfb(iy,jx,lrz,mrfn))
      pointer(urfcptr,urfc(iy,jx,lrz,mrfn))
      pointer(urfeptr,urfe(iy,jx,lrz,mrfn))
      pointer(urffptr,urff(iy,jx,lrz,mrfn))

      pointer(cosmzptr,cosmz(iya,lza,lrzmax))

      pointer(g_ptr,g_(0:iyp1a,0:jxp1a,ngena,lrors))

      pointer(alfagptr,alfag(jxa))
      pointer(argmtptr,argmnt(jxa))
      pointer(ilim1ptr,ilim1d(jxa))
      pointer(ilim2ptr,ilim2d(jxa))
      pointer(ilimdptr,ilim1dd(jxa))
      pointer(ili_ptr,ilim2dd(jxa))
      pointer(sxptr,sx(jxa))
      pointer(xmdxptr,xmdx(jxa))

c     complex*16 arrays
      pointer(cosz1ptr,cosz1(iya))
      pointer(sinz1ptr,sinz1(iya))
      pointer(sinz2ptr,sinz2(iya))

      pointer(thtf1ptr,thtf1(iya))
      pointer(thtf2ptr,thtf2(iya))
      pointer(alfiptr,alfi(iya))
      pointer(alfaptr,alfa(iya))

      pointer(ilin1ptr,ilim1(jjxa))
      pointer(ilin2ptr,ilim2(jjxa))
      pointer(ifct1ptr,ifct1(jjxa))
      pointer(ifct2ptr,ifct2(jjxa))

      pointer(urftmpptr,urftmp(nrayelta*nraya*2))
      pointer(urfpwptr,urfpwr(nrayelta,nraya,mrfn))
      pointer(urfpcptr,urfpwrc(nrayelta,nraya,mrfn))
      pointer(urfplptr,urfpwrl(nrayelta,nraya,mrfn))
      pointer(jminptr,jminray(nrayelta,nraya,mrfn))
      pointer(jmaxptr,jmaxray(nrayelta,nraya,mrfn))
      pointer(llocptr,lloc(nrayelta,nraya,mrfn))
      pointer(llrayptr,llray(nrayelta,nraya,mrfn))
      pointer(psiloptr,psiloc(nrayelta,nraya,mrfn))
      pointer(scal_ptr,scalurf(nrayelta,nraya,mrfn))

c     complex*16 arrays
      pointer(cwexdptr,cwexde(nrayelta,nraya,mrfn))
      pointer(cweydptr,cweyde(nrayelta,nraya,mrfn))
      pointer(cwezdptr,cwezde(nrayelta,nraya,mrfn))

      pointer(delpwptr,delpwr(nrayelta,nraya,mrfn))
      pointer(fluxnptr,fluxn(nrayelta,nraya,mrfn))
      pointer(seikoptr,seikon(nrayelta,nraya,mrfn))
      pointer(spsiptr,spsi(nrayelta,nraya,mrfn))
      pointer(sdpwrptr,sdpwr(nrayelta,nraya,mrfn))
      pointer(sbtotptr,sbtot(nrayelta,nraya,mrfn))
      pointer(seneptr,sene(nrayelta,nraya,mrfn))
      pointer(salp1ptr,salphac(nrayelta,nraya,mrfn))
      pointer(salp2ptr,salphal(nrayelta,nraya,mrfn))
      pointer(wsptr,ws(nrayelta,nraya,mrfn))
      pointer(wrptr,wr(nrayelta,nraya,mrfn))
      pointer(wzptr,wz(nrayelta,nraya,mrfn))
      pointer(wnparptr,wnpar(nrayelta,nraya,mrfn))
      pointer(wdnpaptr,wdnpar(nrayelta,nraya,mrfn))
      pointer(wnperptr,wnper(nrayelta,nraya,mrfn))
      pointer(wphiptr,wphi(nrayelta,nraya,mrfn))

      pointer(ilowptr,ilowp(ipacka,mrfn))
      pointer(iupptr,iupp(ipacka,mrfn))
      pointer(if1ptr,ifct1_(ipack16a,mrfn))
      pointer(if2ptr,ifct2_(ipack16a,mrfn))

      pointer(nrayeptr,nrayelt(nraya,mrfn))
      pointer(jslosptr,jslofas(nraya,mrfn))
      pointer(nurefptr,nurefls(nraya,mrfn))
      pointer(keiksptr,keiks(nraya,mrfn))
      pointer(jpesptr,jpes(nraya,mrfn))
      pointer(jpisptr,jpis(nraya,mrfn))
      pointer(istatptr,istarts(nraya,mrfn))
      pointer(iprmtptr,iprmt5(nraya,mrfn))
      pointer(jhlfsptr,jhlfs(nraya,mrfn))
      pointer(sxxrtptr,sxxrt(nraya,mrfn))
      pointer(skpsiptr,skpsi(nraya,mrfn))
      pointer(skthptr,skth(nraya,mrfn))
      pointer(skphiptr,skphi(nraya,mrfn))
      pointer(lrayeptr,lrayelt(nraya,mrfn))
      pointer(delp0ptr,delpwr0(nraya,mrfn))
      pointer(nray0ptr,nrayelt0(nraya,mrfn))
CPTR<<<END PTR-URFALLOC

c..................................................................
c     Allocatable arrays allocated in subroutine rdc_multi,
c     used after subroutine execution.
c     Here, we introduce f90 pointers, as they are easier
c     to allocate.

      pointer rdcb
      dimension rdcb(:,:,:)
      common /dptr95/ rdcb
      pointer rdcc
      dimension rdcc(:,:,:)
      common /dptr95/ rdcc
      pointer rdce
      dimension rdce(:,:,:)
      common /dptr95/ rdce
      pointer rdcf
      dimension rdcf(:,:,:)
      common /dptr95/ rdcf


c..................................................................

      common anecc(nrada),tekev(nrada),tikev(nint1a,nrada),
     1  anicc(nint1a,nrada),
     1  amass(nint1a),achrg(nint1a),elecf(nrada),
     1  rho_(nrada),names(10),psiar_(nrada),
     1  nspc,ialign9,fpsiar_(nrada),pary(nrada),ppary(nrada),
     1  gpary(nrada),ztop_,zbot_,rleft,rright,nx,nz,
     1  npsitm

c-----------------------------------------------------------------------
c     BEGIN variables for WP... modules for CQLP case
c-----------------------------------------------------------------------

      character*8  
     1  analegco

      common /wpscal/
     1  analegco,
     1  iymax,
     1  nsleft,nsrigt,numclas,numindx

      common /wpvec/
     1  cofdfds(0:lsa1,2,ntrmdera),
     1  enrgypa(ntotala,0:lsa1),
     1  vthpar(ntotala,0:lsa1),
     1  ilpm1ef(0:iya+1,0:lsa1,-1:+1),
     1  lsbtopr(0:lsa1),lsprtob(0:lsa1),lpm1eff(0:lsa1,-1:+1),
     1  sz(0:lsa1),dsz(0:lsa1),dszm5(0:lsa1),dszp5(0:lsa1),
     1  eszm5(0:lsa1),eszp5(0:lsa1),
     1  psis(0:lsa1),psisp(0:lsa1),psipols(0:lsa1),
     1  solrs(0:lsa1),solzs(0:lsa1),
     1  elparol(0:lsa1),elparnw(0:lsa1),
     1  l_upper(iya),flux1(0:lsa1),flux2(0:lsa1)

c.......................................................................
c     Arrays allocated in subroutine wpalloc for CQLP
c.......................................................................

CPTR>>>REPLACE PTR-WPALLOC
      common /dptr/ wpdumptr
      common /dptr/ fnhalfptr,fnp0ptr,fnp1ptr,dlsptr
      common /dptr/ fhptr
      common /dptr/ fedgeptr
      common /dptr/ rhsparptr,bndmatsptr
      common /dptr/ wcqlbptr,wcqlcptr,wcqleptr,wcqlfptr

      pointer(wpdumptr,wpdum(1))

      pointer(fnhalfptr,fnhalf(0:iyp1a,0:jxp1a,ngena,0:ls+1))
      pointer(fnp0ptr,fnp0(0:iyp1a,0:jxp1a,ngena,0:ls+1))
      pointer(fnp1ptr,fnp1(0:iyp1a,0:jxp1a,ngena,0:ls+1))
      pointer(dlsptr,dls(0:iyp1a,0:jxp1a,ngena,0:ls+1))
      pointer(fhptr,fh(0:iyp1a,0:jxp1a,ngena,0:ls+1))

      pointer(fedgeptr,fedge(iya,jxa,ngena,4))

c     second dimension: max of iya and jxa
      pointer(rhsparptr,rhspar(0:lsa+nbanda+1,jxa,2))

c     first dimension: max of iya and jxa
      pointer(bndmatsptr,bndmats(jxa,0:lsa+1,nbanda,2))

      pointer(wcqlbptr,wcqlb(iya,jxa,nmodsa,ls))
      pointer(wcqlcptr,wcqlc(iya,jxa,nmodsa,ls))
      pointer(wcqleptr,wcqle(iya,jxa,nmodsa,ls))
      pointer(wcqlfptr,wcqlf(iya,jxa,nmodsa,ls))
CPTR<<<END PTR-WPALLOC

