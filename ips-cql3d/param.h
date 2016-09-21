c     param.h
c**********************************************************************
c**********************************************************************

c     PARAMETERS CHOSEN BY USER FOLLOW

c..................................................................
c     version is simply a mnemonic printed out to indicate
c     the version of the source being run.  Further below, the non-
c     user set parameter precursr indicates from which author-released
c     version this source is descended.
c
c_cray machinea is =1 for 64-bit integers (e.g., Cray), and
c_pc               =2 for 32-bit integers (e.g., PCs).
cBH081218:  Present usage, machinea=2 works with 32- and 64-bit machines
c
c     ngena is the maximum number of general (time advanced) species.
c
c     nmaxa is the maximum number of background species which are
c     not time advanced and which are assumed to be Maxwellian at a
c     fixed temperature and density. These species exist only
c     to contribute to the collision integral.
c
c     mx is the order of Legendre terms used in Rosenbluth
c     expansions of the collision integral.
c     NOTE:   If relativ="enabled" mx => 2,
c             BUT, 2*jx*(mx+3)*(mx+3) <= (iyp1+1)*(jxp1+1) to
c                  avoid overwrite by tamt1 and tamt2.
c
c     lza is the number of poloidal mesh points for purposes of
c     computing bounce averages (number of "z" mesh points).
c     (if cqlpmod="enabled", lza=lsa=number of mesh points
c     along the magnetic field line)
c
c     jx is the number of speed (momentum) mesh points.
c
c     iy is the number of pitch angle mesh points.
c
c     noncha is the number of elements in arrays carrying time plot
c     information.
c
c     nbctimea is max number of points in arrays supplying time
c     dependent profile information.
c
c     ndtr1a is maximum number of time step intervals dtr1().
c
c     nplota is maximum number of plot times for 2d and 3d plots.
c
c     nefitera is the maximum number of iterations permitted for
c     the electric field per time step (to obtain a target current).
c
c..................................................................

c     PARAMETERS CHOSEN BY USER FOLLOW

      character version*64
      parameter(version="cql3d_cswim_131030_mpi")
      character precursr*64
      parameter(precursr="cql3d_cswim_130928_mpi")
      parameter(machinea=2)
cBH081218:  Present usage, machinea=2 works with 32- and 64-bit machines
      parameter(ngena=4)
      parameter(nmaxa=8)
ccc      parameter(mxa=3)     YuP-101216: Not used anymore
ccc      parameter(jxa=300)   YuP-101216: Not used anymore
ccc      parameter(iya=160)   YuP-101216: Not used anymore
ccc      parameter(noncha=2000) YuP-101221: Not used anymore
      parameter(nbctimea=101)
      parameter(ndtr1a=10)
      parameter(nplota=10)
      parameter(nefitera=10)

c*******************************************************************
c     BEGIN parameters for SOURCE routines
c*******************************************************************

c..................................................................
c     nsoa is the maximum number of sources per species allowed.
c..................................................................

      parameter (nsoa=3)
c..................................................................
c     jfla is maximum length of grid for the reduced parallel
c     distribution function used for the Besedin-Pankratov 
c     knock-on source. Use odd number.
c..................................................................
ccc      parameter (jfla=jxa)    YuP-101216: Not used anymore


c*******************************************************************
c     Parameter for knock-on source soureko routine
c*******************************************************************

c..................................................................
c     
c     i0param is a table length in pitch angle, for translation
c     of sin(theta) a given poloidal position into the corresponding 
c     equatorial pitch angle.
c
c..................................................................

      parameter (i0param=1001)


c*******************************************************************
c     BEGIN parameters for 3-D (td...) routines
c*******************************************************************

c..................................................................
c     lrza is the number of radial mesh points
c
c     lsa is the number of parallel (spatial) mesh points
c     lsa should be .ge.lrza [BH081106:  perhaps this the usual
c       case before big memory and bigger lrza cases.   Its
c       contravention is causing an out-of-bounds ref for
c       cqlpmod.ne.enabled.   Could check code out further
c       to see how essential it is.]
c
c     Note: Limitations on relative sizes of lrza, lsa, and lza!
c
c     nbanda is the maximum bandwidth of the parallel transport equation
c
c     msxr is the order of the highest polynomial used in the Legendre
c     expansion for computing SXR spectra (It is set to mx in aindflt.f).
c
c     nena is the maximum number of energies in the SXR spectra calc.
c
c     nva is the maximum number of SXR sightlines.
c
c     ntrmdera is the maximum of node for computing df/ds
c
c     njenea is the number of nodes for spline profiles
c     (ryain, tein, tiin, enein, elecin)
c..................................................................

      parameter(lrza=128)
      parameter(lsa=max(128,lrza), lsa1=lsa+1)
c     lza should be .ge. lsa. if cqlpmod, lz is set equal to ls
      parameter(lza=lsa)
c     lrorsa should be equal to max of lrza and lsa
      parameter(lrorsa=lsa)
      
      parameter(nbanda=5)
      
      parameter(nena=60,nva=32)
      parameter(ntrmdera=4)
      parameter(njenea=256)

c*******************************************************************
c     BEGIN parameters for the NON-CIRCULAR CROSS SECTION (eq..)
c     routines.
c*******************************************************************

c..................................................................
c     nnra (nnza) is the number of R (Z) meshpoints used in
c     MHD equilibrium calculation
c
c     lfielda is the maximum number of poloidal points allowed
c     in the orbit calculation. (250 is a good number)
c..................................................................

c      parameter(nnra=201)
c      parameter(nnza=201)
      parameter(nnra=257)
      parameter(nnza=257)
      parameter(lfielda=250)
c-YuP-> Grid size for storing the values 
c     of equilibrium field (Beqr,Beqz,Beqphi, etc.)
c     on req(1:nreqa),zeq(1:nzeqa) grid
c     for finite-orbit-width (FOW) tracing.
      parameter(nreqa=128, nzeqa=196)
c-YuP-> Could be set to (nnra,nnza) ?


c*******************************************************************
c     BEGIN parameters for the NFREYA (beam deposition) (fr..) code.
c*******************************************************************

c..................................................................
c     Parameters for the "fr" routines (NFREYA code) as it exists in
c     the transport code ONETWO. Parameters and common blocks
c     appearing will be defined as in ONETWO.
c     Where possible the coding has been lifted from ONETWO.
c
c..................................................................

ccc      parameter(maxp=1500)  YuP-101211: maxp is not used anymore.

c*******************************************************************
c     BEGIN parameters for LH, FW,  and ECH (urf..) routines
c*******************************************************************

c..................................................................
c     nraya is the maximum allowable number of LH rays
c
c     nrayelta is the maximum allowable number of ray elements
c     per ray. A ray element is a small subsection of a ray
c     characterized by a position, a length, to which the
c     ray tracing code assigns a number of wave characteristics.
c
cBH060314  nmodsa is the maximum number of either number of wave modes,
cBH060314  or number of harmonics for a single wave mode (presently
cBH060314  there is some hardwiring to values .le.3).(and some hardwiring
cBH060314  to values .ge. 3. E.G., nrfstep1,powrf,powurf,powurfi).
cBH060314  Need to check this out before changing nmodsa.ne.3! 
cBH060314  (BobH, 990702).  This has been upgraded as follows.

c     nmodsa is maximum of the sum over{ wave types (mrf) times
c     harmonics for each wave type}.  (BobH, 060314).
c
c     nharma is maximum harmonic number for cyclotron interactions.
c..................................................................

c      parameter (nraya=1) 
c      parameter (nrayelta=1) 
c     YuP 101122: nraya and nrayelta are not used anymore.
c     YuP 101122: Instead, nrayn and nrayelts are determined 
c     YuP 101122: in urfsetup by reading rays\' data files.
      
      parameter (nmodsa=20)  !Suggest not using .le.3, unless check
                             !cqlinput that some vestigial inputs
                             !are not set for index larger than nmodsa.
ccc      parameter (nharma=1,nharm2a=nharma+2)
c     YuP 101208: nharma is not used anymore.
c..................................................................
c     NPA related:
      parameter(npaproca=5)
c..................................................................


c*******************************************************************
c     Parameters for sigma-v routines
c*******************************************************************

c..................................................................
c     mtaba is a table length for passage of sigma values.
c
c     mmsv is the order of the highest polynomial used in the Legendre
c     expansion for computing sigma-v (set to mx in aindflt.f)
c
c..................................................................

      parameter (mtaba=1000)
ccc      parameter (mmsva=mxa)    YuP-101216: Not used anymore


c     END OF PARAMETER SET WHICH IS NORMALLY USER CHOSEN FOLLOW
c***********************************************************************



c..................................................................
c     Parameters defined just below should not be altered
c     (with the possible exception of negyrga (used in plots)).
c..................................................................

      parameter(negyrga=3)
      parameter(mbeta=10)
ccc      parameter(mxp1a=mxa+1)  YuP-101216: Not used anymore
ccc      parameter(ngenpa=ngena+1)  Not used anymore
ccc      parameter(iyjxa=iya*jxa) YuP-101216: Not used anymore
ccc      parameter(iyjx2a=(iya+2)*(jxa+2)) YuP-101216: Not used anymore
ccc      parameter(iyjxnga=iyjxa*ngena) YuP-101216: Not used anymore
ccc      parameter(jxp1a=jxa+1) YuP-101216: Not used anymore
ccc      parameter(iyp1a=iya+1) YuP-101216: Not used anymore
      parameter(ntotala=ngena+nmaxa)
      parameter(ift07a=01,ift21a=01)

ccc      parameter(ipxya=iya)    YuP-101216: Not used anymore
ccc      parameter(jpxya=jxa+1)  YuP-101216: Not used anymore
ccc      parameter(iyjxua=iya*(jxa+1)) YuP-101216: Not used anymore
ccc      parameter(iyujxa=(iya+1)*jxa) YuP-101216: Not used anymore
ccc      parameter(miyjxa=6*iyjxa)     YuP-101216: Not used ?
      parameter(incza=301,inczpa=incza+1)

      integer tlfld1a
      parameter(tlfld1a=3*lfielda+1)
c     parameter(nconteqa=nnra)
      parameter(nnrab2=(nnra+1)/2)
      parameter(nconteqa=nnrab2)
      parameter(nrz3p1a=3*(nnza+nnra)+1)


      parameter(ki=nnra,kix2=2,kj=nnza,kjx2=2)
      parameter(kikj=ki*kj,kwork=3*(kj-1)+kj)

c..................................................................
c     kb is the maximum number of neutral injectors
c     ke (=3) is the number of beam energy states
c     kf is the maximum number of flux zones (or volumes) used
c     by FREYA, and HEX** routines (reaction rate routines).
c     nap is the number of source aperatures.
c     kprim is the maximum allowable number of primary species
c     kimp is the maximum number of impurities.
c     kion is the maximum number of ionic species
c..................................................................

      parameter(kprim=ntotala)
      parameter(kimp=ntotala)
      parameter(kion=ntotala)
      parameter(kkq=kion+3,kbctim=11,kb=8,ke=3,kf=nconteqa+3,
     1  kz=kf-1,kzm1=kz-1)

c..................................................................
c     ONETWO DIVERGENCE and source of confusion. Cray32 uses kj
c     as the number of R mesh points, cray30 and the rest of the
c     ONETWO code uses kj as he number of radial mesh points.
c     We will retain the meaning as in cray32. In some of the
c     common blocks below which are lifted from ONETWO kj
c     has been changed to k_. These arrays are not used in CQL3d,
c     but are retained to maintain continuity with ONETWO.
c..................................................................

      parameter(k_=3,k_m1=k_-1)
c     WARNING: DO NOT ALTER nap UNLESS IT IS ALTERED ALSO IN SUB ROTATE
      parameter(nap=10)
ccc      parameter(kjp=maxp+1) YuP-101211: Not used?

c     ibytes is bytes per integer word, for 64 or 32-bit integers.
c     (1 byte contains 8 bits).  
c     It is used for packing data for urf subroutines.
c     jjxa is 1st multiple of 8 greater than jxa.
      parameter(ibytes=8/machinea)      
ccc      parameter(jjxa=((jxa-1)/ibytes)*ibytes+ibytes)
ccc      parameter(ipacka=jjxa/ibytes*nrayelta*nraya) !moved to urfalloc

c     Set up new dimensioning for ifct1_,ifct2_ (from previous
c     1 byte per word to 2 bytes per word, improving accuracy.
c     BH, 050812).
c     ibytes16 is number of 16-bit words per integer word.
c     ipack16a is number of integer words required to store 1 set
c     of ray data (*jjxa) in the ifct1_,ifct2_ 16-bit-word arrays.
      parameter(ibytes16=ibytes/2)
ccc      parameter(ipack16a=jjxa/ibytes16*nrayelta*nraya)!->to urfalloc

cBH070118      parameter (nrada=129,ninta=8,nint1a=ninta+1)
      parameter (nrada=nnra,ninta=8,nint1a=ninta+1)

c.......................................................................
c     maximum number of options in tdoutput (for nlotp1,.. arrays)
c.......................................................................

      parameter(noutpta=10)
