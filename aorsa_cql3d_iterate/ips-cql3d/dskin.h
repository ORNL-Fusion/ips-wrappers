
c..................................................................
c
c   Include file for dskin.f, 
c     which reads diskf distribution function file from CQL3D.
c
c   parameters iya,jxa,lrza,ngena need to be set 
c     in accord with iy,jx,lrz,ngen in diskf file.
c
c   Note: You might need to know that storage is slightly 
c         different here than in CQL3D.
c
c..................................................................



c      parameter(iya=80)
      parameter(iya=120)
c      parameter(jxa=80)
      parameter(jxa=120)
c      parameter(lrza=15)
      parameter(lrza=20)
      parameter(ngena=1)


      real*8 x(jxa),y(iya,lrza),rovera(lrza),elecfld(lrza),
     +       bthr(lrza),btoru(lrza),bthr0(lrza),btor0(lrza),
     +       reden(lrza,ngena),temp(lrza,ngena)
      real*8 radmin,vnorm,vmaxdvt,eovedd
      real*8 bnumb(ngena),fmass(ngena)
      integer*4 iy(lrza),itl(lrza),itu(lrza)
      real*8 f(iya,jxa,lrza,ngena)
 
      COMMON /dskincomm/x,y,rovera,elecfld,bthr,btoru,bthr0,btor0
      COMMON /dskincomm/reden,temp,radmin,vnorm,vmaxdvt,eovedd
      COMMON /dskincomm/bnumb,fmass,itl,itu
      COMMON /dskincomm/f

c      IDL related test of pointers:
c      common/dptr/fptr,roveraptr,radminptr
c      pointer(fptr,f_idl(1:iya,1:jxa,lrza,ngena))
c      pointer(roveraptr,rovera_idl(1:lrza))
c      pointer(radminptr,radmin_idl)
