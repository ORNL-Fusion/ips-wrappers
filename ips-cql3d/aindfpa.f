c
c
      subroutine aindfpa
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     Set defaults for input of variables in first namelist setup,
c     now (as of 2007) also explicitly named setup0.
c     As of 2017, also can denote the first setup as namelist fsetup.
c     (fsetup worked around a problem with the pathscale compiler.)
c     Warning: In this subroutine, should set only variables read 
c     in first namelist setup/setup0/fsetup.
c     This also shows which variables belong to the 1st namelist setup.
c..................................................................

      include 'param.h'
      include 'name_decl.h'
c.......................................................................

      ibox(1)="unset"
      ibox(2)="unset"
      ibox(3)="unset"
      iuser="unset"
      ioutput(1)=0 !YuP[2020-10-21] was 6, but was not used anywhere
      ! After 2020, the usage of ioutput(1):
      ! ioutput(1)=0  means the printout to screen is reduced to a minimum,
      !               leaving namelist printout, warning messages,
      !               and physics-related values used or computed by code.
      !               (default value ioutput(1)=0)
      ! ioutput(1)=1  The above, and more printout - 
      !               values of arrays for diagnostic purpose.
      ! ioutput(1)=2  The above, and even more diagnostic printout.
      ioutput(2)=0 ! No usage at present [2020]
      
      mnemonic="mnemonic"

c     Avoid some special calls, if special_calls=disabled in 
c     &setup0 namelist.
c     [System calls not supported on some machines.]
c     [Could use this for other system dependent branching, in future....]

      special_calls="enabled"

c     main parameters (used to allocate memory and determine model)
      lrz=0
      lrzmax=0
      lrzdiff="disabled"
      cqlpmod="disabled"
      ls=0
      lsmax=0
      lsdiff="disabled"
      do 100 ll=0,lrorsa
        lrindx(ll)=ll
        lsindx(ll)=ll
 100  continue

      noplots="disabled"
      nmlstout="trnscrib"
      lnwidth=3  ! Plot line width in units of 0.005 inches 

c      nlwritf=.false.
c      nlrestrt=.false.

      nlwritf="disabled"
      nlrestrt="disabled"

      return
      end
c
c
c==================================================================
      subroutine ainadjnl(kopt)
      implicit integer (i-n), real*8 (a-h,o-z)
CMPIINSERT_INCLUDE
      Save inlmod

c..................................................................
c     kopt=0: Test if two &setup namelist sections in the cqlinput
c             If so, make copy of cqlinput to cqlinput_tmp
c             Modify cqlinput to &setup0 and &setup sections
c     kopt=1: Restore original cqlinput (if modification was made)
c
cBH070414:    This adjustment of namelist structure is to facilitate
c             writes and reads of namelist for the SWIM Integrated
c             Plasma Simulator (IPS), and to maintain backwards
c             compatibility.
c             Previously, there were two namelist sections both
c             designated &setup.  For SWIM, it has been necessary
c             to distinguish these (we use &setup0 and &setup),
c             for purposes of writing out the complete input namelist.
c..................................................................

      character*128  command
      character*8 line
      character*9 line0
      parameter (long_enough=100000)  !To accomodate namelist sections 
                                       !written on one line
      character(len=long_enough) :: line1   !Automatic array
      character(len=long_enough+1) :: line2   !Automatic array

c.......................................................................

CMPIINSERT_IF_RANK_NE_0_RETURN

      if (kopt.eq.0) then

         inlmod=0
         open(unit=4,file="cqlinput",delim='apostrophe',status="old") 
 1       read(unit=4,fmt=1003,end=2) line
         if (line.eq." &setup" .or. line.eq." &SETUP"
     +       .or. line.eq."&setup" .or. line.eq."&SETUP") then
            inlmod=inlmod+1
            if (inlmod.eq.2) go to 2
         endif
         go to 1
 2       continue
         WRITE(*,*)'inlmod = ',inlmod
         close(4)
         
 1003    format(a8)
         
c.......................................................................
c     Make copy if inlmod=2, and adjust nl structure
c.......................................................................
         
         if (inlmod.le.1) then
            goto 999            !Implies new namelist &setup0/&setup
                                !regimen, i.e., only one &setup.
            
         else
            
c     Copy cqlinput to cqlinput_tmp
            open(unit=4,file="cqlinput",delim='apostrophe',status="old")
            open(unit=5,file="cqlinput_tmp",delim='apostrophe',
     +           status="replace")
 3          read(unit=4,fmt='(a)',end=4) line1
            len_line1=len_trim(line1)
            if (len_line1.ge.(long_enough-1)) then
               WRITE(*,*)'len_line1,long_enough',len_line1,long_enough
               STOP 'Adjust long_enough'
            endif
            write(unit=5, fmt='(a)') trim(line1)
            go to 3
 4          continue
            close(4)
            close(5)
            
c     Reopen cqlinput as new file, and put adjusted namelist into it.
c     There are 2 &setup namelist sections.
            open(unit=4,file="cqlinput",delim='apostrophe',
     +           status="replace")
            open(unit=5,file="cqlinput_tmp",delim='apostrophe',
     +           status="old")
            ifirst=0
 5          read(unit=5,fmt='(a)',end=6) line1
            if ((line1(1:8).eq." &setup" .or. line1(1:8).eq." &SETUP"
     .          .or.line1(1:7).eq."&setup" .or. line1(1:7).eq."&SETUP")
     .           .and. (ifirst.eq.0)) then
               ifirst=1
cBH080118               line2(1:9)=" &setup0"
               line2(1:8)=" &setup0"
cBH080118                line2(10:len(line1)+1)=line1(9:len(line1))
               line2(9:len(line1)+1)=line1(9:len(line1))
            else
               line2=line1(1:len(line1))
            endif
            write(unit=4, fmt='(a)') trim(line2)
            go to 5
 6          continue
            
         endif
         
         close(4)
         close(5)
         go to 999
         
      elseif (kopt.eq.1) then
         
c.......................................................................
c     If cqlinput changed, copy cqlinput_tmp back to cqlinput
c.......................................................................
         if (inlmod.eq.2) then
            open(unit=5,file="cqlinput",delim='apostrophe',
     +           status="replace")
            open(unit=4,file="cqlinput_tmp",delim='apostrophe',
     +           status="old")
 7          read(unit=4,fmt='(a)',end=8) line1
            if (len_trim(line1).ge.(long_enough-1)) 
     .           STOP 'Adjust long_enough'
            write(unit=5, fmt='(a)') trim(line1)
            go to 7
 8          continue
            close(4)
            close(5)
            
         endif
         
      endif                     !on kopt
      
 999  return
      end

c
c
c==================================================================
      subroutine ainadjnl_fsetup_setup0(kopt)
      !YuP[2017] 
      implicit integer (i-n), real*8 (a-h,o-z)
CMPIINSERT_INCLUDE
      Save inlmod

c..................................................................
c     kopt=0: Test if &fsetup namelist sections in the cqlinput
c             is present (instead of setup0, as in svn versions of CQL3D)
c             If so, make copy of cqlinput to cqlinput_tmp
c             Modify/rename $fsetup to &setup0 
c     kopt=1: Restore original cqlinput (if modification was made)
c     
c     The usage of &fsetup was introduced by MIT and PPPL
c     because of problems with Pathscale compiler 
c(cannot read more than five letters in the name of namelist group)
c..................................................................

      character*128  command
      character*8 line
      character*9 line0
      parameter (long_enough=100000)  !To accomodate namelist sections 
                                       !written on one line
      character(len=long_enough) :: line1   !Automatic array
      character(len=long_enough+1) :: line2   !Automatic array

c.......................................................................

CMPIINSERT_IF_RANK_NE_0_RETURN

      if (kopt.eq.0) then

         inlmod=0
         open(unit=4,file="cqlinput",delim='apostrophe',status="old") 
 1       read(unit=4,fmt=1003,end=2) line
         if (line.eq." &fsetup" .or. line.eq." &FSETUP" .or. 
     +       line.eq."&fsetup"  .or. line.eq."&FSETUP"  ) then
             inlmod=inlmod+1
           WRITE(*,*)'FSETUP namelist is present.'
           WRITE(*,*)'cqlinput is adjusted to change FSETUP to SETUP0.'
           WRITE(*,*)'After reading from the adjusted file,'
           WRITE(*,*)'  it will be restored back.'
           if (inlmod.eq.1) go to 2
         endif
         go to 1
 2       continue
         
         close(4)
         
 1003    format(a8)
         
c.......................................................................
c     Make copy if inlmod=1, and adjust nl structure
c.......................................................................
         
c     Copy cqlinput to cqlinput_tmp

         if (inlmod.le.0) then
cBH180613           WRITE(*,*)'FSETUP namelist is absent (probably using SETUP0)'
           goto 999 !Implies new namelist &setup0 is present (not $fsetup)
            
         else ! inlmod=1

            !Save the original file as cqlinput_tmp
            open(unit=4,file="cqlinput",delim='apostrophe',status="old")
            open(unit=5,file="cqlinput_tmp",delim='apostrophe',
     +           status="replace")
 3          read(unit=4,fmt='(a)',end=4) line1
            len_line1=len_trim(line1)
            if (len_line1.ge.(long_enough-1)) then
               WRITE(*,*)'len_line1,long_enough',len_line1,long_enough
               STOP 'Adjust long_enough'
            endif
            write(unit=5, fmt='(a)') trim(line1)
            go to 3
 4          continue
            close(4)
            close(5)
            
c     Reopen cqlinput as new file, and put adjusted namelist into it.
c     There are 2 &setup namelist sections.
            open(unit=4,file="cqlinput",delim='apostrophe',
     +           status="replace")
            open(unit=5,file="cqlinput_tmp",delim='apostrophe',
     +           status="old")
            ifirst=0
 5          read(unit=5,fmt='(a)',end=6) line1
            if ((line1(1:8).eq." &fsetup" .or. 
     +          line1(1:8).eq." &FSETUP"  .or.
     .          line1(1:7).eq."&fsetup"   .or. 
     +          line1(1:7).eq."&FSETUP")
     .           .and. (ifirst.eq.0)) then
               ifirst=1
               line2(1:8)=" &setup0"
               line2(9:len(line1)+1)=line1(9:len(line1))
            else
               line2=line1(1:len(line1))
            endif
            write(unit=4, fmt='(a)') trim(line2)
            go to 5
 6          continue
            close(4)
            close(5)
            
         endif ! inlmod
         
         go to 999
         
      elseif (kopt.eq.1) then
         
c.......................................................................
c     If cqlinput changed, copy cqlinput_tmp back to cqlinput
c.......................................................................
         if (inlmod.eq.1) then
            open(unit=5,file="cqlinput",delim='apostrophe',
     +           status="replace")
            open(unit=4,file="cqlinput_tmp",delim='apostrophe',
     +           status="old")
 7          read(unit=4,fmt='(a)',end=8) line1
            if (len_trim(line1).ge.(long_enough-1)) 
     .           STOP 'Adjust long_enough'
            write(unit=5, fmt='(a)') trim(line1)
            go to 7
 8          continue
            close(4)
            close(5)
         endif
         
      endif                     !on kopt
      
 999  return
      end

c

c==================================================================
      subroutine ain_transcribe(filename)
      implicit integer (i-n), real*8 (a-h,o-z)
CMPIINSERT_INCLUDE

c..................................................................
c     Transcibe filename contents to the standard output unit.
c..................................................................

      character(len=*), intent(in) :: filename
      parameter(long_enough=1000000)
      character(len=long_enough) :: line1          !Automatic array
      logical logic1

CMPIINSERT_IF_RANK_NE_0_RETURN

      WRITE(*,*)
      WRITE(*,*)  'ain_transcribe write of nl to stdout:'
      WRITE(*,*)  '[set nmlstout="disabled" to turn it off]'
      WRITE(*,*)

      max_length=0
      WRITE(*,*)'ain_transcribe: filename =',filename
      inquire(file=filename,iostat=kiostat,opened=logic1,number=inumber)
      WRITE(*,*)'ain_transcribe: inquire on ',filename,
     1   ' opened=',logic1,'iostat=',kiostat,'unit=',inumber
      open(unit=20,file=filename,delim='apostrophe',status="old",
     +     iostat=kiostat)
      if (kiostat.ne.0) WRITE(*,*)'ain_transcribe: kiostat=',kiostat
      if (kiostat.ne.0) STOP 'ain_transcribe: prblm with filename'
 3    read(unit=20,fmt='(a)',end=4) line1
      len_line1=len_trim(line1)
      if (len_line1.gt.max_length) max_length=len_line1
      WRITE(*,*) trim(line1)
      if (len_line1.ge.(long_enough-1)) then
         WRITE(*,*)'ain_transcribe,long_enough',len_line1,long_enough
         STOP 'Adjust long_enough'
      endif
      go to 3
 4    continue
      WRITE(*,*)'ain_transcribe:  max line length =',max_length
      close(20)
      WRITE(*,*)
 
      return
      end

      
