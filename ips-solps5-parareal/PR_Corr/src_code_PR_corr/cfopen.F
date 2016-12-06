      subroutine cfopen (ioun, chname, chstat, chform)
      use b2mod_types
      implicit none
      integer ioun
      character chname*(*), chstat*(*), chform*(*)
*     ------------------------------------------------------------------
*     CFOPEN opens an external file. The file may be formatted or
*     unformatted.
*!!   This routine may contain system dependencies.
*     The routine performs an operation akin to
*       open (ioun,file=chname,status=chstat,form=chform) ,
*     but using edited arguments; see the code for details.
*     ------------------------------------------------------------------
      logical ltst
      integer ios
      character fform0*(*), fname*1686, fstat*4, fform*12
      parameter (fform0='formatted')
      intrinsic len
      logical streql
      external xertst, xerrab, chcase, streql
*     ------------------------------------------------------------------
*!!   Note 0: set fform0 to the desired mode for un*formatted output.
*!!   Note 1, on porting this routine to another operating system:
*     The issues are: lower case or upper case for filenames and for
*     the other character arguments; maximum length of filenames;
*     restrictions on special characters in filenames. Check the
*     calls to this routine to see which special characters require
*     consideration (perhaps only the period (.), maybe the slash (/)).
*     ------------------------------------------------------------------
*   ..construct fname
*     (max length is len(fname) characters; make lowercase)
      fname = chname
      call chcase (0, fname)
      call xertst (len(chname).le.len(fname),
     &  'cfopen--faulty argument chname, too long: '//fname)
*   ..construct fstat
      if (streql(chstat,'old')) then
       fstat = 'old'
      else if (streql(chstat,'new')) then
       fstat = 'new'
      else
       fstat = chstat
       call xerrab ('cfopen--faulty argument chstat: '//fstat)
      endif
*   ..construct fform
*     ('formatted' or 'unformatted')
      if (streql(chform,'formatted')) then
       fform = 'formatted'
      else if (streql(chform,'unformatted')) then
       fform = 'unformatted'
      else if (streql(chform,'un*formatted')) then
       fform = fform0
      else
       fform = chform
       call xerrab ('cfopen--faulty argument chform: '//fform)
      endif
*   ..test file existence
      inquire (file=fname,exist=ltst)
      if (streql(fstat,'old')) then
       call xertst (ltst, 'cfopen--supposed old file does not '//
     &   'exist: file='//fname)
      else
       call xertst (.not.ltst, 'cfopen--supposed new file already '//
     &   'exists: file='//fname)
      endif
*   ..open file
      open (ioun,iostat=ios,file=fname,status=fstat,form=fform)
      call xertst (ios.eq.0, 'cfopen--open failure: file='//fname//
     &  ', status='//fstat//', form='//fform)
      return
*     ------------------------------------------------------------------
      end
