      subroutine xerset (nprt)
      use b2mod_types
      implicit none
      integer nprt, nprt1
*     ------------------------------------------------------------------
*     XERSET sets the unit number for error messages.
*     XERGET gets the unit number for error messages.
*     ------------------------------------------------------------------
      integer npdflt, npsave
      parameter (npdflt=0)
      save npsave
*     npdflt is the default unit number.
*     npsave is the designated unit number.
      data npsave/npdflt/
*     ------------------------------------------------------------------
*     This routine is used to specify a unit number for output from the
*     routine xerrab and from other error handling routines. If xerset
*     is called with a zero or positive argument nprt, then the error
*     handling routines will use unit number nprt.
*     If xerset is never called, or is called with a negative argument,
*     then a default unit number will be used.
*     A call to the alternate entry xerget returns the unit number to
*     be used for output by the error handling routines.
*     ------------------------------------------------------------------
      if (0.le.nprt) then
       npsave = nprt
      else
       npsave = npdflt
      endif
      return
*     ------------------------------------------------------------------
      entry xerget (nprt1)
      nprt1 = npsave
      return
*     ------------------------------------------------------------------
      end
