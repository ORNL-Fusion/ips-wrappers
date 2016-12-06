      subroutine subsys (lbot, srname, srnam1)
      use b2mod_types
      implicit none
      logical lbot
      character srname*(*), srnam1*8
*     ------------------------------------------------------------------
*     SUBSYS is a group of routines (entry points in SUBSYS) that
*     perform certain functions related to workspace management and
*     error handling. The entry points are subini, subend, and subxnd.
*     SUBINI performs any tasks required at the start of subprogram
*     execution.
*     SUBEND performs any tasks required at the end of subprogram
*     execution.
*     SUBXND is employed to build down the stack when the program
*     is terminated because of an error condition.
*     ------------------------------------------------------------------
      integer nchmx, nchs, lcr0, lci0
      parameter (nchmx=100)
      character chstck(0:nchmx-1)*8
      integer lcstck(0:1,0:nchmx-1)
      save nchs, chstck, lcstck
      external lwget, lwiget, lwset, lwiset
      data nchs/0/
*     ------------------------------------------------------------------
*!!   Note: nchmx may be changed. nchmx=100 is more than enough to
*     handle any nonrecursive call chain that one might meet.
*     ------------------------------------------------------------------
      return
*     ------------------------------------------------------------------
      entry subini (srname)
#ifdef PERFMON
      call perfon (srname)
#endif
#ifdef _OPENMP
      return !!FR do not do the stack tracing if running in parallel
#else
*   ..test nchs
      if (nchs.eq.nchmx) then
       stop 'subini--error, nchs.eq.nchmx'
      endif
*   ..get base addresses for free workspace
      call lwget (0, lcr0)
      call lwiget (0, lci0)
*   ..push data onto the stack
      chstck(nchs) = srname
      lcstck(0,nchs) = lcr0
      lcstck(1,nchs) = lci0
      nchs = nchs+1
c-dpc      write(*,*) 'subini: ',srname,nchs
#endif
      return
*     ------------------------------------------------------------------
      entry subend ()
#ifdef _OPENMP
      return  !!FR do not do the stack tracing if running in parallel
#else
*   ..test nchs
      if (nchs.eq.0) then
       stop 'subend--error, nchs.eq.0'
      endif
*   ..pop data from the stack
*     (lcr0 and lci0 are the base addresses at invocation of the present
*     subprogram)
c-dpc      write(*,*) 'subend: ',chstck(nchs)
      nchs = nchs-1
      lcr0 = lcstck(0,nchs)
      lci0 = lcstck(1,nchs)
*   ..reset base addresses for free workspace
      call lwset (lcr0)
      call lwiset (lci0)
#ifdef PERFMON
      call perfoff
#endif
#endif

      return
*     ------------------------------------------------------------------
      entry subxnd (lbot, srnam1)
*   ..build down the stack
      if (nchs.eq.0) then
       lbot = .true.
       srnam1 = ' '
      else
       nchs = nchs-1
       lbot = .false.
       srnam1 = chstck(nchs)
       lcr0 = lcstck(0,nchs)
       lci0 = lcstck(1,nchs)
       call lwset (lcr0)
       call lwiset (lci0)
      endif
      return
*     ------------------------------------------------------------------
      end
