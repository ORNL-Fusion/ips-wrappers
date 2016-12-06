      subroutine lwmain ()
      use b2mod_types
      implicit none
      integer lw0, lw1, nw
#     include "stackw.h"
*     ------------------------------------------------------------------
*     LWMAIN is a group of routines (entry points in lwmain) that
*     maintain a base address for free space in the workspace area
*     stackw%w. The entry points are lwset, lwget, and lwmax.
*     LWSET frees workspace by resetting the base address.
*     LWGET reserves workspace by incrementing the base address.
*     LWMAX returns the maximum workspace use.
*     ------------------------------------------------------------------
      integer lwrk, lwrkmx
      save lwrk, lwrkmx
      intrinsic max
      external prvrt
      data lwrk/0/, lwrkmx/0/
*     Invariant: 0.le.lwrk.le.lwrkmx.le.nstckw
*     ------------------------------------------------------------------
      return
*     ------------------------------------------------------------------
      entry lwset (lw0)
*     (revert base address to lw0)
      if (lw0.lt.0.or.lwrk.lt.lw0) then
       stop 'lwset--faulty argument lw0'
      endif
*     (set newly-free space to 'not-a-number')
      call prvrt (lwrk-lw0, w(lw0), 1)
      lwrk = lw0
      return
*     ------------------------------------------------------------------
      entry lwget (nw, lw1)
      if (nw.lt.0) then
       stop 'lwget--faulty argument nw'
      endif
      lw1 = lwrk
      lwrk = lwrk+nw
      lwrkmx = max(lwrk,lwrkmx)
      if (nstckw.lt.lwrk) then
       stop 'lwget--insufficient workspace stackw%w(0:*)'
      endif
      return
*     ------------------------------------------------------------------
      entry lwmax (lw1)
      lw1 = lwrkmx
      return
*     ------------------------------------------------------------------
      end subroutine lwmain
