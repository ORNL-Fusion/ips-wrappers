      subroutine lwimai ()
      use b2mod_types
      implicit none
      integer lwi0, lwi1, nwi
#     include "stackw.h"
*     ------------------------------------------------------------------
*     LWIMAI is a group of routines (entry points in lwimai) that
*     maintain a base address for free space in the workspace area
*     stackw%wi. The entry points are lwiset, lwiget, and lwimax.
*     LWISET frees workspace by resetting the base address.
*     LWIGET reserves workspace by incrementing the base address.
*     LWIMAX returns the maximum workspace use.
*     ------------------------------------------------------------------
      integer lwrki, lwrkimx
      save lwrki, lwrkimx
      intrinsic max
      external prvrti
      data lwrki/0/, lwrkimx/0/
*     Invariant: 0.le.lwrki.le.lwrkimx.le.nstcwi
*     ------------------------------------------------------------------
      return
*     ------------------------------------------------------------------
      entry lwiset (lwi0)
*     (revert base address to lwi0)
      if (lwi0.lt.0.or.lwrki.lt.lwi0) then
       stop 'lwiset--faulty argument lwi0'
      endif
*     (set newly-free space to 'not-a-number')
      call prvrti (lwrki-lwi0, wi(lwi0), 1)
      lwrki = lwi0
      return
*     ------------------------------------------------------------------
      entry lwiget (nwi, lwi1)
      if (nwi.lt.0) then
       stop 'lwiget--faulty argument nwi'
      endif
      lwi1 = lwrki
      lwrki = lwrki+nwi
      lwrkimx = max(lwrki,lwrkimx)
      if (nstcwi.lt.lwrki) then
       stop 'lwiget--insufficient workspace stackw%wi(0:*)'
      endif
      return
*     ------------------------------------------------------------------
      entry lwimax (lwi1)
      lwi1 = lwrkimx
      return
*     ------------------------------------------------------------------
      end subroutine lwimai
