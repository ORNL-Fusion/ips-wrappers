      subroutine sfill (n, sa, sx, incx)
      use b2mod_types
      implicit none
      integer n, incx
      real (kind=R8) ::
     *   sa, sx(0:*)
*     ------------------------------------------------------------------
*     SFILL (BLAS-like routine) fills a vector with a constant value.
*     sfill(n,sa,sx,incx) sets all elements of the s.p. real*8 vector x
*     to sa. x is described by (n,sx,incx). BLAS storage conventions.
*     Return without effects if n.lt.1.
*     Undefined if incx.eq.0.
*     ------------------------------------------------------------------
      integer i, ix
      intrinsic max
      integer ibase
      ibase(n,incx) = max(0,(1-n)*incx)
*     ------------------------------------------------------------------
      if (incx.eq.0) then
       stop 'sfill--error, incx.eq.0'
      else if (incx.eq.1) then
       do i = 0, n-1
        sx(i) = sa
       enddo
      else
       ix = ibase(n,incx)
       do i = 0, n-1
        sx(ix) = sa
        ix = ix+incx
       enddo
      endif
      return
*     ------------------------------------------------------------------
      end
