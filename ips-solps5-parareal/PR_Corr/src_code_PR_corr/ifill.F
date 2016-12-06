      subroutine ifill (n, ka, kx, incx)
      use b2mod_types
      implicit none
      integer n, ka, kx(0:*), incx
*     ------------------------------------------------------------------
*     IFILL (BLAS-like routine) fills a vector with a constant value.
*     ifill(n,ka,kx,incx) sets all elements of the integer vector x
*     to ka. x is described by (n,kx,incx). BLAS storage conventions.
*     Return without effects if n.lt.1.
*     Undefined if incx.eq.0.
*     ------------------------------------------------------------------
      integer i, ix
      intrinsic max
      integer ibase
      ibase(n,incx) = max(0,(1-n)*incx)
*     ------------------------------------------------------------------
      if (incx.eq.0) then
       stop 'ifill--error, incx.eq.0'
      else if (incx.eq.1) then
       do i = 0, n-1
        kx(i) = ka
       enddo
      else
       ix = ibase(n,incx)
       do i = 0, n-1
        kx(ix) = ka
        ix = ix+incx
       enddo
      endif
      return
*     ------------------------------------------------------------------
      end
