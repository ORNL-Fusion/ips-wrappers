      subroutine prvrti (n, kx, incx)
      use b2mod_types
      implicit none
      integer n, kx(0:*), incx
*     ------------------------------------------------------------------
*     PRVRTI sets all elements of an integer array to a value that is
*     expected to cause an error when used in an arithmetic operation.
*     This is used to protect against certain programming errors.
*!!   This routine may contain system dependencies.
*     The array is represented by the triplet (n,kx,incx), using BLAS
*     storage conventions.
*     The present implementation of prvrti operates by passing the array
*     kx to the routine ifill, which sets all elements to a constant
*     value, treated as integer in ifill. However, the constant value
*     passed to ifill is the real*8 value -1.0.
*     This does not actually cause the elements of kx to acquire an
*     illegal value on all computer systems--the result may just be
*     some large negative integer.
*!!   I do not know a machine-independent way to achieve the desired
*     result.
*     ------------------------------------------------------------------
      real (kind=R8) ::
     *   rnan
      integer inan(0:1)
      equivalence (rnan,inan)
      external ifill
*     ------------------------------------------------------------------
      rnan = -1.0_R8
      call ifill (n, inan(1), kx, incx)
*     (The second argument to ifill is deliberately of wrong type.)
*     ------------------------------------------------------------------
*     Note: When porting this routine to another operating system,
*     check if the IEEE not-a-number is available; if so, consider to
*     use it as the second argument to ifill.
*     ------------------------------------------------------------------
      end
