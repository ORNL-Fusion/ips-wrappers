      subroutine prvrt (n, sx, incx)
      use b2mod_types
      implicit none
      integer n, incx
      real (kind=R8) ::
     *   sx(0:*)
*     ------------------------------------------------------------------
*     PRVRT sets all elements of a real*8 array to a value that is
*     expected to cause an error when used in an arithmetic operation.
*     This is used to protect against certain programming errors.
*!!   This routine may contain system dependencies.
*     The array is represented by the triplet (n,sx,incx), using BLAS
*     storage conventions.
*     The present implementation of prvrt operates by passing the array
*     sx to the routine sfill, which sets all elements to a constant
*     value, treated as real*8 in sfill. However, the constant value
*     passed to sfill is the integer -1.
*     At least on the NERSC Crays and on Sparc workstations, this
*     causes the elements of sx to acquire a value that is not a valid
*     real*8 number.
*     ------------------------------------------------------------------
      integer nan(0:1)
      real (kind=R8) ::
     *   rnan
      equivalence (nan,rnan)
      external sfill
*     ------------------------------------------------------------------
      nan(0) = -1
      nan(1) = 0
      call sfill (n, rnan, sx, incx)
*     (The second argument to sfill is deliberately of wrong type.)
*     ------------------------------------------------------------------
*     Note: When porting this routine to another operating system,
*     check if the IEEE not-a-number is available; if so, consider to
*     use it as the second argument to sfill.
*     ------------------------------------------------------------------
      end
