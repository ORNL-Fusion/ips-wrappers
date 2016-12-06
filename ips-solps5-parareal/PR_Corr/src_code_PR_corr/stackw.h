*     (/stackw/ contains two workspace stacks.)
      integer nstckw, nstcwi
      parameter (nstckw=10, nstcwi=10)
      real (kind=R8) ::
     *   w(0:nstckw-1)
      integer wi(0:nstcwi-1)
      common /stackw/ w, wi
