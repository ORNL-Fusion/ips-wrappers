      include 'mpif.h'

      integer mpirank, mpisize, mpiierr, mpitag, mpisrc
      integer mpicmd, mpimode, mpiworker, mpill, mpil_
      integer mpistatus(MPI_STATUS_SIZE)

      real*8 mpitime, mpitime1

      common /mpi/ mpirank, mpisrc, mpisize, mpiierr, mpistatus, 
     +  mpicmd, mpimode, mpiworker, mpill,
     +  mpitime, mpitime1
