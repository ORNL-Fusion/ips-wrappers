F90 = ftn
F90FLAGS = -x f95-cpp-input -g -fbounds-check
F90_FREEFORM_FLAG = -ffree-form
F90_FIXEDFORM_FLAG = -ffixed-form
F90_MOD_INCLUDE_PREFIX = -I
F90_INCLUDE_PREFIX = -I

F77 = ftn
F77FLAGS = -x f95-cpp-input
F77_FREEFORM_FLAG = -ffree-form
F77_FIXEDFORM_FLAG = -ffixed-form
F77_INCLUDE_PREFIX = -I

IPS_PHYS_ROOT=/project/projectdirs/m876/phys-bin/phys
NTCCHOME =/project/projectdirs/atom/atom-install-edison/plasma-state-source/ntcc-gnu
NETCDFHOME=$(NETCDF_DIR)

MDSPLUS_ROOT=$(NTCCHOME)
MKLHOME = /usr

LAPACKHOME=/usr
LAPACKLIB= -L$(LAPACKHOME)/lib64 -llapack -lblas
LAPACKLIB= 

NETCDFLIB= -L$(NETCDFHOME)/lib -lnetcdff -lnetcdf

NTCCLIB = -L$(NTCCHOME)/lib -lxplasma_debug -lxplasma2 -lgeqdsk_mds -lvaxonly -lnscrunch -lfluxav -lr8bloat -lpspline -lezcdf -llsode -llsode_linpack -lcomput -lportlib -lureadsub -lsmlib
PLASMA_STATE_LIB = -L$(IPS_ROOT)/lib -lplasma_state -lps_xplasma2 -lplasma_state_kernel

MDSPLUSLIB = -L$(MDSPLUS_ROOT)/lib -lmdstransp
LIBS = $(PLASMA_STATE_LIB) $(NTCCLIB) $(MDSPLUSLIB) $(NETCDFLIB) $(LAPACKLIB) -L$(IPS_ROOT)/lib -lswim-utils

INSTALL = /usr/bin/install
