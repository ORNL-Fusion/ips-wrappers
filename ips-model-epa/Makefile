VPATH = $(PWD)
OBJDIR= .

include Makefile.include

F90_SRC = $(shell ls *.f90)

OBJS = $(F90_SRC:%.f90=$(OBJDIR)/%.o) $(F95_SRC:%.f95=$(OBJDIR)/%.o)
F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(IPS_ROOT)/include $(F90_MOD_INCLUDE_PREFIX)$(NETCDFHOME)/include

#MDSPLUSLIB = -L$(MDSPLUS_ROOT)/lib -L$(GLOBUS_LOCATION)/lib -lMdsLib -lglobus_xio_gcc64
#LAPACKLIB= -L$(MKLHOME)/lib/64 -lmkl_lapack -lmkl_ipf -lguide -lpthread
#NETCDFLIB= -L$(NETCDFHOME)/lib -lnetcdf
#NTCCLIB = -L$(NTCCHOME)/lib -lxplasma2 -lgeqdsk_mds -lmdstransp -lvaxonly -lnscrunch -lfluxav -lr8bloat -lpspline -lezcdf -llsode -llsode_linpack -lcomput -lportlib
#PLASMA_STATE_LIB = -L$(IPS_ROOT)/lib -lplasma_state -lplasma_state_kernel -lps_xplasma2

#LIBS = $(PLASMA_STATE_LIB) $(NTCCLIB) $(MDSPLUSLIB) $(NETCDFLIB) $(LAPACKLIB) -L$(IPS_ROOT)/lib -lswim-utils

EXECUTABLES = model_epa model_epa_ps_file_init model_epa_ps_file_init_no_RF model_EPA_mdescr
SCRIPTS = model_epa.py model_epa_ps_init.py model_epa_ps_file_init.py  py_model_epa_ps_init.py model_epa_ps_file_init_no_RF.py model_EPA_mdescr.py

TARGETS = $(EXECUTABLES)

all: .depend
	- $(MAKE) $(TARGETS)
	@ for i in `echo $(EXECUTABLES)`; do \
	   if [ -f $$i ]; then \
                echo "SUCCESS building $$i" >> make-result;\
        else \
                echo "***FAILURE building $$i" >> make-result;\
        fi;\
	done;

model_epa:model_epa.o
	$(F90) -o $@ $@.o $(LIBS)

model_epa_ps_file_init:model_epa_ps_file_init.o
	$(F90) -o $@ $@.o $(LIBS)

model_epa_ps_file_init_no_RF:model_epa_ps_file_init_no_RF.o
	$(F90) -o $@ $@.o $(LIBS)

model_EPA_mdescr:model_EPA_mdescr.o
	$(F90) -o $@ $@.o $(LIBS)

install:
	$(INSTALL) $(EXECUTABLES) $(SCRIPTS) $(IPS_ROOT)/bin

clean:
	$(RM) $(OBJS) *.a *.mod $(TARGETS)

distclean:
	$(MAKE) clean
	$(RM) .depend

.depend: $(F90MOD_SRC) $(F90_SRC) $(F95_SRC)
	touch .depend
	$(MAKEDEPF90) -b . $(F90MOD_SRC) $(F90_SRC) $(F95_SRC) > $@

-include .depend
