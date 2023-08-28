VPATH = $(PWD)
OBJDIR= .

include Makefile.config.perlmutter 

F90_SRC = $(shell ls *.f90)
F_SRC = $(shell ls *.f)

OBJS = $(F90_SRC:%.f90=$(OBJDIR)/%.o) $(F_SRC:%.f=$(OBJDIR)/%.o)
F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(IPS_ROOT)/include $(F90_INCLUDE_PREFIX)$(NETCDFHOME)/include

EXECUTABLES = prepare_cql3d_input process_cql3d_output
SCRIPTS = fp_cql3d_general.py

TARGETS = $(EXECUTABLES) 

all: .depend 
	- $(MAKE) $(TARGETS)
	@ for i in `echo $(EXECUTABLES)`; do \
	   if [ -f $$i ]; then \
                echo "SUCCESS building $$i" >> $(IPS_ROOT)/.make-result;\
        else \
                echo "***FAILURE building $$i" >> $(IPS_ROOT)/.make-result;\
        fi;\
	done;

prepare_cql3d_input: prepare_cql3d_input.o aindflt.o  aindfpa.o  eqindflt.o  frinitl.o urfindfl.o  zcunix.o
	$(F90) -o $@ $^ $(LIBS)

process_cql3d_output: process_cql3d_output.o
	$(F90) -o $@ $^ $(LIBS)

install:
	$(INSTALL) $(EXECUTABLES) $(SCRIPTS) $(IPS_ROOT)/bin

clean:
	$(RM) $(OBJS) *.a *.mod $(TARGETS)

distclean: clean 
	$(RM) .depend

.depend: $(F90MOD_SRC) $(F90_SRC)
	$(MAKEDEPF90) -b . $(F90MOD_SRC) $(F90_SRC) > $@

-include .depend
