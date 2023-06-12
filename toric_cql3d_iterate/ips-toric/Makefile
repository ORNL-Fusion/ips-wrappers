VPATH = $(PWD)
OBJDIR= .

include Makefile.config.perlmutter

F90_SRC = $(shell ls *.f90)

OBJS = $(F90_SRC:%.f90=$(OBJDIR)/%.o) $(F95_SRC:%.f95=$(OBJDIR)/%.o)
F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(NETCDFHOME)/include

EXECUTABLES = prepare_toric_input process_toric_output prepare_abj_input

SCRIPTS = rf_ic_toric.py

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

toric_utils_mod.o:
	$(F90) -c toric_utils_mod.f90

prepare_toric_input: prepare_toric_input.o
	$(F90) -o $@ $@.o $(LIBS)

process_toric_output: process_toric_output.o
	$(F90) -o $@ $@.o $(LIBS)

prepare_abj_input: prepare_abj_input.o
	$(F90) -o $@ $@.o $(LIBS)

install:
	$(INSTALL) $(EXECUTABLES) $(SCRIPTS) $(IPS_ROOT)/bin

clean:
	$(RM) $(OBJS) *.a *.mod $(TARGETS) make-result

distclean:
	$(MAKE) clean
	$(RM) .depend

.depend: $(F90MOD_SRC) $(F90_SRC)
	touch .depend

-include .depend
# DO NOT DELETE
