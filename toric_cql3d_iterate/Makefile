VPATH = $(PWD)
OBJDIR= .

include Makefile.include

F90_SRC = $(shell ls *.f90)

OBJS = $(F90_SRC:%.f90=$(OBJDIR)/%.o) $(F95_SRC:%.f95=$(OBJDIR)/%.o)
F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(IPS_ROOT)/include $(F90_MOD_INCLUDE_PREFIX)$(NETCDFHOME)/include

EXECUTABLES = do_toric_init_abr prepare_toric_input_abr

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

do_toric_init_abr.o: toric_utils_mod.o

do_toric_init_abr: do_toric_init_abr.o toric_utils_mod.o toric_namelists.h
	$(F90) -g -traceback -o $@ $@.o toric_utils_mod.o $(LIBS)

prepare_toric_input_abr: prepare_toric_input_abr.o
	$(F90) -g -traceback -o $@ $@.o $(LIBS)

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
