VPATH = $(PWD)
OBJDIR= .

include Makefile.include

F90_SRC = $(shell ls *.f90)

OBJS = $(F90_SRC:%.f90=$(OBJDIR)/%.o) $(F95_SRC:%.f95=$(OBJDIR)/%.o)
F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(IPS_ROOT)/include $(F90_MOD_INCLUDE_PREFIX)$(NETCDFHOME)/include

EXECUTABLES = do_torlh_init_abr prepare_torlh_input_abr

SCRIPTS = rf_lh_torlh.py

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

do_torlh_init_abr.o: toric_utils_mod.o

do_torlh_init_abr: do_torlh_init_abr.o toric_utils_mod.o torlh_namelists.h
	$(F90) -g -o $@ $@.o toric_utils_mod.o $(LIBS)

prepare_torlh_input_abr: prepare_torlh_input_abr.o
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
