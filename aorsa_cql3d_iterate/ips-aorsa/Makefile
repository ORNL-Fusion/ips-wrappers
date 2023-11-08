VPATH = $(PWD)
OBJDIR= .

include Makefile.config.perlmutter

F90_SRC = $(shell ls *.f90)

OBJS = $(F90_SRC:%.f90=$(OBJDIR)/%.o) $(F95_SRC:%.f95=$(OBJDIR)/%.o)
F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(NETCDFHOME)/include

EXECUTABLES = prepare_aorsa_input

SCRIPTS = rf_ic_aorsa.py

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

prepare_aorsa_input: aorsa2din_mod.o profile_mod.o  prepare_aorsa_input.o 
	$(F90) -o $@ $@.o aorsa2din_mod.o profile_mod.o $(LIBS) 

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
