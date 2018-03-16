VPATH = $(PWD)
OBJDIR= .

include Makefile.include

F90_SRC = $(shell ls *.f90)

OBJS = $(F90_SRC:%.f90=$(OBJDIR)/%.o) $(F95_SRC:%.f95=$(OBJDIR)/%.o) $(F77_SRC:%.f=$(OBJDIR)/%.o)
F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(IPS_ROOT)/include

LIBS = $(PLASMA_STATE_LIB) $(NTCCLIB) $(MDSPLUSLIB) $(NETCDFLIB) $(LAPACKLIB) -L$(IPS_ROOT)/lib #-lswim-utils

PREPARE_INPUT_SRCS=
PREPARE_INPUT_OBJS = $(PREPARE_INPUT_SRCS:.f=.o)

EXECUTABLES = set_ps_beam_power
SCRIPTS =   nb_nubeam_pll.py nb_nubeam_mcmd.py nb_nubeam_dbb.py nb_nubeam.py

TARGETS = $(EXECUTABLES) 

all: .depend 
	$(RM) $(EXECUTABLES)
	- $(MAKE) $(TARGETS)
	@ for i in `echo $(EXECUTABLES)`; do \
	   if [ -f $$i ]; then \
                echo "SUCCESS building $$i" >> make-result;\
        else \
                echo "***FAILURE building $$i" >> make-result;\
        fi;\
	done;

set_ps_beam_power: set_ps_beam_power.o
	$(F90) -o $@  $@.o $(LIBS) 

install:
	$(INSTALL) $(EXECUTABLES) $(SCRIPTS) $(IPS_ROOT)/bin

clean:
	$(RM) $(OBJS) *.a *.mod $(EXECUTABLES)

distclean: clean 
	$(RM) .depend
	
.depend: $(F90MOD_SRC) $(F90_SRC)
	$(MAKEDEPF90) -b . $(F90MOD_SRC) $(F90_SRC) > $@

-include .depend
