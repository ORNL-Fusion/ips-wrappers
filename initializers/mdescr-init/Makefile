VPATH = $(PWD)
OBJDIR= .

include Makefile.include

F90_SRC = $(shell ls *.f90)

OBJS = $(F90_SRC:%.f90=$(OBJDIR)/%.o) 
F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(IPS_ROOT)/include

EXECUTABLES = mdescr_state_init
SCRIPTS = mdescr_state_init.py

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
	
#init_empty_state_1: init_empty_state_1.o
#	$(F90) -o $@ $@.o $(LIBS)
	
#init_empty_state: init_empty_state.o
#	$(F90) -o $@ $@.o $(LIBS)
	
minimal_state_init: minimal_state_init.o
	$(F90) -o $@ $@.o $(LIBS)

mdescr_state_init: mdescr_state_init.o
	$(F90) -o $@ $@.o $(LIBS)
	
install:
	$(INSTALL) $(SCRIPTS) $(EXECUTABLES) $(IPS_ROOT)/bin

clean:
	$(RM) $(OBJS) *.a *.mod $(TARGETS)

distclean: 
	$(MAKE) clean
	$(RM) .depend
	
.depend: $(F90_SRC) $(F95_SRC)
	$(MAKEDEPF90) -b $(OBJDIR)  $(F90_SRC) $(F95_SRC) > $@

-include .depend
