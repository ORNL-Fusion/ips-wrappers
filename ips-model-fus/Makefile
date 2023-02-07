VPATH = $(PWD)
OBJDIR= .

include Makefile.config.perlmutter

F90_SRC = $(shell ls *.f90)

OBJS = $(F90_SRC:%.f90=$(OBJDIR)/%.o) $(F95_SRC:%.f95=$(OBJDIR)/%.o)

F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(NETCDFHOME)/include

EXECUTABLES = model_FUS
SCRIPTS = model_FUS.py  py_model_FUS.py

TARGETS =  $(EXECUTABLES)

all: .depend 
	- $(MAKE) $(TARGETS)
	@ for i in `echo $(EXECUTABLES)`; do \
	   if [ -f $$i ] ; then \
                echo "SUCCESS building $$i" >> make-result;\
        else \
                echo "***FAILURE building $$i" >> make-result;\
        fi;\
	done;
model_FUS: model_FUS.o
	$(F90) -o $@  $@.o $(LIBS)

install:
	-$(INSTALL) $(EXECUTABLES) $(SCRIPTS) $(IPS_ROOT)/bin

clean:
	$(RM) *.a *.mod *.o $(EXECUTABLES)

distclean: clean
	$(RM) .depend
	
.depend: $(F90MOD_SRC) $(F90_SRC) $(F95_SRC)
	touch .depend
	$(MAKEDEPF90) -b . $(F90MOD_SRC) $(F90_SRC) $(F95_SRC) > $@

-include .depend
