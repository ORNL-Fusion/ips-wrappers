VPATH = $(PWD)
OBJDIR= .

include Makefile.include

F90_SRC = $(wildcard *.f90)
OBJS = $(F90_SRC:%.f90=%.o)

F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(IPS_ROOT)/include

EXECUTABLES = model_FUS model_FUS_2_mcmd 
SCRIPTS = model_FUS.py  py_model_FUS.py model_FUS_2_mcmd.py

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

model_FUS_2_mcmd: model_FUS_2_mcmd.o
	$(F90) -o $@  $@.o $(LIBS)

install:
	-$(INSTALL) $(EXECUTABLES) $(SCRIPTS) $(IPS_ROOT)/bin

clean:
	$(RM) *.a *.mod *.o $(EXECUTABLES)

distclean: clean
	$(RM) .depend
	
.depend: $(F90MOD_SRC) $(F90_SRC)
	$(MAKEDEPF90) -b ../obj $(F90MOD_SRC) $(F90_SRC) > $@

-include .depend
