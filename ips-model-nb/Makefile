VPATH = $(PWD)
OBJDIR= .

include Makefile.include

F90_SRC = $(wildcard *.f90)
OBJS = $(F90_SRC:%.f90=%.o)

F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(IPS_ROOT)/include

EXECUTABLES = model_NB model_NB_2 model_NB_2_mcmd py_model_NB_2_mcmd 
SCRIPTS = model_NB.py model_NB_2.py model_NB_2_mcmd.py py_model_NB_2_mcmd.py 

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
model_NB: model_NB.o
	$(F90) -o $@  $@.o $(LIBS)
model_NB_2: model_NB_2.o
	$(F90) -o $@  $@.o $(LIBS)
model_NB_2_mcmd: model_NB_2_mcmd.o
	$(F90) -o $@  $@.o $(LIBS)
py_model_NB_2_mcmd: py_model_NB_2_mcmd.o
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
