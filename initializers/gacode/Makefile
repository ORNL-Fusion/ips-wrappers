VPATH = $(PWD)
OBJDIR= .

include Makefile.config.perlmutter

F90_SRC = $(wildcard *.f90)
OBJS = $(F90_SRC:%.f90=%.o)

F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX).

EXECUTABLES = gacode_init
SCRIPTS = gacode_init.py

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

gacode_init: geo.o expro.o expro_comm.o expro_util.o gacode_init.o 
	$(F90) -o $@ $+ $(LIBS)

install:
	-$(INSTALL) $(EXECUTABLES) $(SCRIPTS) $(IPS_ROOT)/bin

clean:
	$(RM) *.a *.mod *.o $(EXECUTABLES)

distclean: clean
	$(RM) .depend

.depend: $(F90MOD_SRC) $(F90_SRC)
	$(MAKEDEPF90) -b ../obj $(F90MOD_SRC) $(F90_SRC) > $@

-include .depend
