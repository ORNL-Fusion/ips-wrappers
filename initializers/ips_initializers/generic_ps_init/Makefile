VPATH = $(PWD)
OBJDIR= .

include Makefile.include

F90_SRC = $(wildcard *.f90)
OBJS = $(F90_SRC:%.f90=%.o)

F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(IPS_ROOT)/include

EXECUTABLES = generic_ps_init

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

generic_ps_init: generic_ps_init.o
	$(F90) -o $@ $@.o $(LIBS)

install:
	-$(INSTALL) $(EXECUTABLES) $(DESTDIR)/bin

clean:
	$(RM) *.a *.mod *.o $(EXECUTABLES)

distclean: clean
	$(RM) .depend

.depend: $(F90MOD_SRC) $(F90_SRC)
	$(MAKEDEPF90) -b ../obj $(F90MOD_SRC) $(F90_SRC) > $@

-include .depend
