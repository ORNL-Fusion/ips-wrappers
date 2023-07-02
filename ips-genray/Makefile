VPATH = $(PWD)
OBJDIR= .

include Makefile.config.perlmutter 

F90_SRC = $(shell ls *.f90)
F_SRC = $(shell ls *.f)

OBJS = $(F90_SRC:%.f90=$(OBJDIR)/%.o) $(F_SRC:%.f=$(OBJDIR)/%.o)
F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(IPS_ROOT)/include $(F90_INCLUDE_PREFIX)./genray_includes $(F90_INCLUDE_PREFIX)$(NETCDFHOME)/include

#BH130927:EXECUTABLES = prepare_genray_input process_genray_output process_genray_output_mcmd 
EXECUTABLES = prepare_genray_input process_genray_output
SCRIPTS = rf_genray.py genray_xmhd.py rf_genray_LH.py rf_genray_EC.py rf_genray_EC_p.py rf_genray_EC_dbb.py

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

prepare_genray_input: prepare_genray_input.o read_write_genray_input_prep.o bcast.o
	$(F90) -o $@ $^ $(LIBS)

process_genray_output: process_genray_output.o bcast.o
	$(F90) -o $@ $^ $(LIBS)

#BH130927:process_genray_output_mcmd: process_genray_output_mcmd.o bcast.o
#BH130927:	$(F90) -o $@ $^ $(LIBS)

install:
	$(INSTALL) $(EXECUTABLES) $(SCRIPTS) $(IPS_ROOT)/bin

clean:
	$(RM) $(OBJS) *.a *.mod $(TARGETS)

distclean: clean 
	$(RM) .depend
	
.depend: $(F90MOD_SRC) $(F90_SRC)
	$(MAKEDEPF90) -b . $(F90MOD_SRC) $(F90_SRC) > $@

-include .depend
