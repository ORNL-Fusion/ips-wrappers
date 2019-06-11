VPATH = $(PWD)
OBJDIR= .

include Makefile.include

F90_SRC = $(shell ls *.f90)

OBJS = $(F90_SRC:%.f90=$(OBJDIR)/%.o) $(F95_SRC:%.f95=$(OBJDIR)/%.o)
F90_MOD_INCLUDE =  $(F90_MOD_INCLUDE_PREFIX). $(F90_MOD_INCLUDE_PREFIX)$(IPS_ROOT)/include $(F90_MOD_INCLUDE_PREFIX)$(NETCDFHOME)/include

EXECUTABLES = model_RF_IC model_RF_IC_2 model_RF_IC_2_mcmd py_model_RF_IC_2_mcmd zero_RF_IC_power zero_RF_EC_power model_RF_IC_3 load_mdescr_sconfig
SCRIPTS = model_RF_IC.py model_RF_IC_2.py model_RF_IC_2_mcmd.py py_model_RF_IC_2_mcmd.py model_RF_IC_3.py

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
model_RF_IC: model_RF_IC.o
	$(F90) -o $@  $@.o $(LIBS)
model_RF_IC_2: model_RF_IC_2.o
	$(F90) -o $@  $@.o $(LIBS)
model_RF_IC_2_mcmd: model_RF_IC_2_mcmd.o
	$(F90) -o $@  $@.o $(LIBS)
model_RF_IC_3: model_RF_IC_3.o
	$(F90) -o $@  $@.o $(LIBS)
py_model_RF_IC_2_mcmd: py_model_RF_IC_2_mcmd.o
	$(F90) -o $@  $@.o
zero_RF_IC_power: zero_RF_IC_power.o
	$(F90) -o $@  $@.o $(LIBS)
zero_RF_EC_power: zero_RF_EC_power.o
	$(F90) -o $@  $@.o $(LIBS)
load_mdescr_sconfig: load_mdescr_sconfig.o
	$(F90) -o $@  $@.o $(LIBS)

install:
	$(INSTALL) $(EXECUTABLES) $(SCRIPTS) $(IPS_ROOT)/bin

clean:
	$(RM) $(OBJS) *.a *.mod $(TARGETS)

distclean:
	$(MAKE) clean
	$(RM) .depend

.depend: $(F90MOD_SRC) $(F90_SRC) $(F95_SRC)
#	touch .depend
	$(MAKEDEPF90) -b . $(F90MOD_SRC) $(F90_SRC) $(F95_SRC) > $@

-include .depend
