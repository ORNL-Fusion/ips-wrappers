SRC:=src_code_PR_conv

NETCDFDIR = ${ITMLIBDIR}/netcdf/4.3.0
NETCDFFDIR = ${ITMLIBDIR}/netcdf-fortran/4.2/amd64_intel_12

FC:= ifort

ThisMachine:=$(shell uname -n)

ifneq (,$(findstring edison,$(ThisMachine)))
ThisMachine:=edison
endif

include Makefile.$(ThisMachine)

INCLUDEFLAGS:= -I${NETCDFDIR} -I${NETCDFFDIR}
LIBS:= -L${NETCDFFDIR} -lnetcdff -L ${NETCDFDIR} -lnetcdf

PR_GS2_conv: ${SRC}/PR_GS2_conv.f90
	${FC} ${INCLUDEFLAGS} $^ -o $@ ${LIBS}

PR_GS2_conv_oldVer: ${SRC}/PR_GS2_conv_oldVer.f90
	${FC} ${INCLUDEFLAGS} $^ -o $@ ${LIBS} 

clean:
	rm *.o test_output     

test:
	./PR_GS2_conv ${SRC}/b2time.nc ${SRC}/b2time2.nc test_output    
