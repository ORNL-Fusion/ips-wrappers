#! /bin/tcsh -f

echo "We have now switched to the cell-centred version of the code"
echo "You should probably do a"
echo " 'sbb ; gmake clean ; gmake listobj ; gmake depend ; gmake'"
echo "If you experience problems, the original version is in"
echo " 'src_staggered'"
echo "one level up from sbb"
echo "New: use your web browser to look at http://www.rzg.mpg.de/~dpc/solps.html"

setenv SOLPSTOP $PWD

if (-e whereami) then
  set iamat=`./whereami|tail -1`
  echo Running at $iamat
else
  set iamat="unknown"
endif

if (-e setup.csh.SOLPSMASTER) source setup.csh.SOLPSMASTER
if (-e setup.csh.SOLPSMASTER.local) source setup.csh.SOLPSMASTER.local

if(! $?SOLPSMASTER) then 
  switch ($iamat)
  case "IPP":
  case "EU_UNKNOWN":
    setenv SOLPSMASTER /afs/ipp-garching.mpg.de/u/dpc
    breaksw
  case "PPPL":
  case "USA_UNKNOWN":
    setenv SOLPSMASTER /afs/pppl.gov/u/dcoster
    breaksw
  default:
    setenv SOLPSMASTER /afs/ipp-garching.mpg.de/u/dpc
    breaksw
  endsw
endif

if($1 == "") then
  if (-e setup.csh.OBJECTCODE) source setup.csh.OBJECTCODE
  if (-e setup.csh.OBJECTCODE.local) source setup.csh.OBJECTCODE.local
  if(! $?OBJECTCODE) setenv OBJECTCODE `$SOLPSMASTER/@sys/machine`
else
  setenv OBJECTCODE $1
endif
if(! $?OBJECTCODE) then
  echo "OBJECTCODE not defined !"
endif

if (-e setup.csh.NAG) source setup.csh.NAG
if (-e setup.csh.NAG.$OBJECTCODE) source setup.csh.NAG.$OBJECTCODE
if (-e setup.csh.NAG.local) source setup.csh.NAG.local
if (-e setup.csh.NAG.local.$OBJECTCODE) source setup.csh.NAG.local.$OBJECTCODE

if(! $?NAG) then 
  source setup.csh.NAG.guess
endif

if (-e setup.csh.NCARG) source setup.csh.NCARG
if (-e setup.csh.NCARG.$OBJECTCODE) source setup.csh.NCARG.$OBJECTCODE
if (-e setup.csh.NCARG.local) source setup.csh.NCARG.local
if (-e setup.csh.NCARG.local.$OBJECTCODE) source setup.csh.NCARG.local.$OBJECTCODE

if(! $?NCARG_ROOT) then
  if (-e $SOLPSTOP/src/NCARG/$OBJECTCODE) setenv NCARG_ROOT $SOLPSTOP/src/NCARG/$OBJECTCODE
endif

if(! $?NCARG_ROOT) then 
  source setup.csh.NCARG.guess
endif
if(! $?NCARG) then 
  if($?NCARG_ROOT) then
    setenv NCARG_PATH `echo $NCARG_ROOT | sed -e "s:$SOLPSTOP/::"`
    switch ($NCARG_PATH)
    case '*3.*':
      echo 'Found NCAR Version 3.*'
      setenv NCARG '-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lncarg_loc -lX11 -lm'
      setenv NCAR_VERSION 3
      breaksw
    case '*4.*':
      echo 'Found NCAR Version 4.*'
      setenv NCARG '-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm'
      setenv NCAR_VERSION 4
      if (! $?SOLPS_CPP) then
        setenv SOLPS_CPP "-DNCAR4"
      else
        setenv SOLPS_CPP "$SOLPS_CPP -DNCAR4"
      endif
      breaksw
    default:
      echo 'Found NCAR Version ?.?; Assume 4'
      setenv NCARG '-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm'
      setenv NCAR_VERSION 4
      if (! $?SOLPS_CPP) then
        setenv SOLPS_CPP "-DNCAR4"
      else
        setenv SOLPS_CPP "$SOLPS_CPP -DNCAR4"
      endif
    endsw
  endif
endif

if (! $?GRAPHCAP) setenv GRAPHCAP X11

switch ($OBJECTCODE)
case "IBMaix":
  setenv nOBJECTCODE Aix
  breaksw
case "DECalpha":
  setenv nOBJECTCODE Alpha
  breaksw
case "SGIirix":
  setenv nOBJECTCODE Iris
  breaksw
case "linux.ifort64":
  setenv nOBJECTCODE Intel
  breaksw
case "sun4c":
  setenv nOBJECTCODE SunOS
  breaksw
case "sun5":
  setenv nOBJECTCODE Solaris
  breaksw
case "unicos":
  setenv nOBJECTCODE Unicos
  breaksw
default:
  setenv nOBJECTCODE Unknown
  breaksw
endsw

if (! $?DEVICE) setenv DEVICE iter

setenv GLI_HOME $SOLPSTOP/src/lib
# setenv GLI_WSTYPE 210
setenv GRSOFT_DEVICE "211 62"
setenv WSTYPE ${OBJECTCODE}
setenv SonnetTopDirectory ${SOLPSTOP}/src/Sonnet
setenv EscapeSonnet `echo ${SonnetTopDirectory} | sed 's:\/:\\\/:g'`

setenv DG ${SOLPSTOP}/src/DivGeo/dg
setenv CARRE_STOREDIR $SOLPSTOP/src/Carre/SAVE 

alias sb2 'cd ${SOLPSTOP}/src/Braams/b2'
alias sbb 'cd ${SOLPSTOP}/src/Braams/b2/base'
alias sbr 'cd ${SOLPSTOP}/src/Braams/b2/runs'
alias sei 'cd ${SOLPSTOP}/src/Eirene'
alias ssw 'cd ${SOLPSTOP}/src/Sonnet'
alias ssd 'cd ${SOLPSTOP}/src/DivGeo'
alias ssc 'cd ${SOLPSTOP}/src/Carre'
alias sbin 'cd ${SOLPSTOP}/src/bin/${OBJECTCODE}'
alias sbinc 'cd ${SOLPSTOP}/src/bin/common'
alias slib 'cd ${SOLPSTOP}/src/lib/${OBJECTCODE}'
alias srun 'cd ${SOLPSTOP}/src/Braams/b2/runs'
alias stop 'cd ${SOLPSTOP}'

alias sdg 'cd ${SOLPSTOP}/src/DivGeo/dg/class/${DEVICE}'
alias ssf 'cd ${SOLPSTOP}/src/Sonnet/device/${DEVICE}'
alias ssu 'cd ${SOLPSTOP}/src/DivGeo/uinp'

alias xyplot plot xyplot
alias xyplot2 plot xyplot2
alias xyplot3 plot xyplot3
alias xyplot4 plot xyplot4
alias xyplot5 plot xyplot5
alias xyplot6 plot xyplot6
alias xyplot7 plot xyplot7
alias xyplot8 plot xyplot8
alias xyplot8 plot xyplot8
alias xyplot9 plot xyplot9
alias xlyplot plot xlyplot
alias xlyplot2 plot xlyplot2
alias xlyplot3 plot xlyplot3
alias xlyplot4 plot xlyplot4
alias xlyplot5 plot xlyplot5
alias xlyplot6 plot xlyplot6
alias xlyplot7 plot xlyplot7
alias xlyplot8 plot xlyplot8
alias xlyplot8 plot xlyplot8
alias xlyplot9 plot xlyplot9
alias xylplot plot xylplot
alias xylplot2 plot xylplot2
alias xylplot3 plot xylplot3
alias xylplot4 plot xylplot4
alias xylplot5 plot xylplot5
alias xylplot6 plot xylplot6
alias xylplot7 plot xylplot7
alias xylplot8 plot xylplot8
alias xylplot8 plot xylplot8
alias xylplot9 plot xylplot9
alias xlylplot plot xlylplot
alias xlylplot2 plot xlylplot2
alias xlylplot3 plot xlylplot3
alias xlylplot4 plot xlylplot4
alias xlylplot5 plot xlylplot5
alias xlylplot6 plot xlylplot6
alias xlylplot7 plot xlylplot7
alias xlylplot8 plot xlylplot8
alias xlylplot8 plot xlylplot8
alias xlylplot9 plot xlylplot9

alias set_debug 'source $SOLPSTOP/debug'
alias unset_debug 'source $SOLPSTOP/nodebug'

if (-e setup.csh.$OBJECTCODE) source setup.csh.$OBJECTCODE
if (-e setup.csh.local) source setup.csh.local
if (-e setup.csh.local.$OBJECTCODE) source setup.csh.local.$OBJECTCODE
if (-e setup.csh.mdsplus) source setup.csh.mdsplus
if (-e setup.csh.mdsplus.$OBJECTCODE) source setup.csh.mdsplus.$OBJECTCODE
if (-e setup.csh.EIRENE) source setup.csh.EIRENE
if (-e setup.csh.EIRENE.$OBJECTCODE) source setup.csh.EIRENE.$OBJECTCODE
if (-e setup.csh.$USER) source setup.csh.$USER
if (-e setup.csh.$USER.$OBJECTCODE) source setup.csh.$USER.$OBJECTCODE

if (! $?OLDPATH) then
  setenv OLDPATH ${PATH}
endif
setenv PATH ${SOLPSTOP}/src/bin.local/${OBJECTCODE}:${SOLPSTOP}/src/bin.local/common:${SOLPSTOP}/src/bin/${OBJECTCODE}:${SOLPSTOP}/src/bin/common:${SOLPSTOP}/src/lib/python:${SonnetTopDirectory}/bin/${nOBJECTCODE}:${SOLPSTOP}/src/Braams/b2/base/${OBJECTCODE}:${PATH}
if ($?MANPATH) then
  if (! $?OLDMPATH) then
    setenv OLDMPATH ${MANPATH}
  else
    setenv MANPATH ${OLDMPATH}
  endif
  setenv MANPATH ${MANPATH}:${SonnetTopDirectory}/man:${DG}/equtrn/doxygen/man
else
  setenv MANPATH ${SonnetTopDirectory}/man:${DG}/equtrn/doxygen/man
endif

setenv SOLPS_LIB ${SOLPSTOP}/src/lib/${OBJECTCODE}

if ($?LD_LIBRARY_PATH) then
  if (! $?OLDLDPATH) then
    setenv OLDLDPATH ${LD_LIBRARY_PATH}
  else
    setenv LD_LIBRARY_PATH ${OLDLDPATH}
  endif
  setenv LD_LIBRARY_PATH ${SOLPS_LIB}:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH ${SOLPS_LIB}
endif

setenv PATH $NCARG_ROOT/bin:$PATH
setenv MANPATH $NCARG_ROOT/man:$MANPATH

if (! $?IDL_PATH) setenv IDL_PATH
if (! $?OLDIPATH) then
  setenv OLDIPATH ${IDL_PATH}
else
  setenv IDL_PATH ${OLDIPATH}
endif
setenv IDL_PATH +$SOLPSTOP/data/IDL:${IDL_PATH}

if ($?PYTHONPATH) then
  setenv PYTHONPATH ${PYTHONPATH}:${SOLPSTOP}/src/lib/python
else
  setenv PYTHONPATH ${SOLPSTOP}/src/lib/python
endif

setenv PLOT_SET_PATH '..:../..:${SOLPSTOP}/data/plot_set'

