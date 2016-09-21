c
c
      subroutine eqindflt
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'name_decl.h'


c..................................................................
c     This routine set defaults for the "eq" module.
c..................................................................

      atol=1.e-8
      bsign=+1.0
      ellptcty=0.
      eqdskalt="disabled"
      eqdskin="eqdsk"
      eqmod="disabled"
      eqmodel="power"
      eqpower=1.
      eqsource="eqdsk"
      eqsym="average"
      fpsimodl="constant"
      lfield=lfielda
      methflag=10
      nconteq="psigrid"
      nconteqn=0
      povdelp=.2
      rbox=100.
      rboxdst=50.
      rmag=100.
      rtol=1.e-8
      zbox=100.
      return
      end
