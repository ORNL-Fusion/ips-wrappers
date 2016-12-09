# Building the Convergence and Corrector binaries
First the corrector ...
Default: corrects all species in SOLPS including neutrals.
To avoid correcting neutrals (such as in a case where coarse solver uses fluid neutrals while fine solver uses Eirene):
set CORRECT_NEUTRALS=0 in buildAndTestCorr.sh
```
cd ips-solps5-parareal
source env.edison.sh
cd PR_Corr
./buildAndTestCorr.sh
```
Which should give output as ...
```
 ... 
 as zamin is 0 at is=           8 no corr
 as zamin is 0 at is=           8 no corr
 VERSION: 01.001.029
  Normal end of execution.
 Time for run is  0.200011000000000
```
Then the convergence ...
```
cd ../PR_Conv
make clean
make
make test
```
Which should give output as ...
```
./conv_solps src_code_PR_conv/b2time.nc src_code_PR_conv/b2time2.nc test_output
 pwmxip id found
 var_namepwmxip             xtype           6 ndims           2
 IDs of dimensions of flux variable           9 &          11
 species & time_index lengths:nc                 :           2
 time               :           2
 Error =   2.045924170498220E-003 tolerance  0.000000000000000E+000
```
